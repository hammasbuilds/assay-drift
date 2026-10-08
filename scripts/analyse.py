"""Score every assay against every downloaded genome, grouped by collection date.

    python scripts/analyse.py [--data data/sequences.jsonl] [--json results.json]

Reads what `scripts/fetch.py` downloaded and answers, per assay and per quarter:
of the sequences collected then, what share would this test still have detected?

The output deliberately carries four numbers per period rather than one rate:

    n        sequences with this assay evaluated
    usable   ... and with no unknown bases under any oligo
    intact   ... and still expected to amplify
    excl%    share dropped for unknown bases

A single rate hides which of "the assay failed" and "the data could not say"
is responsible, and those call for opposite actions.
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift import primers  # noqa: E402
from assaydrift.analyze import evaluate, summarise, trend  # noqa: E402
from assaydrift.mutations import observe, places, rank  # noqa: E402
from assaydrift.ncbi import Record  # noqa: E402

try:
    from importlib.metadata import version
except ImportError:  # pragma: no cover - 3.8 and older
    version = None

ROOT = Path(__file__).resolve().parent.parent


def load(path: Path):
    if not path.exists():
        raise SystemExit(f"no data at {path} - run scripts/fetch.py first")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            yield Record(
                accession=row["accession"],
                organism=row["organism"],
                sequence=row["sequence"],
                collected=row["collected"],
                country=row.get("country", ""),
                length=row["length"],
            )


def baseline_for(assay) -> int:
    """The assay's own mismatch count against the 2019 reference.

    Charite RdRp shipped with one. Scoring it against zero would report it as
    broken for every sequence ever collected, including the reference itself.
    """
    return sum(primers.expected_reference_mismatches(o.name) for o in assay.oligos)


def _tool_version() -> str | None:
    """The installed version, or None when running from source.

    The README states that the demo and the tests need no install step, so this must not
    raise when the distribution metadata is absent.
    """
    if version is None:
        return None
    try:
        return version("assay-drift")
    except Exception:  # noqa: BLE001 - absence is the normal case here
        return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    ap.add_argument("--json", type=Path, default=ROOT / "results" / "drift.json")
    ap.add_argument("--granularity", default="quarter", choices=["quarter", "year", "month"])
    ap.add_argument(
        "--mutations-json",
        type=Path,
        default=ROOT / "results" / "mutations.json",
        help="where the named-mutation table (C28311T, ...) is written",
    )
    ap.add_argument(
        "--top-mutations",
        type=int,
        default=5,
        help="how many named mutations to keep per assay, most frequent first",
    )
    args = ap.parse_args()

    records = list(load(args.data))
    print(f"{len(records):,} sequences\n")

    per_assay = defaultdict(list)
    undated = 0
    for record in records:
        for assay in primers.ASSAYS:
            result = evaluate(record, assay, args.granularity)
            if result is None:
                undated += 1
                continue
            per_assay[assay.name].append(result)

    payload = {}
    for assay in primers.ASSAYS:
        results = per_assay[assay.name]
        if not results:
            continue
        base = baseline_for(assay)
        summaries = summarise(results, baseline=base)
        movement = trend(summaries)

        label = f"{assay.name}  ({assay.target_gene})"
        if assay.retired:
            label += "   [RETIRED BY CDC]"
        if base:
            label += f"   [baseline {base} mismatch vs reference]"
        print(label)
        print(
            f"  {'period':<9}{'n':>6}{'usable':>8}{'eff n':>7}{'perfect':>9}"
            f"{'failing':>9}{'excl%':>8}{'3prime':>8}{'lost':>6}"
        )
        for summary in summaries.values():
            pr, fr = summary.perfect_rate, summary.failing_rate
            print(
                f"  {summary.period:<9}{summary.total:>6}{summary.usable:>8}"
                f"{('    n/a' if summary.effective_n is None else f'{summary.effective_n:>7.1f}')}"
                f"{('      n/a' if pr is None else f'{pr:>8.1%}')}"
                f"{('      n/a' if fr is None else f'{fr:>8.1%}')}"
                f"{summary.excluded_rate:>8.1%}{summary.three_prime:>8}"
                f"{summary.lost_oligo:>6}" + ("" if summary.reliable else "   (too few)")
            )
        if movement["enough_to_say"]:
            print(
                f"  perfect match   {movement['first_period']} "
                f"{movement['first_rate']:.1%} (n={movement['first_n']})  ->  "
                f"{movement['last_period']} {movement['last_rate']:.1%} "
                f"(n={movement['last_n']})   change {movement['change']:+.1%}"
            )
            print(
                f"  likely failing  {movement['first_period']} "
                f"{movement['first_failing_rate']:.1%}  ->  "
                f"{movement['last_period']} {movement['last_failing_rate']:.1%}"
                f"   worst {movement['worst_failing_period']} "
                f"{movement['worst_failing_rate']:.1%}"
            )
        else:
            print(
                f"  not enough reliable periods to claim a trend "
                f"({movement['periods']} with >= 25 usable sequences)"
            )
        print()

        payload[assay.name] = {
            "target_gene": assay.target_gene,
            "source": assay.source,
            "retired": assay.retired,
            "baseline_mismatches": base,
            "trend": movement,
            "periods": {
                p: {
                    "total": s.total,
                    "usable": s.usable,
                    # How many INDEPENDENT observations the rates below rest on. `usable`
                    # counts sequences; consecutive GenBank accessions come from one
                    # submission, so a quarter can be 97% a single batch. See
                    # PeriodSummary.effective_n.
                    "clusters": s.clusters,
                    "effective_n": (None if s.effective_n is None else round(s.effective_n, 2)),
                    "perfect": s.perfect,
                    "failing": s.failing,
                    "perfect_rate": s.perfect_rate,
                    "failing_rate": s.failing_rate,
                    "excluded_rate": s.excluded_rate,
                    "three_prime": s.three_prime,
                    "lost_oligo": s.lost_oligo,
                    "median_mismatches": s.median_mismatches,
                    "reliable": s.reliable,
                }
                for p, s in summaries.items()
            },
        }

    # Which corpus produced these numbers. Without it the file is a set of rates with no
    # way to tell what they are rates over, and the corpus is gitignored.
    manifest = ROOT / "data" / "accessions.tsv"
    identity = ""
    manifest_records = 0
    if manifest.is_file():
        for line in manifest.read_text(encoding="utf-8").splitlines():
            if line.startswith("# identity sha256:"):
                identity = line.split(":", 1)[1].strip()
            elif line.startswith("# records:"):
                manifest_records = int(line.split(":", 1)[1].strip() or 0)
    payload["_run"] = {
        "generated": datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "granularity": args.granularity,
        "sequences_read": len(records),
        "manifest": manifest.relative_to(ROOT).as_posix() if manifest.is_file() else None,
        "manifest_records": manifest_records or None,
        "corpus_identity_sha256": identity or None,
        "tool_version": _tool_version(),
        "reproduce": ("python scripts/fetch.py --from-manifest && python scripts/analyse.py"),
        "note": (
            "Period rates carry `clusters` and `effective_n` as well as `total` and "
            "`usable`. The first two count independent submitter groups; the last two "
            "count sequences, and consecutive GenBank accessions come from one "
            "submission. A period can pass the 25-sequence floor on an effective n near "
            "1, so read effective_n before any rate."
        ),
    }

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"written to {args.json}")

    print("\nNamed mutations under each oligo (most frequent first)")
    observations, mut_denominators = observe(records, primers.ASSAYS, args.granularity)
    mutations_payload = {}
    for assay in primers.ASSAYS:
        assay_observations = {
            key: obs for key, obs in observations.items() if obs.site.assay == assay.name
        }
        ranked = rank(assay_observations, mut_denominators)[: args.top_mutations]
        if not ranked:
            continue
        print(f"  {assay.name}")
        entries = []
        for obs, overall_rate in ranked:
            window_flag = " [3' window]" if obs.site.in_three_prime_window else ""
            top_places = places(obs, mut_denominators)
            print(
                f"    {obs.name:<10} {obs.site.oligo:<10} base {obs.site.base_in_oligo} of "
                f"{obs.site.oligo_length}  overall {overall_rate:.1%}{window_flag}"
            )
            if obs.site.in_three_prime_window:
                # Where a 3'-window mutation was collected decides how much a
                # per-quarter rate can carry: a batch from one laboratory moves
                # a quarter a long way without telling you anything about what
                # is circulating anywhere else.
                shown = ", ".join(f"{place} {n}/{total}" for place, n, total, _rate in top_places)
                print(f"      collected in: {shown}")
            entries.append(
                {
                    "name": obs.name,
                    "oligo": obs.site.oligo,
                    "role": obs.site.role,
                    "base_in_oligo": obs.site.base_in_oligo,
                    "oligo_length": obs.site.oligo_length,
                    "bases_from_3_prime_end": obs.site.bases_from_3_prime_end,
                    "in_three_prime_window": obs.site.in_three_prime_window,
                    "genome_position": obs.site.genome_position,
                    "overall_rate": overall_rate,
                    "total": obs.total,
                    "by_period": dict(sorted(obs.by_period.items())),
                    "top_places": [
                        {"place": place, "carrying": n, "sequenced": total, "rate": rate}
                        for place, n, total, rate in top_places
                    ],
                }
            )
        mutations_payload[assay.name] = entries
    args.mutations_json.parent.mkdir(parents=True, exist_ok=True)
    args.mutations_json.write_text(json.dumps(mutations_payload, indent=1), encoding="utf-8")
    print(f"written to {args.mutations_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
