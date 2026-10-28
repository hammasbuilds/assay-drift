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


def default_manifest(data: Path) -> Path:
    """The manifest that goes with a corpus file.

    fetch.py writes data/sequences.jsonl beside data/accessions.tsv, and
    data/sequences_big.jsonl beside data/accessions_big.tsv, so the name follows the data
    rather than being fixed.
    """
    name = data.name.replace("sequences", "accessions")
    if name.endswith(".jsonl"):
        name = name[: -len(".jsonl")] + ".tsv"
    return data.parent / name


def show(path: Path) -> str:
    """A path as a reader should see it, for any path at all.

    `path.relative_to(ROOT)` raises for a relative path and for anything outside the
    repository, and this is called while building the metadata block - after the whole
    corpus has been read. A run over 291,827 genomes reached the last step and died on
    `--manifest data/accessions_big.tsv`, which is relative, so relative_to had no
    absolute path to subtract ROOT from. Hours of reading were discarded by a line whose
    only job was to make a filename look tidy. Nothing presentational may raise here.

    `absolute()`, not `resolve()`: data/ in this checkout is a Windows junction to a
    folder elsewhere on the drive, and resolve() follows it, so every path under data/
    resolves to somewhere outside the repository and relative_to fails for the one
    directory whose files this script exists to read.
    """
    try:
        return path.absolute().relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def baseline_for(assay) -> int:
    """The assay's own mismatch count against the 2019 reference.

    Charite RdRp shipped with one. Scoring it against zero would report it as
    broken for every sequence ever collected, including the reference itself.
    """
    return sum(primers.expected_reference_mismatches(o.name) for o in assay.oligos)


def _manifest_lines(manifest: Path) -> list[str]:
    """The manifest's lines, plain or gzipped - the same reader fetch.py uses.

    The large manifest is committed gzipped (13.4 MB of text, 1.54 MB compressed),
    because a results file published from a corpus the repository cannot rebuild is not
    reproducible, whatever its `reproduce` field says.
    """
    if manifest.suffix == ".gz":
        import gzip

        with gzip.open(manifest, "rt", encoding="utf-8") as handle:
            return handle.read().splitlines()
    return manifest.read_text(encoding="utf-8").splitlines()


def _source_revision() -> str | None:
    """The git revision of this checkout, or None outside one.

    `tool_version` is None whenever the package is not installed, which the README makes
    the normal way to run this - so a published results file recorded nothing at all
    about the code that produced it. The corpus is identified by a digest; the code
    should be identified too, or a rate cannot be attributed to a version of the
    analysis.

    Never raises: git may be absent, and this is metadata on a file that took hours.
    """
    try:
        import subprocess

        out = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        rev = out.stdout.strip()
        if out.returncode or not rev:
            return None
        dirty = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain"],
            capture_output=True, text=True, timeout=10,
        )
        return f"{rev}-dirty" if dirty.stdout.strip() else rev
    except Exception:  # noqa: BLE001 - absence is a normal case here
        return None


def _reproduce(args: argparse.Namespace, manifest: Path) -> str:
    """The command that rebuilds this file, with this run's own arguments in it.

    Paths go through show(), so a manifest outside the repository prints as an absolute
    path rather than raising - and data/ in this checkout is a Windows junction, which
    is why show() uses absolute() and not resolve().
    """
    parts = ["python scripts/fetch.py --from-manifest"]
    if manifest.is_file():
        parts[0] += f" --manifest {show(manifest)}"
    run = "python scripts/analyse.py"
    if args.data is not None:
        run += f" --data {show(args.data)}"
    if manifest.is_file():
        run += f" --manifest {show(manifest)}"
    if args.granularity:
        run += f" --granularity {args.granularity}"
    parts.append(run)
    return " && ".join(parts)


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
    ap.add_argument(
        "--manifest",
        type=Path,
        help="the manifest identifying this corpus (default: derived from --data, so "
        "sequences_big.jsonl pairs with accessions_big.tsv)",
    )
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
    # Resolved here, at the boundary. Downstream code compares these against ROOT and
    # writes them into the report; a relative path works for opening a file and then
    # fails on anything that reasons about where the file is.
    for name in ("data", "json", "manifest", "mutations_json"):
        value = getattr(args, name)
        if value is not None:
            # absolute(), not resolve(): see show() - data/ is a junction here.
            setattr(args, name, Path(value).absolute())

    # Two streaming passes, not one list. `load` is already a generator and the records
    # are needed twice - here and by observe() below - so the old `list(load(...))` held
    # the whole corpus in memory. That cost nothing at 2,765 sequences and does not
    # survive 291,827: the corpus is 8.15 GB on disk and far more as Python objects.
    # Re-reading the file is the cheaper resource.
    sequences_read = 0
    per_assay = defaultdict(list)
    undated = 0
    for record in load(args.data):
        sequences_read += 1
        # Streaming removed the one line that said how big the corpus was, and a run over
        # 291,827 sequences then printed nothing at all for its first hour. Progress every
        # 25,000 keeps a long run legible without flooding the log.
        if sequences_read % 25_000 == 0:
            print(f"  read {sequences_read:,} sequences", flush=True)
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
                f"(n={movement['last_n']}, eff {movement['last_effective_n']:.0f})"
                f"   change {movement['change']:+.1%}"
                + (
                    f"   [{len(movement['periods_after_the_endpoint'])} later period(s) "
                    f"too lopsided to end on]"
                    if movement["periods_after_the_endpoint"]
                    else ""
                )
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
    # Follows --data. Hardcoded to data/accessions.tsv, a run over sequences_big.jsonl
    # recorded manifest_records: 2,765 beside sequences_read: 291,827 - the field whose
    # whole purpose is to say which corpus produced these rates was naming a different one.
    manifest = args.manifest or default_manifest(args.data)
    identity = ""
    manifest_records = 0
    if manifest.is_file():
        for line in _manifest_lines(manifest):
            if line.startswith("# identity sha256:"):
                identity = line.split(":", 1)[1].strip()
            elif line.startswith("# records:"):
                manifest_records = int(line.split(":", 1)[1].strip() or 0)
    payload["_run"] = {
        "generated": datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "granularity": args.granularity,
        "sequences_read": sequences_read,
        "manifest": show(manifest) if manifest.is_file() else None,
        "manifest_records": manifest_records or None,
        # Stated rather than left for a reader to notice. The manifest identifies the
        # corpus, so one that counts a different number of records is not this corpus's.
        "manifest_matches_data": (
            None if not manifest_records else manifest_records >= sequences_read
        ),
        "corpus_identity_sha256": identity or None,
        "tool_version": _tool_version(),
        "source_revision": _source_revision(),
        # Built from this run's arguments. It was a constant, so a file produced from the
        # 291,827-sequence corpus carried the command that rebuilds the 2,765-sequence one
        # - the same defect as `manifest` naming data/accessions.tsv whatever was read,
        # one field further down, and the field a reader trusts most because it claims to
        # be executable.
        "reproduce": _reproduce(args, manifest),
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
    # stream=, not a list: observe() exposes it precisely so the corpus need not be held
    # in memory, and the previous caller passed a materialised list anyway.
    observations, mut_denominators = observe(
        assays=primers.ASSAYS, granularity=args.granularity, stream=load(args.data)
    )
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
