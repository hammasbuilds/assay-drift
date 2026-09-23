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
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift import primers  # noqa: E402
from assaydrift.analyze import evaluate, summarise, trend  # noqa: E402
from assaydrift.ncbi import Record  # noqa: E402

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


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    ap.add_argument("--json", type=Path, default=ROOT / "results" / "drift.json")
    ap.add_argument("--granularity", default="quarter", choices=["quarter", "year", "month"])
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
            f"  {'period':<9}{'n':>6}{'usable':>8}{'intact':>8}{'rate':>8}"
            f"{'excl%':>8}{'3prime':>8}{'lost':>6}"
        )
        for summary in summaries.values():
            rate = summary.intact_rate
            print(
                f"  {summary.period:<9}{summary.total:>6}{summary.usable:>8}"
                f"{summary.intact:>8}"
                f"{('  n/a' if rate is None else f'{rate:>7.1%}')}"
                f"{summary.excluded_rate:>8.1%}{summary.three_prime:>8}"
                f"{summary.lost_oligo:>6}" + ("" if summary.reliable else "   (too few)")
            )
        if movement["enough_to_say"]:
            print(
                f"  {movement['first_period']} {movement['first_rate']:.1%} "
                f"(n={movement['first_n']})  ->  "
                f"{movement['last_period']} {movement['last_rate']:.1%} "
                f"(n={movement['last_n']})   change {movement['change']:+.1%}"
                f"   worst {movement['worst_period']} {movement['worst_rate']:.1%}"
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
                    "intact": s.intact,
                    "intact_rate": s.intact_rate,
                    "excluded_rate": s.excluded_rate,
                    "three_prime": s.three_prime,
                    "lost_oligo": s.lost_oligo,
                    "median_mismatches": s.median_mismatches,
                    "reliable": s.reliable,
                }
                for p, s in summaries.items()
            },
        }

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"written to {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
