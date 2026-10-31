"""How much of the `likely_failing` rate is the rule rather than the data.

A 2025 wet-lab study (Frontiers in Cellular and Infection Microbiology, 10.3389/
fcimb.2025.1524025) tested assays against mismatched templates in vitro and found most
kept working - concluding that in-silico predictions of false negatives may **overstate**
the risk. `likely_failing` is exactly such a prediction, so the honest response is not to
defend it but to measure how much of it depends on choices nobody validated against a
bench.

The rule has three clauses (`analyze.AssayResult.likely_failing`):

  1. an oligo with no binding site left at all,
  2. any mismatch within `THREE_PRIME_WINDOW` bases of a primer's 3' end (default 5),
  3. `total_mismatches >= SEVERE_MISMATCH_LOAD` across the assay (default 3).

Clauses 2 and 3 carry a threshold each, and neither number came from a wet-lab
calibration - 5 is a primer-design rule of thumb and 3 is a judgement. So this sweeps
them:

  * 3' window: 3, 5, 7
  * mismatch load: 2, 3, 4, 5
  * and three rule shapes: all three clauses, clause 2 dropped, clause 1 alone

Clause 1 alone is the conservative floor - the only prediction a bench would not dispute,
because an oligo with no binding site cannot prime. The distance between that floor and
the published rate is the size of the claim being made on top of the data.

Run on `data/sequences.jsonl` (2,765 genomes, minutes). The published headline uses the
291,827-genome corpus at window 5 / load 3; one pass over that takes about three hours,
so the grid is measured on the small corpus and its direction reported, not extrapolated
to the large one.

    .venv-research/Scripts/python research/sensitivity.py
    .venv-research/Scripts/python research/sensitivity.py --data data/sequences_big.jsonl
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

WINDOWS = (3, 5, 7)
LOADS = (2, 3, 4, 5)
SHAPES = ("all three clauses", "clustered 3' clause", "without the 3' clause", "oligo lost only")
# `clustered 3' clause` is the only shape here with wet-lab evidence behind its
# threshold, and it was added after @wetlab2025 was read in full rather than quoted from
# a search summary. Otwell et al. tested, among others:
#
#   * a reverse-primer mismatch ~1 nt from the 3' end (BVP 501Y, Omicron) - about a
#     1-cycle Ct shift, no false negative;
#   * two forward-primer mismatches in the same assay - a 9 to 12 cycle shift and a
#     genuine false negative;
#   * three mismatches within the last four 3' nucleotides (C3 ORF3a) - failure below
#     5,000 copies, Ct shifts above 15.
#
# A single 3'-proximal mismatch is therefore tolerated and a cluster is not, which is
# precisely what the published clause - `three_prime_mismatches > 0` - cannot express.
# This shape requires two.
CLUSTERED_MINIMUM = 2


def wilson(successes: int, n: int, z: float = 1.96) -> dict | None:
    if n == 0:
        return None
    p = successes / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "rate": round(p, 4),
        "lo": round(max(0.0, centre - half), 4),
        "hi": round(min(1.0, centre + half), 4),
        "n": n,
        "failing": successes,
    }


def features(data: Path, window: int) -> dict[str, list[tuple[bool, int, int, bool]]]:
    """Per sequence per assay: (any_oligo_lost, 3' mismatches, total mismatches, usable).

    Aligned once per window, because the 3' count is computed during alignment against
    `match.THREE_PRIME_WINDOW` and `Hit` keeps the count, not the positions - so the
    window cannot be varied after the fact. The mismatch load can, which is why only the
    three window values cost a pass each.
    """
    import assaydrift.match as match

    match.THREE_PRIME_WINDOW = window  # before the first alignment
    from analyse import load  # noqa: PLC0415

    from assaydrift.analyze import AssayResult  # noqa: PLC0415
    from assaydrift.match import find_oligo  # noqa: PLC0415
    from assaydrift.primers import ASSAYS  # noqa: PLC0415

    out: dict[str, list] = {assay.name: [] for assay in ASSAYS}
    seen = 0
    for record in load(data):
        sequence = record.sequence
        if not sequence:
            continue
        seen += 1
        for assay in ASSAYS:
            result = AssayResult(
                accession=record.accession,
                assay=assay.name,
                period="",
                sequence_length=record.length or len(sequence),
            )
            for oligo in assay.oligos:
                result.hits[oligo.name] = find_oligo(sequence, oligo.sequence, oligo.role)
            if result.region_absent:
                # A fragment that does not reach the target region would be counted as
                # `lost` below and so would land in every shape including "oligo lost
                # only" - the floor the whole comparison rests on. Excluded by the same
                # predicate the published analysis uses, imported rather than restated.
                continue
            out[assay.name].append(
                (
                    result.any_oligo_lost,
                    result.three_prime_mismatches,
                    result.total_mismatches,
                    result.usable,
                )
            )
        if seen % 500 == 0:
            print(f"    window {window}: {seen:,} sequences", flush=True)
    print(f"    window {window}: {seen:,} sequences aligned", flush=True)
    return out


def score(rows: list[tuple[bool, int, int, bool]], load: int, shape: str) -> dict:
    """Re-score already-aligned sequences under one rule shape and load.

    `unknown` sequences are excluded exactly as the published analysis excludes them: an
    alignment with unknown bases under an oligo is not evidence either way, and folding
    them in would move every rate for a reason that has nothing to do with the rule.
    """
    failing = usable_n = 0
    for lost, three_prime, total, usable in rows:
        if not usable and not lost:
            continue
        usable_n += 1
        if lost:
            failing += 1
        elif shape == "oligo lost only":
            pass
        elif shape == "without the 3' clause":
            failing += total >= load
        elif shape == "clustered 3' clause":
            failing += three_prime >= CLUSTERED_MINIMUM or total >= load
        else:
            failing += three_prime > 0 or total >= load
    return wilson(failing, usable_n)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    ap.add_argument("--out", type=Path, default=ROOT / "research" / "results" / "sensitivity.json")
    args = ap.parse_args()
    if not args.data.exists():
        raise SystemExit(f"no sequence store at {args.data}")

    print(
        f"sensitivity grid over {args.data.name}: "
        f"{len(WINDOWS)} window(s) x {len(LOADS)} load(s) x {len(SHAPES)} shape(s)\n"
    )
    grid: dict = {}
    for window in WINDOWS:
        per_assay = features(args.data, window)
        for assay, rows in per_assay.items():
            for load in LOADS:
                for shape in SHAPES:
                    key = f"{assay}|w{window}|l{load}|{shape}"
                    grid[key] = score(rows, load, shape)

    published = {}
    floor = {}
    for assay in {k.split("|")[0] for k in grid}:
        published[assay] = grid[f"{assay}|w5|l3|all three clauses"]
        floor[assay] = grid[f"{assay}|w5|l3|oligo lost only"]

    print("\npublished rule (window 5, load 3, all clauses) against the conservative floor:")
    print(f"  {'assay':16}{'published':>12}{'lost only':>12}{'the claim':>12}")
    for assay in sorted(published):
        p, f = published[assay], floor[assay]
        print(f"  {assay:16}{p['rate']:>12.4f}{f['rate']:>12.4f}{p['rate'] - f['rate']:>+12.4f}")

    print("\nrange across the whole grid, per assay:")
    for assay in sorted(published):
        rates = [v["rate"] for k, v in grid.items() if k.startswith(assay + "|") and v]
        print(
            f"  {assay:16} {min(rates):.4f} to {max(rates):.4f} "
            f"(spread {max(rates) - min(rates):.4f})"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(
            {
                "corpus": args.data.name,
                "windows": list(WINDOWS),
                "loads": list(LOADS),
                "shapes": list(SHAPES),
                "published_rule": {"window": 5, "load": 3, "shape": "all three clauses"},
                "challenge": (
                    "Frontiers Cell Infect Microbiol 2025, 10.3389/fcimb.2025.1524025: assays "
                    "tested in vitro mostly kept working despite primer and probe mismatches, so "
                    "in-silico predictions of false negatives may overstate the risk. "
                    "`likely_failing` is such a prediction; this grid measures how much of it is "
                    "the rule."
                ),
                "conservative_floor": "oligo lost only - the clause a bench would not dispute",
                "grid": grid,
            },
            indent=1,
        ),
        encoding="utf-8",
    )
    print(f"\nwritten to {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
