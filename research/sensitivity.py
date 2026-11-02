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

Only the **window** costs a pass. The 3' count is produced during alignment and `Hit`
keeps the count rather than the positions, so a window cannot be varied after the fact;
the load and the shape can, and are re-scored from the aligned rows in memory. The grid
is therefore three passes, not twelve - a figure the manuscript's section 5 got wrong
until this file was read again.

One pass over the 291,827-genome corpus is about three hours, and the three windows are
independent, so `--window` grids one of them and `--merge` joins the parts. That is how
the large corpus was gridded: three concurrent processes, one wall-clock pass.

    .venv-research/Scripts/python research/sensitivity.py
    .venv-research/Scripts/python research/sensitivity.py --data data/sequences_big.jsonl
    # or, concurrently, one process per window, then a merge:
    ... --data data/sequences_big.jsonl --window 3 --out part-w3.json
    ... --merge part-w3.json,part-w5.json,part-w7.json --out sensitivity_big.json
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
    """Per sequence per assay: (any_oligo_lost, 3' mismatches, EXCESS mismatches, usable).

    The third element is `total_mismatches - baseline_for(assay)`, not the raw total, and
    that is the whole point of this docstring. The production predicate is

        self.total_mismatches - baseline >= SEVERE_MISMATCH_LOAD

    and `baseline_for` is the assay's own mismatch count against the 2019 reference:
    **Charite RdRp shipped with one**, every other assay with none. `baseline_for`'s own
    docstring says scoring RdRp against zero "would report it as broken for every sequence
    ever collected, including the reference itself".

    This file re-implemented the clause as `total >= load` and so did exactly that. Over
    the full corpus it put RdRp's load clause at **406** genomes where the published
    pipeline's `severe_load` is **4** - a hundredfold, on the one assay with a non-zero
    baseline, and invisible on the other four. It was found by comparing the two corpora's
    grids against each other and then against `results/drift_big.json`, not by reading
    either one: the published rate is dominated by the 3' clause, so the error never
    reached the headline column and nothing disagreed until the clause was isolated.

    Aligned once per window, because the 3' count is computed during alignment against
    `match.THREE_PRIME_WINDOW` and `Hit` keeps the count, not the positions - so the
    window cannot be varied after the fact. The mismatch load can, which is why only the
    three window values cost a pass each.
    """
    import assaydrift.match as match

    match.THREE_PRIME_WINDOW = window  # before the first alignment
    from analyse import baseline_for, load  # noqa: PLC0415

    from assaydrift.analyze import AssayResult  # noqa: PLC0415
    from assaydrift.match import find_oligo  # noqa: PLC0415
    from assaydrift.primers import ASSAYS  # noqa: PLC0415

    # One lookup per assay, not per sequence. Imported from the analysis script rather
    # than hard-coded: an assay's baseline is a property of its published oligos, and a
    # constant here would be a second copy of it.
    baselines = {assay.name: baseline_for(assay) for assay in ASSAYS}

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
                    result.total_mismatches - baselines[assay.name],
                    result.usable,
                )
            )
        if seen % 500 == 0:
            print(f"    window {window}: {seen:,} sequences", flush=True)
    print(f"    window {window}: {seen:,} sequences aligned", flush=True)
    return out


def score(rows: list[tuple[bool, int, int, bool]], load: int, shape: str) -> dict:
    """Re-score already-aligned sequences under one rule shape and load.

    The third element of each row is the mismatch load **in excess of the assay's own
    baseline**, so `excess >= load` here is the production predicate
    `total_mismatches - baseline >= SEVERE_MISMATCH_LOAD` with the load varied. See
    `features`.

    `unknown` sequences are excluded exactly as the published analysis excludes them: an
    alignment with unknown bases under an oligo is not evidence either way, and folding
    them in would move every rate for a reason that has nothing to do with the rule.
    """
    failing = usable_n = 0
    for lost, three_prime, excess, usable in rows:
        if not usable and not lost:
            continue
        usable_n += 1
        if lost:
            failing += 1
        elif shape == "oligo lost only":
            pass
        elif shape == "without the 3' clause":
            failing += excess >= load
        elif shape == "clustered 3' clause":
            failing += three_prime >= CLUSTERED_MINIMUM or excess >= load
        else:
            failing += three_prime > 0 or excess >= load
    return wilson(failing, usable_n)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    ap.add_argument("--out", type=Path, default=ROOT / "research" / "results" / "sensitivity.json")
    ap.add_argument(
        "--window",
        type=int,
        default=None,
        choices=WINDOWS,
        help="grid only this 3' window, so the three passes can run concurrently over "
        "a corpus where one pass is hours. Parts join with --merge",
    )
    ap.add_argument(
        "--merge",
        default=None,
        help="comma-separated part files written by --window runs, joined into one "
        "grid. A set that misses a window, or whose parts disagree about the "
        "corpus or about a shared cell, is refused rather than merged",
    )
    args = ap.parse_args()
    if args.merge:
        return merge([Path(q.strip()) for q in args.merge.split(",") if q.strip()], args.out)
    if not args.data.exists():
        raise SystemExit(f"no sequence store at {args.data}")
    windows = WINDOWS if args.window is None else (args.window,)

    print(
        f"sensitivity grid over {args.data.name}: "
        f"{len(windows)} window(s) x {len(LOADS)} load(s) x {len(SHAPES)} shape(s)\n"
    )
    grid: dict = {}
    for window in windows:
        per_assay = features(args.data, window)
        for assay, rows in per_assay.items():
            for load in LOADS:
                for shape in SHAPES:
                    key = f"{assay}|w{window}|l{load}|{shape}"
                    grid[key] = score(rows, load, shape)

    if 5 in windows:
        _tables(grid)
    else:
        print(
            f"\nwindow {windows[0]} only: the published cell is at window 5, so the "
            "comparison tables come from the merge rather than from this part"
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(_payload(args.data.name, windows, grid), indent=1), encoding="utf-8"
    )
    print(f"\nwritten to {args.out.name}")
    return 0


def _tables(grid: dict) -> None:
    """The two comparisons the grid exists for, printed from whatever grid holds them."""
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


def _payload(corpus: str, windows, grid: dict) -> dict:
    return {
        "corpus": corpus,
        "windows": list(windows),
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
    }


def merge(parts: list[Path], out: Path) -> int:
    """Join per-window parts into one grid.

    A set that does not cover every window is refused rather than merged. The headline
    this file carries is a *range* across the grid, and a range over two of three windows
    is narrower than the truth while looking like the result. Two parts measured over
    different corpora would be worse: the rates would be comparable in shape and about
    different populations, so the corpus name has to agree as well.
    """
    loaded = []
    for part in parts:
        if not part.exists():
            raise SystemExit(f"missing part: {part}")
        loaded.append(json.loads(part.read_text(encoding="utf-8")))
    corpora = {q["corpus"] for q in loaded}
    if len(corpora) != 1:
        raise SystemExit(f"parts disagree about the corpus: {sorted(corpora)}")
    covered = sorted({w for q in loaded for w in q["windows"]})
    if covered != sorted(WINDOWS):
        raise SystemExit(f"parts cover windows {covered}, not {sorted(WINDOWS)}")
    grid: dict = {}
    for q in loaded:
        for key, value in q["grid"].items():
            if key in grid and grid[key] != value:
                raise SystemExit(f"parts disagree about cell {key}")
            grid[key] = value
    _tables(grid)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = _payload(loaded[0]["corpus"], WINDOWS, grid)
    payload["merged_from"] = sorted(q.name for q in parts)
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nmerged {len(parts)} part(s) -> {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
