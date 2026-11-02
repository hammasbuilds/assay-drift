"""Does the sensitivity grid's shape survive the full corpus? Two grids, compared.

    .venv-research/Scripts/python research/grid_corpora.py \
        --small research/results/sensitivity.json \
        --large research/results/sensitivity_big.json

Section 4 reports the grid over 2,765 genomes and section 5 states the limitation: one pass
over 291,827 genomes is about three hours, the grid is three passes, and only the published
setting had been re-scored at full scale. The direction was reported and the magnitudes were
the subset's.

This compares the two. The question is **not** whether the rates moved - they will, over a
corpus a hundred times larger drawn over a different period - but whether the *comparisons*
the paper makes from the grid still hold. Four of them carry section 4:

  * the published rule against the conservative floor, per assay: the size of the claim
    being made on top of the data;
  * the clustered 3' clause against dropping clause 2 altogether, which coincide exactly on
    the small corpus and which section 4.5 says part company at scale;
  * the ordering of the five assays by their published rate;
  * the spread across the grid per assay, which is what the paper quotes as the rule's
    sensitivity to its own parameters.

A cell that moves is reported as a number. A comparison that **reverses** is reported as a
reversal and listed by name, because that is the thing a reader of section 4 would be
misled by and a mean absolute difference would hide it completely.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "research" / "results"

PUBLISHED_WINDOW = 5
PUBLISHED_LOAD = 3
PUBLISHED_SHAPE = "all three clauses"
FLOOR_SHAPE = "oligo lost only"
CLUSTERED_SHAPE = "clustered 3' clause"
DROPPED_SHAPE = "without the 3' clause"


def _fmt(value: float | None, spec: str) -> str:
    """A number or a dash, in one place, so a missing cell cannot print as 0."""
    return "-" if value is None else spec.format(value)


def _rate(grid: dict, key: str) -> float | None:
    cell = grid.get(key)
    return None if cell is None else cell["rate"]


def _assays(grid: dict) -> list[str]:
    return sorted({key.split("|")[0] for key in grid})


def _cell(assay: str, shape: str, window: int = PUBLISHED_WINDOW, load: int = PUBLISHED_LOAD):
    return f"{assay}|w{window}|l{load}|{shape}"


def published_against_floor(small: dict, large: dict) -> dict:
    """The size of the claim, per assay, on each corpus."""
    out = {}
    for assay in _assays(small):
        row = {}
        for name, grid in (("small", small), ("large", large)):
            published = _rate(grid, _cell(assay, PUBLISHED_SHAPE))
            floor = _rate(grid, _cell(assay, FLOOR_SHAPE))
            row[name] = {
                "published": published,
                "floor": floor,
                "claim": (
                    None if published is None or floor is None else round(published - floor, 6)
                ),
            }
        claims = (row["small"]["claim"], row["large"]["claim"])
        row["claim_moved"] = None if None in claims else round(claims[1] - claims[0], 6)
        # The floor is "an oligo has no binding site", which cannot be zero on one corpus
        # and positive on the other for a reason about the rule. Flagged rather than
        # averaged away: it would mean the two corpora differ in how truncated their
        # records are, which is a fact about the corpus and not about the claim.
        row["floor_changed_from_zero"] = (
            None
            if row["small"]["floor"] is None or row["large"]["floor"] is None
            else (row["small"]["floor"] == 0.0) != (row["large"]["floor"] == 0.0)
        )
        out[assay] = row
    return out


def clustered_against_dropped(small: dict, large: dict) -> dict:
    """Does requiring two 3'-proximal mismatches differ from dropping the clause?

    On the small corpus the two columns are identical for every assay, which section 4 says
    is not a coincidence *at that corpus size*. Section 4.5 predicts they part company at
    scale. This is where that prediction is checked.
    """
    out = {}
    for assay in _assays(small):
        row = {}
        for name, grid in (("small", small), ("large", large)):
            clustered = _rate(grid, _cell(assay, CLUSTERED_SHAPE))
            dropped = _rate(grid, _cell(assay, DROPPED_SHAPE))
            row[name] = {
                "clustered": clustered,
                "dropped": dropped,
                "identical": (
                    None
                    if clustered is None or dropped is None
                    else abs(clustered - dropped) < 1e-12
                ),
            }
        out[assay] = row
    identical_small = [a for a, r in out.items() if r["small"]["identical"]]
    identical_large = [a for a, r in out.items() if r["large"]["identical"]]
    return {
        "per_assay": out,
        "identical_on_the_small_corpus": sorted(identical_small),
        "identical_on_the_large_corpus": sorted(identical_large),
        "parted_company_at_scale": sorted(set(identical_small) - set(identical_large)),
        "became_identical_at_scale": sorted(set(identical_large) - set(identical_small)),
    }


def ordering(small: dict, large: dict) -> dict:
    """The assays ranked by their published rate, on each corpus."""

    def rank(grid: dict) -> list[str]:
        pairs = [(assay, _rate(grid, _cell(assay, PUBLISHED_SHAPE))) for assay in _assays(grid)]
        return [a for a, r in sorted(pairs, key=lambda kv: -(kv[1] or -1.0))]

    small_order, large_order = rank(small), rank(large)
    # Named pairs rather than a rank correlation. "Spearman 0.9" over five assays says
    # nothing a reader can act on; "CDC N1 and CDC N2 swapped" does.
    swapped = [
        [a, b]
        for i, a in enumerate(small_order)
        for b in small_order[i + 1 :]
        if large_order.index(a) > large_order.index(b)
    ]
    return {
        "small": small_order,
        "large": large_order,
        "identical": small_order == large_order,
        "pairs_that_swapped": swapped,
    }


def spreads(small: dict, large: dict) -> dict:
    """Range of rates across the whole grid, per assay - the rule's own sensitivity."""
    out = {}
    for assay in _assays(small):
        row = {}
        for name, grid in (("small", small), ("large", large)):
            rates = [
                cell["rate"] for key, cell in grid.items() if key.startswith(assay + "|") and cell
            ]
            row[name] = (
                {
                    "min": round(min(rates), 6),
                    "max": round(max(rates), 6),
                    "spread": round(max(rates) - min(rates), 6),
                }
                if rates
                else None
            )
        out[assay] = row
    return out


def compare(small_payload: dict, large_payload: dict) -> dict:
    small, large = small_payload["grid"], large_payload["grid"]
    shared = sorted(set(small) & set(large))
    moves = [
        {
            "cell": key,
            "small": small[key]["rate"],
            "large": large[key]["rate"],
            "moved": round(large[key]["rate"] - small[key]["rate"], 6),
        }
        for key in shared
        if small[key] and large[key]
    ]
    biggest = sorted(moves, key=lambda row: -abs(row["moved"]))[:10]
    return {
        "corpora": {
            "small": small_payload["corpus"],
            "large": large_payload["corpus"],
        },
        "cells": {
            "in_both": len(shared),
            "only_in_small": sorted(set(small) - set(large)),
            "only_in_large": sorted(set(large) - set(small)),
            "comparable": len(moves),
        },
        "published_against_floor": published_against_floor(small, large),
        "clustered_against_dropped": clustered_against_dropped(small, large),
        "ordering": ordering(small, large),
        "spreads": spreads(small, large),
        "largest_moves": biggest,
        "mean_absolute_move": (
            round(sum(abs(row["moved"]) for row in moves) / len(moves), 6) if moves else None
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, default=RESULTS / "sensitivity.json")
    parser.add_argument("--large", type=Path, default=RESULTS / "sensitivity_big.json")
    parser.add_argument("--out", type=Path, default=RESULTS / "grid_corpora.json")
    args = parser.parse_args()
    for path in (args.small, args.large):
        if not path.exists():
            raise SystemExit(f"no grid at {path}")
    small = json.loads(args.small.read_text(encoding="utf-8"))
    large = json.loads(args.large.read_text(encoding="utf-8"))
    if large.get("partial"):
        raise SystemExit(f"{args.large.name} is a partial grid")
    report = compare(small, large)

    print(f"{report['corpora']['small']} against {report['corpora']['large']}")
    print(f"  {report['cells']['comparable']} of {report['cells']['in_both']} cells comparable")
    print(f"  mean absolute move {report['mean_absolute_move']}")
    print()
    print("published rule against the conservative floor:")
    print(f"  {'assay':16}{'small claim':>14}{'large claim':>14}{'moved':>12}")
    for assay, row in report["published_against_floor"].items():
        # Formatted into plain strings first. Nesting an f-string inside an f-string with
        # the same quote character is a syntax error on Python 3.10, which this package
        # supports - and the research interpreter here is newer, so running the script was
        # no check at all. `ruff check` found it; execution would not have.
        cells = [
            _fmt(row["small"]["claim"], "{:.6f}"),
            _fmt(row["large"]["claim"], "{:.6f}"),
            _fmt(row["claim_moved"], "{:+.6f}"),
        ]
        print(f"  {assay:16}{cells[0]:>14}{cells[1]:>14}{cells[2]:>12}")
    print()
    clustered = report["clustered_against_dropped"]
    print("clustered 3' clause against dropping clause 2 altogether:")
    print(f"  identical on the small corpus: {clustered['identical_on_the_small_corpus']}")
    print(f"  identical on the large corpus: {clustered['identical_on_the_large_corpus']}")
    print(f"  parted company at scale:       {clustered['parted_company_at_scale']}")
    print(f"  became identical at scale:     {clustered['became_identical_at_scale']}")
    print()
    order = report["ordering"]
    print("assays ranked by the published rate:")
    print(f"  small: {' > '.join(order['small'])}")
    print(f"  large: {' > '.join(order['large'])}")
    print(f"  identical ordering: {order['identical']}")
    if order["pairs_that_swapped"]:
        for a, b in order["pairs_that_swapped"]:
            print(f"    swapped: {a} and {b}")
    print()
    print("spread across the grid, per assay:")
    for assay, row in report["spreads"].items():
        small_row, large_row = row["small"], row["large"]
        print(
            f"  {assay:16} small {small_row['spread']:.6f}   large {large_row['spread']:.6f}"
            if small_row and large_row
            else f"  {assay:16} -"
        )
    print()
    print("largest moves:")
    for row in report["largest_moves"][:5]:
        print(f"  {row['cell']:52}{row['small']:.6f} -> {row['large']:.6f}  {row['moved']:+.6f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1), encoding="utf-8")
    print(f"\nwritten to {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
