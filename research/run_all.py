"""Regenerate the sensitivity grid's figures and captions.

    .venv-research/Scripts/python research/run_all.py

`.venv-research` is not in the repository (every `.venv*` is ignored) and holds only
matplotlib, so a fresh clone has to create it before this command can run. A sibling
repository had lost its copy without anything failing, because the committed figures still
looked fine:

    uv venv .venv-research --python 3.12
    uv pip install --python .venv-research/Scripts/python matplotlib==3.11.2

Reads `research/results/sensitivity.json`, written by `research/sensitivity.py`, and the
shipped 2,765-genome corpus for the 3'-mismatch distribution. Nothing here re-runs the
291,827-genome analysis.

This file exists because the figures did not have a driver. They were drawn once, by
hand, and were still on disk months later showing a three-shape grid after a fourth
shape had been added - the same staleness that
`tests/test_results_are_not_stale.py` was written for, one directory across.

Nothing in `research/` ships in the distributed package.
"""

from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

import figures

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "research" / "results"
SENSITIVITY = RESULTS / "sensitivity.json"
CORPUS = ROOT / "data" / "sequences.jsonl"


def three_prime_counts() -> dict[str, dict[int, int]]:
    """Per assay, how many genomes carry 0, 1, 2 ... mismatches in a 3' window.

    Recomputed from the corpus rather than read from the grid, because the grid stores
    rates and the claim in section 4.2 is about the *distribution* behind them.
    """
    sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
    import assaydrift.match as match

    match.THREE_PRIME_WINDOW = 5
    from analyse import load  # noqa: PLC0415

    from assaydrift.analyze import AssayResult  # noqa: PLC0415
    from assaydrift.match import find_oligo  # noqa: PLC0415
    from assaydrift.primers import ASSAYS  # noqa: PLC0415

    out: dict[str, collections.Counter] = {a.name: collections.Counter() for a in ASSAYS}
    for record in load(CORPUS):
        if not record.sequence:
            continue
        for assay in ASSAYS:
            result = AssayResult(accession=record.accession, assay=assay.name, period="")
            for oligo in assay.oligos:
                result.hits[oligo.name] = find_oligo(record.sequence, oligo.sequence, oligo.role)
            if result.usable or result.any_oligo_lost:
                out[assay.name][result.three_prime_mismatches] += 1
    return {k: dict(v) for k, v in out.items()}


def main() -> int:
    if not SENSITIVITY.exists():
        raise SystemExit(f"no {SENSITIVITY}: run research/sensitivity.py first")
    payload = json.loads(SENSITIVITY.read_text(encoding="utf-8"))
    grid = payload.get("grid", payload)

    shapes = {k.split("|")[-1] for k in grid}
    assert "clustered 3' clause" in shapes, (
        f"the grid holds shapes {sorted(shapes)} and not the clustered one, so it "
        f"predates the rule the paper's section 4 turns on - re-run sensitivity.py"
    )

    captions = {}
    for name, (_pdf, png, caption) in {
        "fig1_clause_decomposition": figures.clause_decomposition(grid),
        "fig2_sensitivity_heatmap": figures.sensitivity_heatmap(grid),
    }.items():
        captions[name] = caption
        print(f"  {png.name}")

    if CORPUS.exists():
        counts = three_prime_counts()
        for assay, dist in counts.items():
            print(f"    {assay:16} {dict(sorted(dist.items()))}")
        _, png, caption = figures.three_prime_distribution(counts)
        captions["fig3_three_prime_distribution"] = caption
        print(f"  {png.name}")
    else:
        print(f"  (no {CORPUS.name}; skipping figure 3)")

    path = RESULTS / "captions.json"
    path.write_text(json.dumps(captions, indent=1), encoding="utf-8")
    print(f"\nwritten to {ROOT / 'research' / 'figures'} and {path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
