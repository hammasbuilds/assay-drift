"""The sensitivity grid against the production pipeline, cell by cell.

Two independent pieces of code compute the same rule over the same corpus. `analyse.py`
runs `AssayResult.likely_failing` and writes `results/drift*.json`; `research/sensitivity.py`
re-implements the clauses so it can vary the window, the load and the shape without
re-aligning. Nothing compared them, and they disagreed.

The production predicate is

    self.total_mismatches - baseline >= SEVERE_MISMATCH_LOAD

where `baseline_for(assay)` is the assay's own mismatch count against the 2019 reference:
**Charite RdRp shipped with one**, every other assay with none, and `baseline_for`'s own
docstring warns that scoring RdRp against zero "would report it as broken for every
sequence ever collected, including the reference itself". The grid's `score()` had
`total >= load`, with no baseline, so it did exactly that for RdRp alone.

The size of the error is the reason this file exists. On the 2,765-genome corpus the paper
reports, it was **one genome** - RdRp's load clause read 1 where the pipeline says 0.
Nothing in the published tables moved, because the published rate is dominated by the 3'
clause and the load clause only shows when it is isolated. At 291,827 genomes the same bug
read **406** where the pipeline says **4**. It was found by running the grid at full scale
and comparing the two, which is precisely what section 5 had said was not yet done.

Five identities per assay per corpus, from `figures.py`'s own decomposition:

  * grid `n`                                 == pipeline `usable`
  * `all three clauses` failing              == pipeline `failing`
  * `without the 3' clause` failing          == pipeline `severe_load`
  * `clustered 3' clause` failing            == pipeline `three_prime_pair + severe_load`
  * `oligo lost only` failing                == pipeline `lost_oligo`

They are checked at the published setting (window 5, load 3), which is the only cell the
pipeline computes - it does not vary the window or the load, which is why the grid exists.
"""

from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]

PUBLISHED_WINDOW = 5
PUBLISHED_LOAD = 3
FIELDS = (
    "usable",
    "failing",
    "failing_clustered",
    "three_prime",
    "three_prime_pair",
    "severe_load",
    "lost_oligo",
)

# (grid artefact, pipeline artefact) for each corpus the repository has measured.
CORPORA = (
    ("research/results/sensitivity.json", "results/drift.json"),
    ("research/results/sensitivity_big.json", "results/drift_big.json"),
)


def _load(path: str) -> dict | None:
    full = ROOT / path
    if not full.is_file():
        return None
    return json.loads(full.read_text(encoding="utf-8"))


def _pipeline_totals(drift: dict) -> dict[str, dict[str, int]]:
    """Per-assay totals, summed over the run's periods."""
    out: dict[str, dict[str, int]] = {}
    for assay, body in drift.items():
        if assay == "_run":
            continue
        agg = dict.fromkeys(FIELDS, 0)
        for period in body["periods"].values():
            for field in FIELDS:
                agg[field] += period.get(field) or 0
        out[assay] = agg
    return out


def _cell(grid: dict, assay: str, shape: str) -> dict:
    key = f"{assay}|w{PUBLISHED_WINDOW}|l{PUBLISHED_LOAD}|{shape}"
    assert key in grid, f"the grid has no cell {key}"
    return grid[key]


def _pairs():
    """Every (corpus, grid, pipeline totals) the repository can compare today."""
    found = []
    for grid_path, drift_path in CORPORA:
        grid = _load(grid_path)
        drift = _load(drift_path)
        if grid is None or drift is None:
            continue
        if grid.get("partial"):
            continue
        found.append((grid_path, grid["grid"], _pipeline_totals(drift)))
    return found


def test_at_least_one_corpus_can_be_compared():
    """Otherwise every test below skips and the file guards nothing."""
    pairs = _pairs()
    if not pairs:
        pytest.skip("neither a grid nor a pipeline artefact is present")
    assert pairs


@pytest.mark.parametrize(
    "shape,field",
    [
        ("all three clauses", "failing"),
        ("without the 3' clause", "severe_load"),
        ("oligo lost only", "lost_oligo"),
    ],
)
def test_the_grid_and_the_pipeline_agree_on_every_clause(shape: str, field: str):
    pairs = _pairs()
    if not pairs:
        pytest.skip("nothing to compare")
    for corpus, grid, pipeline in pairs:
        for assay, totals in pipeline.items():
            cell = _cell(grid, assay, shape)
            assert cell["failing"] == totals[field], (corpus, assay, shape, field)


def test_the_clustered_clause_is_the_pipeline_s_own_clustered_count():
    """Against `failing_clustered`, which is the rule, not against the clauses added up.

    This compared the grid to `three_prime_pair + severe_load` when it was written,
    because that is what `figures.py` did. Both were wrong in the same way: the clustered
    rule is a *union* of three clauses, so adding two of them double-counts a genome in
    both and omits one that lost an oligo. The grid itself was always right - it evaluates
    `lost or three_prime >= 2 or excess >= load`, which is the union - so this test was
    asserting the grid equals a quantity the grid correctly does not compute. It passed
    only because the 2,765-genome corpus has no genome in two clauses.

    On 291,827 it would have failed, and that is the one useful thing about how it was
    written: it would have caught figure 4 rather than agreeing with it. Comparing against
    `failing_clustered` catches it on either corpus. See
    `test_clustered_is_the_union_not_the_sum.py`.
    """
    pairs = _pairs()
    if not pairs:
        pytest.skip("nothing to compare")
    for corpus, grid, pipeline in pairs:
        for assay, totals in pipeline.items():
            assert (
                _cell(grid, assay, "clustered 3' clause")["failing"] == totals["failing_clustered"]
            ), (corpus, assay)


def test_both_measure_the_same_population():
    """A rate can only be compared against another rate over the same denominator."""
    pairs = _pairs()
    if not pairs:
        pytest.skip("nothing to compare")
    for corpus, grid, pipeline in pairs:
        for assay, totals in pipeline.items():
            assert _cell(grid, assay, "all three clauses")["n"] == totals["usable"], (
                corpus,
                assay,
            )


def test_the_baseline_is_subtracted_and_only_one_assay_has_one():
    """The fact the bug turned on, asserted rather than remembered.

    If a second assay ever acquires a baseline, every cell of its load clause moves, and
    the test above is what would catch it - but only if the grid still subtracts it.
    """
    import sys

    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(ROOT / "scripts"))
    from analyse import baseline_for

    from assaydrift.primers import ASSAYS

    baselines = {assay.name: baseline_for(assay) for assay in ASSAYS}
    assert baselines["Charite RdRp"] == 1, baselines
    assert sum(1 for value in baselines.values() if value) == 1, baselines

    source = (ROOT / "research" / "sensitivity.py").read_text(encoding="utf-8")
    assert "baseline_for" in source, "the grid no longer subtracts the assay's baseline"
    assert "result.total_mismatches - baselines[assay.name]" in source


def test_the_published_cell_is_the_only_one_the_pipeline_computes():
    """So the grid is not redundant: the pipeline has one cell and the grid has 48."""
    pairs = _pairs()
    if not pairs:
        pytest.skip("nothing to compare")
    for corpus, grid, _pipeline in pairs:
        windows = {key.split("|")[1] for key in grid}
        loads = {key.split("|")[2] for key in grid}
        assert len(windows) == 3, (corpus, windows)
        assert len(loads) == 4, (corpus, loads)
