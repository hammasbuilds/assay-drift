"""The README's two Kish effective-n tables, recomputed from the committed exports.

The README publishes two independence tables - one over the 2,765-genome corpus and one
over the 291,827-genome corpus - and nothing checked either. Every cell turned out to be
correct, but two things about the smaller table were not:

  * **It named no assay.** `clusters` and `effective_n` are computed over the records an
    assay can score, so each assay gives a different count; the larger table says "Counted
    for CDC N2" and explains why, and the smaller one silently used CDC N1.
  * **Its first column was `total`, not `usable`.** The effective n beside it is computed
    over usable records, so each row mixed two denominators - 421 deposits beside the
    effective n of the 420 of them CDC N1 can score. The larger table's column is "usable
    sequences" and does not have this problem.

Both are fixed in the README; this pins the result. The point is not that the numbers are
right today - they were - but that `scripts/analyse.py` can change how a quarter is grouped
and nothing would have said so.
"""

from __future__ import annotations

import json
import pathlib
import statistics

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = ROOT / "README.md"

# (quarter, usable, submitter groups, effective n as the README rounds it)
SMALL_ASSAY = "CDC N1"
SMALL = (
    ("2022-Q4", 420, 42, 19.1),
    ("2025-Q4", 324, 23, 11.2),
    ("2020-Q1", 153, 11, 2.4),
    ("2021-Q4", 330, 20, 5.2),
    ("2024-Q1", 159, 4, 1.1),
    ("2026-Q3", 193, 33, 4.6),
)

BIG_ASSAY = "CDC N2"
BIG = (
    ("2020-Q1", 9070, 745, 119.2),
    ("2021-Q4", 38220, 2281, 232.1),
    ("2023-Q3", 508, 57, 12.4),
    ("2026-Q3", 1860, 173, 34.7),
    ("2026-Q4", 67, 7, 1.6),
)

# The three quarters the README excludes from the endpoint for want of independence.
EXCLUDED = (("2026-Q1", 4.0), ("2026-Q2", 4.7), ("2026-Q3", 4.6))

CASES = (
    ("results/drift.json", SMALL_ASSAY, SMALL),
    ("results/drift_big.json", BIG_ASSAY, BIG),
)


def _table_for(assay: str) -> list[str]:
    """The markdown rows of the table introduced by `Counted for **<assay>**`.

    Scoped rather than searched README-wide: a bare search for a quarter finds the first
    line mentioning it, which is a different table's header - the README has several
    tables whose columns are named after quarters.
    """
    lines = README.read_text(encoding="utf-8").splitlines()
    marker = f"Counted for **{assay}**"
    start = next((i for i, ln in enumerate(lines) if marker in ln), None)
    assert start is not None, f"the README never says `{marker}`"
    rows = []
    seen_table = False
    for line in lines[start:]:
        if line.startswith("|"):
            seen_table = True
            rows.append(line)
        elif seen_table and not line.strip():
            break
    assert rows, f"no table follows `{marker}`"
    return rows


def _periods(name: str, assay: str) -> dict | None:
    path = ROOT / name
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if assay not in data:
        return None
    return data[assay]["periods"]


@pytest.mark.parametrize("artefact,assay,rows", CASES)
def test_every_cell_of_the_table_is_in_the_export(artefact: str, assay: str, rows):
    periods = _periods(artefact, assay)
    if periods is None:
        pytest.skip(f"{artefact} absent")
    for quarter, usable, groups, effective in rows:
        got = periods.get(quarter)
        assert got is not None, (artefact, assay, quarter)
        assert got["usable"] == usable, (quarter, "usable", got["usable"], usable)
        assert got["clusters"] == groups, (quarter, "clusters", got["clusters"], groups)
        assert round(got["effective_n"], 1) == effective, (
            quarter,
            "effective_n",
            got["effective_n"],
            effective,
        )


@pytest.mark.parametrize("artefact,assay,rows", CASES)
def test_the_table_is_in_the_readme_as_written(artefact: str, assay: str, rows):
    """Otherwise the test above pins an export against numbers nobody publishes."""
    if _periods(artefact, assay) is None:
        pytest.skip(f"{artefact} absent")
    table = _table_for(assay)
    for quarter, usable, groups, effective in rows:
        line = next((ln for ln in table if quarter in ln), None)
        assert line is not None, f"no row for {quarter} in the {assay} table"
        # Both spellings: the README writes 2,281 with a separator and 42 without.
        for value in (usable, groups):
            assert f"{value:,}" in line or str(value) in line, (quarter, value, line)
        assert str(effective) in line, (quarter, effective, line)


@pytest.mark.parametrize("artefact,assay", [(a, b) for a, b, _ in CASES])
def test_the_readme_names_the_assay_each_table_was_counted_for(artefact: str, assay: str):
    """The defect this file was written for: the smaller table named none.

    `clusters` and `effective_n` are per-assay, so a table without an assay cannot be
    reproduced - and the two tables used different assays.
    """
    if _periods(artefact, assay) is None:
        pytest.skip(f"{artefact} absent")
    readme = README.read_text(encoding="utf-8")
    assert f"Counted for **{assay}**" in readme, (
        f"the {artefact} independence table does not say which assay it was counted for"
    )


def test_the_two_tables_name_different_assays():
    """Which is the whole reason each has to say so."""
    assert SMALL_ASSAY != BIG_ASSAY


@pytest.mark.parametrize("artefact,assay", [(a, b) for a, b, _ in CASES])
def test_no_other_assay_reproduces_the_table(artefact: str, assay: str):
    """A label is only worth pinning if a different assay would give different numbers.

    If every assay gave the same counts the label would be decoration and this file's
    premise would be wrong, so it is checked rather than assumed.
    """
    path = ROOT / artefact
    if not path.is_file():
        pytest.skip(f"{artefact} absent")
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = next(r for a, _, r in CASES if a == artefact)
    matches = []
    for other in data:
        if other == "_run":
            continue
        periods = data[other]["periods"]
        if all(
            (periods.get(q) or {}).get("usable") == u
            and (periods.get(q) or {}).get("clusters") == g
            for q, u, g, _ in rows
        ):
            matches.append(other)
    assert matches == [assay], (artefact, matches)


def test_the_median_effective_n_over_the_reliable_quarters():
    """`The median across the 28 reliable quarters is 82.3`, on the larger corpus."""
    periods = _periods("results/drift_big.json", BIG_ASSAY)
    if periods is None:
        pytest.skip("results/drift_big.json absent")
    reliable = [p["effective_n"] for p in periods.values() if p.get("reliable")]
    assert len(reliable) == 28, len(reliable)
    assert round(statistics.median(reliable), 1) == 82.3, statistics.median(reliable)
    readme = README.read_text(encoding="utf-8")
    assert "28 reliable quarters is **82.3**" in readme


def test_the_three_quarters_excluded_from_the_endpoint():
    """`2026-Q1 at 4.0, 2026-Q2 at 4.7, 2026-Q3 at 4.6` - the same assay as its table."""
    periods = _periods("results/drift.json", SMALL_ASSAY)
    if periods is None:
        pytest.skip("results/drift.json absent")
    readme = README.read_text(encoding="utf-8")
    for quarter, effective in EXCLUDED:
        got = periods.get(quarter)
        assert got is not None, quarter
        assert round(got["effective_n"], 1) == effective, (
            quarter,
            got["effective_n"],
            effective,
        )
        assert f"{quarter} at {effective}" in readme, quarter
