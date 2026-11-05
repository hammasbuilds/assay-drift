"""The README's by-year table, against `results/drift_by_year.json`.

The README publishes a quarterly view and a yearly view of the same 2,765 sequences. Only
the quarterly one had an artefact in the repository, so the yearly table was the one set of
published rates nothing could check - and four of its cells had gone stale. Charite E's
2026 failing rate read 0.9% where the measurement says 0.6%, and CDC N2's three 2026 cells
read 68.8%, 15.8% and 15.8% where it says 69.4%, 15.3% and 15.3%.

The arithmetic identifies the direction. CDC N2's 2026 figures are now 51 failing of 333
usable; the published 15.8% needs 53 of 335. Two records left the usable population and
both had been counted failing, which is what an ambiguity correction does: a sequence with
a degenerate base under an oligo is not evidence either way, so it leaves the numerator and
the denominator together. 2026 is the smallest year, so two records move it and no other
year visibly.

Running `scripts/analyse.py --granularity year` is what produced the artefact, and that run
has a side effect worth knowing about: it also rewrites `results/mutations.json` from the
same pass, at whatever granularity it was given. Doing it without noticing replaces a
committed quarterly artefact with a yearly one. It was noticed and reverted here; the
`--json` path is the only output that should be redirected casually.
"""

from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
BY_YEAR = ROOT / "results" / "drift_by_year.json"

# assay -> (2020 exact, 2026 exact, 2026 failing, worst failing year)
PUBLISHED = {
    "CDC N1": ("99.0%", "0.0%", "7.0%", "7.0%", "2026"),
    "Charite E": ("98.8%", "0.0%", "0.6%", "0.6%", "2026"),
    "Charite RdRp": ("98.3%", "87.0%", "13.0%", "38.3%", "2021"),
    "CDC N3": ("99.3%", "90.3%", "0.9%", "2.0%", "2023"),
    "CDC N2": ("98.0%", "69.4%", "15.3%", "15.3%", "2026"),
}


def _years() -> dict:
    if not BY_YEAR.is_file():
        pytest.skip("results/drift_by_year.json is not present")
    return json.loads(BY_YEAR.read_text(encoding="utf-8"))


def _readme() -> str:
    if not README.is_file():
        pytest.skip("README.md is not shipped with the tests")
    return README.read_text(encoding="utf-8")


HEADING = "### The same data by year"


def _year_table() -> list[str]:
    """The rows of the by-year table, found by its own heading.

    The README holds several tables whose rows begin with an assay name, so "the first
    line mentioning Charite RdRp" is the quarterly table - which is how the first version
    of this test failed, comparing yearly expectations against quarterly rates.
    """
    lines = _readme().splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.strip() == HEADING), None)
    assert start is not None, f"the README no longer has the heading {HEADING!r}"
    rows, seen = [], False
    for line in lines[start:]:
        if line.startswith("|"):
            seen = True
            rows.append(line)
        elif seen and not line.strip():
            break
    assert rows, "no table follows the by-year heading"
    return rows


def _pct(value: float | None) -> str | None:
    return None if value is None else f"{value * 100:.1f}%"


@pytest.mark.parametrize("assay", sorted(PUBLISHED))
def test_the_published_year_rates_are_the_measured_ones(assay: str):
    data = _years()
    assert assay in data, assay
    periods = data[assay]["periods"]
    first, last, failing, _worst, _year = PUBLISHED[assay]
    assert _pct(periods["2020"]["perfect_rate"]) == first, (assay, "2020 exact")
    assert _pct(periods["2026"]["perfect_rate"]) == last, (assay, "2026 exact")
    assert _pct(periods["2026"]["failing_rate"]) == failing, (assay, "2026 failing")


@pytest.mark.parametrize("assay", sorted(PUBLISHED))
def test_the_worst_year_column_is_the_worst_year(assay: str):
    """A maximum, recomputed - not the last year that happened to be noticed."""
    periods = _years()[assay]["periods"]
    rates = {y: s["failing_rate"] for y, s in periods.items() if s.get("failing_rate") is not None}
    worst_year = max(rates, key=lambda y: rates[y])
    _f, _l, _fail, worst, year = PUBLISHED[assay]
    assert worst_year == year, (assay, worst_year, year)
    assert _pct(rates[worst_year]) == worst, (assay, rates[worst_year], worst)


@pytest.mark.parametrize("assay", sorted(PUBLISHED))
def test_the_readme_row_says_what_the_artefact_says(assay: str):
    """Otherwise the table above pins an artefact to numbers nobody publishes."""
    label = assay.replace("Charite", "Charité")
    row = next((ln for ln in _year_table() if label in ln), None)
    assert row is not None, f"no row for {label} in the by-year table"
    for value in PUBLISHED[assay][:4]:
        assert value in row, (assay, value, row)


def test_the_year_artefact_is_the_same_corpus_as_the_quarterly_one():
    """Two views of one corpus; if they ever diverge, the comparison above is meaningless."""
    quarterly = ROOT / "results" / "drift.json"
    if not quarterly.is_file():
        pytest.skip("results/drift.json is not present")
    a = json.loads(quarterly.read_text(encoding="utf-8"))["_run"]
    b = _years()["_run"]
    assert b["granularity"] == "year", b["granularity"]
    assert a["granularity"] == "quarter", a["granularity"]
    assert a["sequences_read"] == b["sequences_read"], (a["sequences_read"], b["sequences_read"])
    for key in ("corpus_identity_sha256", "manifest_records"):
        if key in a and key in b:
            assert a[key] == b[key], key


def test_every_year_total_is_the_sum_of_its_quarters():
    """The two groupings must partition the same records, or one of them drops some.

    `period_of` returns a year whenever it can find one but `None` for a quarter when no
    month can be read, so a corpus of year-only dates would put records in the yearly view
    and not the quarterly. Over both corpora here that count is zero, and this is the
    assertion that would notice if a future corpus changed it.
    """
    quarterly = ROOT / "results" / "drift.json"
    if not quarterly.is_file():
        pytest.skip("results/drift.json is not present")
    q = json.loads(quarterly.read_text(encoding="utf-8"))
    y = _years()
    for assay in PUBLISHED:
        per_year: dict[str, int] = {}
        for period, summary in q[assay]["periods"].items():
            per_year[period[:4]] = per_year.get(period[:4], 0) + (summary.get("usable") or 0)
        for year, total in per_year.items():
            assert y[assay]["periods"][year]["usable"] == total, (assay, year)


def test_the_readme_records_why_the_table_moved():
    """A corrected number with no account of the correction invites the old one back."""
    readme = re.sub(r"\s+", " ", _readme())
    assert "51 failing of 333 usable" in readme
    assert "needs 53 of 335" in readme
    assert "drift_by_year.json" in readme
