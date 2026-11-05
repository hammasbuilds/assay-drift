"""Section 4.7's generality claims, against `research/results/cohort_generality.json`.

Section 4 derives a structural claim from five assays: the published rate is produced
almost entirely by clause 2 firing on a single 3'-proximal mismatch. Section 4.7 tests it
on eight further published protocols with the same unchanged rule.

The reason this file exists in the shape it does is that **the pooled figure says the
opposite of the per-assay one**. Over the eight new assays together, clause 2 alone
accounts for 7.2% of flagged genomes - which reads as a flat contradiction - because China
CDC N contributes 2,188 of the 2,356 flagged genomes at a 79.4% failing rate. Remove that
one assay and the other seven are 168 flagged, 100% clause 2 alone. A pooled rate over
assays whose denominators differ by three orders of magnitude is a weighted average
dominated by whichever assay is most broken.

So the tests below assert the per-assay structure, and one of them asserts that the pooled
figure is misleading, so that nobody later "simplifies" the section by quoting it.
"""

from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "main.md"
RESULT = ROOT / "research" / "results" / "cohort_generality.json"

# The two assays flagged by mismatch load rather than by 3' position.
LOAD_DOMINATED = {"CDC N1", "China CDC N"}


def _result() -> dict:
    if not RESULT.is_file():
        pytest.skip("run research/cohort_generality.py")
    return json.loads(RESULT.read_text(encoding="utf-8"))


def _rows() -> list[dict]:
    return [r for r in _result()["per_assay"] if r["failing"]]


def _paper() -> str:
    if not PAPER.is_file():
        pytest.skip("the manuscript is not shipped with the tests")
    text = re.sub(r"(?m)^[ \t]*>[ \t]?", "", PAPER.read_text(encoding="utf-8"))
    return re.sub(r"\s+", " ", text)


def test_both_cohorts_were_scored_in_the_same_run():
    """Otherwise the comparison is between two runs, not between two cohorts."""
    payload = _result()
    cohorts = {r["cohort"] for r in payload["per_assay"]}
    assert cohorts == {"primary", "additional"}, cohorts
    assert payload["summary"]["primary"]["assays"] == 5
    assert payload["summary"]["additional"]["assays"] == 8
    assert payload["genomes"] == 2765, payload["genomes"]


def test_eleven_of_thirteen_assays_are_clause_two_alone():
    """The generality claim, per assay rather than pooled."""
    rows = _rows()
    assert len(rows) == 13, [r["assay"] for r in rows]
    clause_two = [r for r in rows if (r["clause_two_share_of_failing"] or 0) >= 0.99]
    assert len(clause_two) == 11, sorted(r["assay"] for r in clause_two)
    # All eleven are exactly 100% on clause 2. The draft said "ten of them exactly 100%",
    # conflating this with the narrower lone-mismatch share below, where HKU N is the one
    # that is not - by two genomes out of 94.
    assert all(r["clause_two_share_of_failing"] == 1.0 for r in clause_two), [
        (r["assay"], r["clause_two_share_of_failing"]) for r in clause_two
    ]
    lone_exact = [r for r in clause_two if r["lone_share_of_failing"] == 1.0]
    assert len(lone_exact) == 10, sorted(r["assay"] for r in lone_exact)
    odd = [r for r in clause_two if r["lone_share_of_failing"] != 1.0]
    assert [r["assay"] for r in odd] == ["HKU N"], odd
    assert odd[0]["failing"] - odd[0]["lone_three_prime"] == 2, odd[0]
    paper = _paper()
    assert "eleven of the thirteen" in paper
    assert "97.9%" in paper and "which is two genomes" in paper


def test_the_two_exceptions_are_load_dominated_not_unexplained():
    """They are clause 3, which §4.1 treats separately - not a counter-example."""
    rows = _rows()
    exceptions = {r["assay"] for r in rows if (r["clause_two_share_of_failing"] or 0) < 0.5}
    assert exceptions == LOAD_DOMINATED, exceptions
    for row in rows:
        if row["assay"] in LOAD_DOMINATED:
            # Nearly every flagged genome carries a severe mismatch load.
            assert row["severe_load"] >= row["failing"] - 1, row


def test_the_pooled_figure_is_misleading_and_the_paper_says_so():
    """The trap this section had to avoid, asserted so it cannot be reintroduced."""
    payload = _result()
    pooled = payload["summary"]["additional"]["clause_two_share"]
    assert pooled is not None and pooled < 0.10, pooled
    paper = _paper()
    assert "7.2%" in paper
    assert "the pooled figure says the opposite" in paper.lower() or (
        "because the pooled figure says the opposite" in paper
    )


def test_one_assay_carries_the_pooled_figure():
    rows = {r["assay"]: r for r in _rows()}
    additional = [r for r in _rows() if r["cohort"] == "additional"]
    total = sum(r["failing"] for r in additional)
    china = rows["China CDC N"]
    assert china["failing"] == 2188, china["failing"]
    assert total == 2356, total
    assert f"{china['failing'] / total * 100:.1f}" == "92.9"
    assert f"{china['failing_rate'] * 100:.1f}" == "79.4"
    paper = _paper()
    assert "2,188 of the 2,356" in paper
    assert "92.9%" in paper and "79.4%" in paper


def test_removing_that_assay_restores_the_result():
    """168 flagged, 100% clause 2 alone - the number the section turns on."""
    rest = [r for r in _rows() if r["cohort"] == "additional" and r["assay"] != "China CDC N"]
    flagged = sum(r["failing"] for r in rest)
    clause_two = sum(r["clause_two_only"] for r in rest)
    assert flagged == 168, flagged
    assert clause_two == flagged, (clause_two, flagged)
    assert "168 flagged, 100% clause 2 alone" in _paper()


def test_the_outlier_is_the_assay_the_literature_already_names():
    """External validation, and the reason to believe the measurement."""
    paper = _paper()
    assert "63.89%" in paper and "86.03%" in paper
    assert "genes2021" in paper


def test_thailand_n_s_small_denominator_is_disclosed():
    """Its single flagged genome must not be read as evidence."""
    rows = {r["assay"]: r for r in _rows()}
    thailand = rows["Thailand N"]
    others = [r["usable"] for k, r in rows.items() if k != "Thailand N"]
    assert thailand["usable"] == 640, thailand["usable"]
    assert thailand["usable"] < min(others) / 4, (thailand["usable"], min(others))
    paper = _paper()
    assert "640 genomes" in paper
    assert "not evidence of anything" in paper


def test_every_per_assay_row_carries_its_denominator():
    """A share with no denominator is the failure mode this whole section is about."""
    for row in _result()["per_assay"]:
        assert row["usable"] > 0, row["assay"]
        if row["failing"]:
            assert row["clause_two_share_of_failing"] is not None, row["assay"]
            assert 0.0 <= row["clause_two_share_of_failing"] <= 1.0, row
