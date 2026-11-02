"""Section 5's full-corpus grid paragraph, against the two grids it describes.

Section 5 used to say the sensitivity grid was the 2,765-genome subset and that only the
direction could be reported. The grid has now been run over all 291,827 genomes as well,
and the paragraph states four things about the comparison. Each is recomputed here from
`research/results/sensitivity.json` and `research/results/sensitivity_big.json` rather than
read off `grid_corpora.json`, so a bug in the comparison script cannot make its own summary
agree with itself.

The fourth claim is the one worth having a test for: **one of the magnitudes was the small
corpus's.** CDC N1's published rate spans a factor of about 981 across the twelve threshold
choices on 2,765 genomes and about 212 on 291,827. The contrast figure 2 exists to show
survives - N1 is the assay the thresholds move, Charité RdRp is not - but the smaller
corpus overstates how far by about fourfold, and figure 2 is drawn on it.
"""

from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "main.md"
SMALL = ROOT / "research" / "results" / "sensitivity.json"
BIG = ROOT / "research" / "results" / "sensitivity_big.json"

WINDOWS = (3, 5, 7)
LOADS = (2, 3, 4, 5)
ASSAYS = ("CDC N1", "CDC N2", "CDC N3", "Charite E", "Charite RdRp")


def _grid(path: pathlib.Path) -> dict | None:
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("partial"):
        return None
    return payload["grid"]


def _both() -> tuple[dict, dict]:
    small, big = _grid(SMALL), _grid(BIG)
    if small is None or big is None:
        pytest.skip("both grids are needed; run research/sensitivity.py over each corpus")
    return small, big


def _published(grid: dict, assay: str) -> list[float]:
    """The assay's published-rule rate in all twelve threshold cells."""
    return [
        grid[f"{assay}|w{w}|l{load}|all three clauses"]["rate"] for w in WINDOWS for load in LOADS
    ]


def _paper() -> str:
    if not PAPER.is_file():
        pytest.skip("the manuscript is not shipped with the tests")
    return PAPER.read_text(encoding="utf-8")


def test_the_ordering_by_published_rate_is_identical_across_the_corpora():
    small, big = _both()

    def rank(grid: dict) -> list[str]:
        return sorted(
            ASSAYS,
            key=lambda a: grid[f"{a}|w5|l3|all three clauses"]["rate"],
            reverse=True,
        )

    assert rank(small) == rank(big), (rank(small), rank(big))
    assert rank(big) == ["Charite RdRp", "CDC N2", "CDC N1", "CDC N3", "Charite E"]


def test_the_clustered_column_parts_company_for_exactly_one_assay_at_scale():
    """§4.2's forward reference, confirmed by the grid rather than the published setting."""
    small, big = _both()

    def identical(grid: dict, assay: str) -> bool:
        clustered = grid[f"{assay}|w5|l3|clustered 3' clause"]["rate"]
        dropped = grid[f"{assay}|w5|l3|without the 3' clause"]["rate"]
        return clustered == dropped

    assert all(identical(small, a) for a in ASSAYS), "the small corpus should separate none"
    parted = sorted(a for a in ASSAYS if not identical(big, a))
    assert parted == ["Charite RdRp"], parted


def test_rdrp_is_insensitive_to_both_thresholds_on_both_corpora():
    """All twelve cells equal on the small corpus, and within a rounding step on the large."""
    small, big = _both()
    assert len(set(_published(small, "Charite RdRp"))) == 1, _published(small, "Charite RdRp")
    large = _published(big, "Charite RdRp")
    assert max(large) - min(large) < 5e-4, large
    assert f"{min(large):.4f}" == "0.1828" and f"{max(large):.4f}" == "0.1829", large


def test_cdc_n1s_threshold_swing_is_several_times_smaller_at_scale():
    """The magnitude that was the small corpus's, and the claim that it still survives."""
    small, big = _both()
    s, b = _published(small, "CDC N1"), _published(big, "CDC N1")
    small_factor = max(s) / min(s)
    big_factor = max(b) / min(b)
    assert round(small_factor) == 981, small_factor
    assert round(big_factor) == 212, big_factor
    # It shrinks by about fourfold, and the contrast it illustrates survives: N1 still
    # moves by orders of magnitude where RdRp does not move at all.
    assert 3.5 < small_factor / big_factor < 5.5, small_factor / big_factor
    rdrp = _published(big, "Charite RdRp")
    assert big_factor > 100 > max(rdrp) / min(rdrp), (big_factor, max(rdrp) / min(rdrp))


def test_the_paper_quotes_these_numbers():
    """Otherwise the assertions above pin a grid against a paragraph nobody wrote."""
    _both()
    paper = _paper()
    for value in ("0.1828", "0.1829", "0.1693", "0.0008", "0.3924", "0.0004", "981", "212"):
        assert value in paper, f"§5 no longer quotes {value}"
    assert "240" in paper, "§5 no longer says how many cells were comparable"


def test_the_paper_no_longer_says_the_grid_is_only_the_small_corpus():
    """The claim this measurement retired, and the stale pass count beside it."""
    paper = _paper()
    assert "**The sensitivity grid is the small corpus.**" not in paper
    assert "the grid is twelve passes" not in paper


def test_every_cell_of_both_grids_is_comparable():
    """`240 of 240` - the figure §5 quotes for the comparison's coverage."""
    small, big = _both()
    shared = set(small) & set(big)
    assert len(shared) == 240, len(shared)
    assert len(small) == len(big) == 240, (len(small), len(big))
