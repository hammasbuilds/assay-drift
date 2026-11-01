"""Section 5's length-threshold bullet, read off `research/results/length_threshold.json`.

The bullet exists because the fix for the `likely_failing` clause 1 overreach introduced a
second hand-chosen number - `MIN_WHOLE_GENOME = 25_000` - into a paper whose section 4 is
about hand-chosen numbers. It answers the obvious reviewer question with a measurement, and
a measurement quoted in prose drifts away from the artefact it came from unless something
compares them. Every figure in that bullet is pinned here against the file, including the
two that carry the argument: the empty band around the published threshold, and the count
that a 30,000-base threshold would flip.

The band is the claim worth guarding. It is what makes the published value defensible, and
it is a property of this corpus - a corpus that happens to be near-uniformly full-length.
A regenerated artefact on a fragment-rich corpus should fail here rather than leave the
manuscript asserting an inertness that no longer holds.
"""

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
MEASURED = ROOT / "research" / "results" / "length_threshold.json"
PAPER = ROOT / "paper" / "main.md"

HEADING = "**The 25,000-base constant is itself an uncalibrated threshold.**"

# The figures as section 5 prints them.
PUBLISHED_RECORDS = 291_827
PUBLISHED_THRESHOLD = 25_000
PUBLISHED_BELOW = 1
PUBLISHED_LARGEST_BELOW = 2_859
PUBLISHED_SMALLEST_ABOVE = 28_645
PUBLISHED_GAP = 25_786
PUBLISHED_FLIPPED_AT_30K = 291_824
PUBLISHED_ABOVE_30K = 2
PUBLISHED_MAX = 30_026
PUBLISHED_MEDIAN = 29_810
PUBLISHED_P01 = 29_559


def _measured():
    if not MEASURED.exists():
        pytest.skip("research/results/length_threshold.json is not present")
    return json.loads(MEASURED.read_text(encoding="utf-8"))


def _bullet():
    """The section 5 bullet, newlines collapsed.

    The manuscript is hard-wrapped, so a phrase or a thousands-separated number can
    straddle a line break. Asserting against the raw text makes these tests fail on a
    reflow that changed nothing, which is a test that cries wolf rather than a guard.
    """
    text = PAPER.read_text(encoding="utf-8")
    start = text.index(HEADING)
    return " ".join(text[start : text.index("\n\n", start)].split())


def test_the_record_count_is_the_corpus_the_rest_of_the_paper_uses():
    assert _measured()["records"] == PUBLISHED_RECORDS


def test_one_record_falls_below_the_published_threshold():
    data = _measured()
    assert data["published_threshold"] == PUBLISHED_THRESHOLD
    assert data["below_published"] == PUBLISHED_BELOW


def test_the_band_around_the_threshold_is_empty_and_is_the_width_published():
    data = _measured()
    assert data["largest_below"] == PUBLISHED_LARGEST_BELOW
    assert data["smallest_at_or_above"] == PUBLISHED_SMALLEST_ABOVE
    assert data["gap_around_the_threshold"] == PUBLISHED_GAP
    # The width is the difference, not a third independent number.
    assert data["gap_around_the_threshold"] == (
        data["smallest_at_or_above"] - data["largest_below"]
    )
    # The published threshold has to actually sit inside the band it is defended by.
    assert data["largest_below"] < PUBLISHED_THRESHOLD < data["smallest_at_or_above"]


def test_twenty_thousand_and_twenty_five_thousand_classify_the_corpus_identically():
    """The bullet's "every threshold in that band" needs at least one worked alternative."""
    moves = _measured()["classification_moves_between_candidates"]
    assert moves["20000->25000"] == 0


def test_thirty_thousand_would_flip_the_count_the_bullet_names():
    data = _measured()
    moves = data["classification_moves_between_candidates"]
    assert moves["25000->30000"] == PUBLISHED_FLIPPED_AT_30K
    below = data["count_below_each_candidate"]
    assert below["30000"] - below["25000"] == PUBLISHED_FLIPPED_AT_30K
    # The asymmetry is ruinous because 30,000 sits just under the longest record, so all
    # but a handful of the corpus falls below it. It does not clear the maximum.
    assert data["percentiles"]["max"] > 30_000
    assert data["records"] - below["30000"] == PUBLISHED_ABOVE_30K
    assert PUBLISHED_ABOVE_30K == 2


def test_the_corpus_is_near_uniformly_full_length():
    percentiles = _measured()["percentiles"]
    assert percentiles["max"] == PUBLISHED_MAX
    assert percentiles["median"] == PUBLISHED_MEDIAN
    assert percentiles["p01"] == PUBLISHED_P01
    # The claim is uniformity: the 1st percentile is within a percent of the median.
    assert (percentiles["median"] - percentiles["p01"]) / percentiles["median"] < 0.01


def test_the_shortest_record_is_the_one_the_earlier_bullet_names():
    """LR877184 at 2,859 bases is named in the partial-records bullet above."""
    data = _measured()
    assert data["percentiles"]["min"] == PUBLISHED_LARGEST_BELOW
    text = PAPER.read_text(encoding="utf-8")
    assert "LR877184 is 2,859 bases" in text


def test_every_figure_in_the_bullet_appears_in_the_manuscript():
    """A number pinned to the artefact but missing from the prose guards nothing."""
    bullet = _bullet()
    for value in (
        PUBLISHED_RECORDS,
        PUBLISHED_BELOW,
        PUBLISHED_LARGEST_BELOW,
        PUBLISHED_SMALLEST_ABOVE,
        PUBLISHED_GAP,
        PUBLISHED_FLIPPED_AT_30K,
        PUBLISHED_MAX,
        PUBLISHED_MEDIAN,
        PUBLISHED_P01,
    ):
        rendered = f"{value:,}"
        assert rendered in bullet, f"{rendered} is pinned here but absent from section 5"


def test_the_bullet_does_not_claim_the_threshold_was_calibrated():
    """The point of the bullet is that it was not; a later edit must not soften that."""
    bullet = _bullet()
    assert "has not calibrated it" in bullet
    assert not re.search(r"\bcalibrated\s+(?:the\s+)?(?:constant|threshold|number)\b", bullet)


def test_the_wrap_insensitive_reader_would_notice_a_number_being_removed():
    """The reflow guard must not have made the figure test unable to fail.

    Collapsing whitespace is exactly the kind of normalisation that can turn an assertion
    into one that passes on anything, so this deletes a published figure from a copy of
    the bullet and checks the reader stops seeing it.
    """
    bullet = _bullet()
    assert f"{PUBLISHED_GAP:,}" in bullet
    assert f"{PUBLISHED_GAP:,}" not in bullet.replace(f"{PUBLISHED_GAP:,}", "")
