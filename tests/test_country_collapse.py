"""A wrapped country qualifier must collapse, because the independence figure depends on it.

`_collapse` is one line and was untested: auditing this package found that replacing its
return with `None` leaves the whole suite green. The docstring records that the defect it
exists to prevent really happened - a newline inside `Record.country` split 14,247 of
291,827 rows across two lines in the accession manifest, which exists so the corpus can
be rebuilt exactly.

The consequence runs further than the manifest. `submitter_block` keys on the country, so
a wrapped country and its unwrapped form are two different submitter groups for the same
laboratory. That inflates `clusters` and `effective_n` - and `effective_n` is what decides
which quarters may anchor a published trend. An untested one-line normaliser sits under
the number that qualifies every headline.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift.analyze import PeriodSummary, submitter_block  # noqa: E402
from assaydrift.ncbi import _collapse  # noqa: E402

# The qualifier exactly as GenBank wraps it, from the docstring's own example.
WRAPPED = "Philippines: National Capital Region,\n Makati City"
FLAT = "Philippines: National Capital Region, Makati City"


def test_a_wrapped_qualifier_becomes_one_line() -> None:
    assert _collapse(WRAPPED) == FLAT


def test_the_result_has_no_newline_at_all() -> None:
    """The manifest is line-oriented; one newline in one field corrupts a row."""
    assert "\n" not in _collapse(WRAPPED)
    assert "\r" not in _collapse("USA:\r\n California")


def test_runs_of_whitespace_collapse_to_single_spaces() -> None:
    assert _collapse("USA:   New   York") == "USA: New York"
    assert _collapse("\tIndia:  Delhi \n") == "India: Delhi"


def test_a_plain_value_is_returned_unchanged() -> None:
    """The case the mutant that survives would also pass, so it is not the whole test."""
    assert _collapse("USA: Testland") == "USA: Testland"


def test_an_empty_value_stays_empty() -> None:
    assert _collapse("") == ""
    assert _collapse("   \n  ") == ""


def test_the_same_laboratory_is_one_submitter_block_either_way() -> None:
    """The consequence, not just the string.

    `submitter_block` keys on the country, so if the wrap survived, one laboratory's
    records would count as two independent groups.
    """
    wrapped = submitter_block("MZ123456", _collapse(WRAPPED))
    flat = submitter_block("MZ123457", _collapse(FLAT))
    assert wrapped == flat, (
        "a wrapped and an unwrapped country must give the same submitter block, or one "
        "laboratory is counted as two independent groups"
    )


def test_an_uncollapsed_country_would_inflate_the_independence_figure() -> None:
    """Why the test above is worth having: the arithmetic it protects.

    Two records from one laboratory, grouped as one, are worth one observation. Grouped
    as two they are worth two - and `effective_n` is what decides whether a period may
    anchor a published trend.
    """
    one_group = PeriodSummary("2021-Q2", total=2, usable=2, perfect=2)
    one_group.submitters = {submitter_block("MZ123456", FLAT): 2}

    split = PeriodSummary("2021-Q2", total=2, usable=2, perfect=2)
    split.submitters = {
        submitter_block("MZ123456", FLAT): 1,
        # What the raw, uncollapsed value would key as.
        submitter_block("MZ123457", WRAPPED): 1,
    }
    assert one_group.effective_n == 1.0
    assert split.effective_n == 2.0, (
        "the uncollapsed form must be a different key - otherwise this test proves "
        "nothing about what collapsing buys"
    )
