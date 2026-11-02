"""Every number in a generated caption must appear in the paper's copy of it.

`research/run_all.py` writes `research/results/captions.json` beside each figure, and the
paper repeats each caption in prose under the image. Those are two copies of the same
sentence, and nothing compared them, so they drifted: Figure 1's caption said the 3' clause
accounts for "0.1112 of a 0.1116 rate" and Figure 2's said the rate "ranges from 0.1116 to
0.1174", both generated before `research/sensitivity.py` was fixed to subtract each assay's
baseline. The figures were regenerated; the paper was not, because nothing said it had to be.

Equality is the wrong test. The paper's copy is deliberately not identical - it carries
markdown emphasis, typographic dashes, and in Figure 2's case a sentence about why the figure
has two panels, none of which belongs in a caption rendered inside the image. What must match
is every **number**, because a number is what a reader checks against the figure and what goes
stale when the underlying artefact is rebuilt.

So: extract the numeric literals from each generated caption and require each to appear in the
paper's caption for that figure. A number the paper adds of its own is fine; a number the
figure states and the paper contradicts is not.
"""

from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
CAPTIONS = ROOT / "research" / "results" / "captions.json"
PAPER = ROOT / "paper" / "main.md"

# A decimal rate, a thousands-separated count, or a bare integer. Deliberately not matching
# inside a word, so `fig2_sensitivity_heatmap` does not contribute a `2`.
NUMBER = re.compile(r"(?<![\w.])\d{1,3}(?:,\d{3})*(?:\.\d+)?(?![\w.])")

# Numbers too small to be a measurement: panel counts, clause numbers, the window sizes that
# label an axis. They appear in both copies anyway, but requiring them would make the test
# about phrasing rather than about data.
PROSE = {"1", "2", "3", "4", "5", "7", "12"}

FIGURES = {
    "fig1_clause_decomposition": "Figure 1",
    "fig2_sensitivity_heatmap": "Figure 2",
    "fig3_three_prime_distribution": "Figure 3",
    "fig4_big_corpus_decomposition": "Figure 4",
}


def _generated() -> dict[str, str]:
    if not CAPTIONS.is_file():
        pytest.skip("captions.json absent: run research/run_all.py")
    return json.loads(CAPTIONS.read_text(encoding="utf-8"))


def _paper_captions() -> dict[str, str]:
    """The paragraph beginning `**Figure n.**`, up to the next blank line or image."""
    text = PAPER.read_text(encoding="utf-8")
    out = {}
    for label in FIGURES.values():
        match = re.search(rf"\*\*{re.escape(label)}\.\*\*(.*?)(?=\n\n|\n!\[)", text, re.DOTALL)
        if match:
            out[label] = match.group(1)
    return out


def _numbers(text: str) -> set[str]:
    return {n for n in NUMBER.findall(text) if n not in PROSE}


# The one equivalence worth allowing, and no wider than itself: a generated caption writes
# "0 carry two or more" where the paper writes "none carries two or more". That is the same
# measurement in English, and the paper is the better sentence. Anything else - a different
# zero, a rounded value, a count the paper rephrases - still has to match, because a reader
# checking the paper against the figure is checking the digits.
ZERO_IN_WORDS = ("none", "no sequence", "no genome", "zero")


def _satisfied(number: str, paper: str) -> bool:
    if number in _numbers(paper):
        return True
    low = paper.lower()
    return number == "0" and any(word in low for word in ZERO_IN_WORDS)


def test_the_paper_has_a_caption_for_every_generated_figure():
    generated = _generated()
    paper = _paper_captions()
    for key, label in FIGURES.items():
        if key in generated:
            assert label in paper, f"{key} has a generated caption but the paper has no {label}"


@pytest.mark.parametrize("key,label", sorted(FIGURES.items()))
def test_every_number_in_the_generated_caption_is_in_the_paper_s(key: str, label: str):
    generated = _generated()
    if key not in generated:
        pytest.skip(f"{key} not generated")
    paper = _paper_captions().get(label)
    assert paper is not None, f"the paper has no caption for {label}"
    missing = sorted(n for n in _numbers(generated[key]) if not _satisfied(n, paper))
    assert not missing, (
        f"{label}: the figure states {missing} and the paper's caption does not. "
        f"Regenerate with `python research/run_all.py` and copy the caption across."
    )


def test_the_extractor_finds_numbers_and_ignores_the_ones_it_should():
    """Otherwise the test above passes by extracting nothing."""
    found = _numbers("the rate runs from 0.3924 to 0.0004, a factor of 981, over 2,759 usable")
    assert found == {"0.3924", "0.0004", "981", "2,759"}
    assert _numbers("fig2_sensitivity_heatmap and clause 2 of 3") == set()


@pytest.mark.parametrize("key,label", sorted(FIGURES.items()))
def test_the_paper_states_no_number_the_figure_does_not(key: str, label: str):
    """The other direction, which is the one that caught Figure 1.

    Figure 1's stale caption read "0.1112 of a 0.1116 rate" where the figure now says
    "0.1116 of a 0.1116". Checking only that the figure's numbers appear in the paper let
    that through, because 0.1116 was present - the wrong number was the *extra* one. So the
    two sets have to be equal, and they can be: measured over all four figures, the paper's
    captions currently carry no number of their own beyond the typographic ones in `PROSE`.
    If a caption ever needs to cite a figure from elsewhere in the paper, this is the test
    that will say so, and the right fix then is to name it in `PROSE` with a reason rather
    than to drop the direction.
    """
    generated = _generated()
    if key not in generated:
        pytest.skip(f"{key} not generated")
    paper = _paper_captions().get(label)
    assert paper is not None, f"the paper has no caption for {label}"
    invented = sorted(_numbers(paper) - _numbers(generated[key]))
    assert not invented, (
        f"{label}: the paper's caption states {invented} and the figure does not. "
        f"Either it is stale - regenerate with `python research/run_all.py` and copy the "
        f"caption across - or it is a number from elsewhere, which belongs in the body."
    )


def test_the_zero_exemption_is_only_for_zero_and_only_for_the_words():
    """Both directions, because an exemption that is wider than its fact is a hole."""
    assert _satisfied("0", "none carries two or more")
    assert _satisfied("0", "zero for all five assays")
    # Not any other number, even next to the word.
    assert not _satisfied("52", "none carries two or more")
    assert not _satisfied("0.0004", "none carries two or more")
    # Not a zero with no such word in sight.
    assert not _satisfied("0", "the rate is small")
    # And a plain digit match still works without the exemption.
    assert _satisfied("0", "0 carry two or more")


def test_at_least_one_caption_carries_a_measurement():
    """A guard on the guard: if every caption were prose, the parametrised test is vacuous."""
    generated = _generated()
    assert any(_numbers(text) for text in generated.values()), generated.keys()
