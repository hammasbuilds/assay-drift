"""What the manuscript says about other people's papers, pinned to what those papers say.

Every claim here was checked against the article text in November 2026, and two of them
were wrong when checked:

  * **§1 said [@genes2021] reported no independence structure.** It does deduplicate, and
    heavily: 154,602 high-quality genomes collapsed with CD-HIT-EST at a 1.0 identity
    cut-off to 84,305 representatives, so 45.5% of its high-quality set was an exact
    duplicate. The paper's contribution claim rested on that negative, which is the one
    kind of claim a search summary cannot establish. §1 now draws the distinction it should
    have drawn - duplicate removal is not independence weighting - and measures the gap
    between the two on this corpus.
  * **§2 said Primer Monitor tracks GISAID sequences.** Its own pages say NCBI GenBank,
    twice, and name GenBank's citation.

The point of the file is the direction of the check. A test cannot fetch these articles,
so what it pins is that the manuscript still makes the *corrected* claim rather than
drifting back to the convenient one, and that the numbers quoted from each source stay
internally consistent with the arithmetic the manuscript does on them.
"""

from __future__ import annotations

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "main.md"


def _paper() -> str:
    """Flattened, because every check is for a phrase and the source hard-wraps.

    Blockquote markers go first. The citation note is a blockquote, so every line of it
    begins `> `, and flattening whitespace alone leaves those markers stranded mid-phrase -
    `**not that > paper**`. Stripping them is the difference between testing the note's
    wording and testing where its lines happen to break.
    """
    if not PAPER.is_file():
        pytest.skip("the manuscript is not shipped with the tests")
    text = PAPER.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^[ 	]*>[ 	]?", "", text)
    return re.sub(r"\s+", " ", text)


# --- what the sources say -----------------------------------------------------------


def test_the_genes2021_figures_are_the_articles():
    """Eleven sets, 84,305 of 154,602, the dates, and the China CDC series."""
    paper = _paper()
    for fact in (
        "eleven primer/probe sets",
        "84,305 genomes",
        "June 2020 to January 2021",
        "154,602",
        "CD-HIT-EST",
        "63.89%",
        "86.03%",
        "74.54%",
    ):
        assert fact in paper, f"the manuscript no longer states {fact!r} for @genes2021"


def test_the_deduplication_share_is_arithmetic_not_assertion():
    """45.5% has to be 1 - 84,305/154,602, or the sentence is decoration."""
    share = (154602 - 84305) / 154602
    assert f"{share * 100:.1f}%" == "45.5%"
    assert "45.5%" in _paper()


def test_primer_monitor_is_credited_to_genbank_not_gisaid():
    """The error this file was written for. NEB's own pages say NCBI GenBank."""
    paper = _paper()
    assert "NCBI GenBank** sequences against user-supplied primer sets" in paper
    assert "Primer Monitor tracks GISAID" not in paper
    assert "Primer Monitor tracks **GISAID**" not in paper


def test_the_validation_websites_three_criteria_are_quoted():
    """§2 leans on them to argue the rule family is shared, not idiosyncratic."""
    paper = _paper()
    for clause in (
        "three or more mismatches in an oligo/target alignment",
        "melting temperature below 40",
        "last two 3′ positions of a primer",
    ):
        assert clause in paper, f"§2 no longer quotes {clause!r}"


# --- the claim that was overstated, and the shape of its replacement ------------------


def test_the_paper_no_longer_claims_prior_work_reports_no_independence_structure():
    """The exact sentence that was false, and the family of ways back to it."""
    paper = _paper()
    assert "Neither reports an **independence structure**" not in paper
    assert "both count sequences" not in paper


def test_the_paper_distinguishes_deduplication_from_independence_weighting():
    """What replaced it, and the distinction the contribution now rests on."""
    paper = _paper()
    assert "duplicate removal is not independence weighting" in paper or (
        "**removing duplicates** and **weighting for independence**" in paper
    )
    assert "does deduplicate" in paper


def test_the_gap_between_the_two_is_measured_on_this_corpus():
    """A distinction worth drawing is worth a number, and these are the numbers."""
    paper = _paper()
    assert "**159 sequences to 153**" in paper
    assert "3.8%" in paper
    assert "effective sample size is **1.1**" in paper
    # And the arithmetic behind the quoted reduction.
    assert f"{(159 - 153) / 159 * 100:.1f}%" == "3.8%"


def test_the_corpus_wide_duplicate_share_is_quoted_so_one_quarter_cannot_carry_it():
    """Per `weak result smells`: one quarter cannot establish a corpus-scale claim."""
    assert "4.0% of sequences are an exact duplicate" in _paper()


def test_the_fifth_failed_absence_claim_is_recorded():
    """The paper counts these deliberately; the count must move when another fails."""
    paper = _paper()
    assert "fourth absence-based" in paper
    assert "A fifth has since failed" in paper


# --- the citation note must not overstate what was read -------------------------------


def test_the_citation_note_says_which_sources_were_read_and_which_were_not():
    paper = _paper()
    assert "@genes2021 **has now been read**" in paper
    assert "@bioinformatics2020 **has now been read**" in paper
    # The two described from primary sources short of a full text, each named as such.
    assert "from its abstract and its authors' own repository" in paper
    assert "The full text remains unread and no number in this paper comes from it" in paper
    assert "described **from the tool's own pages**" in paper


def test_assaym_is_not_described_with_the_detail_no_source_supports():
    """An earlier draft gave it a genome browser; neither the abstract nor the repo does."""
    paper = _paper()
    assert "maps primer and probe regions onto a genome browser" not in paper.replace(
        'assayM "maps primer and probe regions onto a genome browser"', ""
    )
    assert "primer and probe sequences published by the WHO" in paper


def test_the_note_records_the_misattributed_preprint_doi():
    """A DOI circulating for the Primer Monitor preprint resolves to another paper."""
    paper = _paper()
    assert "10.1101/2022.04.22.489186" in paper
    assert "is **not that paper**" in paper
    # And it must not appear as a citation anywhere.
    assert paper.count("10.1101/2022.04.22.489186") == 1


def test_lefever_is_still_flagged_as_abstract_only():
    """Its results section is paywalled and two claims in §4 rest on it."""
    paper = _paper()
    assert "@lefever2013" in paper
    assert "paywalled" in paper


# --- the TaqMAMA provenance -----------------------------------------------------------


def test_the_two_base_rule_is_traced_to_a_genotyping_method():
    """§4.1's strongest piece of provenance, and the one most worth protecting.

    [@bioinformatics2020] justifies calling any mismatch in a primer's last two 3' bases a
    predicted false negative by citing "Li et al. (2004)". That reference is *Genotyping
    with TaqMAMA* - TaqMan plus the Mismatch Amplification Mutation Assay - where a
    3'-terminal mismatch is deliberately engineered so the primer fails on the allele it is
    not meant to detect, and the Ct shift is the discrimination you are trying to maximise.
    Both halves were read: the citing sentence from the Bioinformatics full text, and what
    TaqMAMA is from the Genomics abstract.
    """
    paper = _paper()
    assert "[@taqmama]" in paper, "the TaqMAMA reference is gone"
    assert "Genotyping with TaqMAMA" in paper
    assert "Mismatch Amplification Mutation Assay" in paper
    assert "deliberately engineered" in paper
    # The quoted sentence from the citing paper, and the two from the cited one.
    assert "increasing detection Ct ≥2" in paper
    assert "strength of allelic discrimination" in paper
    assert "directed use of nucleotide mismatches" in paper


def test_the_paper_does_not_overclaim_what_taqmama_measured():
    """Its full text is subscription-only; the argument is about what the study was for."""
    paper = _paper()
    assert "Its full text is subscription-only and was not read" in paper
    assert "only about what the study was *for*" in paper


def test_the_paper_separates_that_rule_from_its_own():
    """`likely_failing` uses five bases, not two. Conflating them would be the easy error."""
    paper = _paper()
    assert "it is not this paper's rule" in paper
    assert "`likely_failing` uses five bases, not two" in paper


def test_the_taqmama_reference_carries_its_doi_and_its_role():
    paper = _paper()
    assert "10.1016/j.ygeno.2003.08.005" in paper
    assert "Genomics 83(2):311–320, 2004" in paper
