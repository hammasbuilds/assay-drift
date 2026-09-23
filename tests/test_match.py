"""Tests for the matcher, where every wrong answer is a plausible-looking number.

Nothing here touches the network. Sequences are written inline so a failure
means the matching logic is wrong, not that GenBank changed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift.match import (  # noqa: E402
    THREE_PRIME_WINDOW,
    find,
    find_oligo,
    reverse_complement,
    score_at,
)

PAD = "GGATCTTACCGATCAGTTACGGATCAGCTTAAGCCATGCATGCATCAGTA"


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("ATGC") == "GCAT"

    def test_ambiguity_codes_complement_correctly(self):
        # R is A-or-G, so its complement is T-or-C, which is Y.
        assert reverse_complement("R") == "Y"
        assert reverse_complement("S") == "S"  # G-or-C is its own complement
        assert reverse_complement("N") == "N"

    def test_double_reverse_complement_is_identity(self):
        for seq in ("ACGTRYSWKM", PAD, "AAGGTTCC"):
            assert reverse_complement(reverse_complement(seq)) == seq


class TestExactMatching:
    def test_a_perfect_match_scores_zero(self):
        primer = "ACGTACGTACGTACGTACGT"
        hit = find(primer, PAD + primer + PAD)
        assert hit.found
        assert hit.mismatches == 0
        assert hit.position == len(PAD)

    def test_a_primer_that_is_not_there_is_not_found(self):
        # Must report "no binding site", not a low-scoring alignment somewhere.
        hit = find("TTTTTTTTTTTTTTTTTTTT", "ACGACGACGACGACGACGACGACGACGACG")
        assert not hit.found

    def test_a_single_middle_mismatch_is_counted_once(self):
        primer = "ACGTACGTACGTACGTACGT"
        broken = primer[:10] + ("A" if primer[10] != "A" else "C") + primer[11:]
        hit = find(primer, PAD + broken + PAD)
        assert hit.found
        assert hit.mismatches == 1


class TestAmbiguityCodes:
    """A primer's ambiguity code is a real mixture, not a literal letter."""

    @pytest.mark.parametrize(("code", "matches"), [("R", "AG"), ("Y", "CT"), ("S", "GC")])
    def test_an_ambiguous_primer_base_matches_each_of_its_bases(self, code, matches):
        for base in matches:
            primer = "ACGTACG" + code + "TACGTACGT"
            target = "ACGTACG" + base + "TACGTACGT"
            hit = find(primer, PAD + target + PAD)
            assert hit.mismatches == 0, f"{code} should match {base}"

    def test_an_ambiguous_base_still_rejects_what_it_excludes(self):
        # R is A-or-G. It must NOT match C, or the code means nothing.
        hit = find("ACGTACGRTACGTACGT", PAD + "ACGTACGCTACGTACGT" + PAD)
        assert hit.mismatches == 1


class TestAmbiguousTargetBases:
    """`N` in the target is missing data. This is the correctness point.

    Counting N as a mismatch would turn poor sequencing into apparent drift, and
    because sequencing quality changed over time it would produce a *trend*.
    """

    def test_an_n_in_the_target_is_not_a_mismatch(self):
        primer = "ACGTACGTACGTACGTACGT"
        gapped = primer[:5] + "N" + primer[6:]
        hit = find(primer, PAD + gapped + PAD)
        assert hit.mismatches == 0, "an unknown base must not be scored as disagreeing"
        assert hit.ambiguous == 1

    def test_an_alignment_with_unknown_bases_is_not_usable(self):
        primer = "ACGTACGTACGTACGTACGT"
        gapped = primer[:5] + "N" + primer[6:]
        assert not find(primer, PAD + gapped + PAD).usable

    def test_a_clean_alignment_is_usable(self):
        primer = "ACGTACGTACGTACGTACGT"
        assert find(primer, PAD + primer + PAD).usable

    def test_a_fully_unknown_window_is_not_reported_as_a_perfect_match(self):
        """Zero mismatches over unknown bases is not a match, it is no evidence.

        The alignment may start a base early, where a real base in the flank
        happens to agree, so this asserts "almost all unknown" rather than an
        exact count - pinning the offset would be testing the search order
        rather than the property that matters.
        """
        primer = "ACGTACGTACGTACGTACGT"
        hit = find(primer, PAD + "N" * len(primer) + PAD)
        assert hit.ambiguous >= len(primer) - 1
        assert hit.mismatches == 0
        assert not hit.usable, "an all-unknown window must never count as evidence"


class TestThreePrimeEnd:
    """Polymerase extends from the 3' end, so position matters, not just count."""

    def test_a_mismatch_in_the_last_bases_is_flagged(self):
        primer = "ACGTACGTACGTACGTACGT"
        broken = primer[:-1] + ("A" if primer[-1] != "A" else "C")
        hit = find(primer, PAD + broken + PAD)
        assert hit.mismatches == 1
        assert hit.three_prime_mismatches == 1

    def test_a_mismatch_at_the_far_end_is_not_flagged(self):
        primer = "ACGTACGTACGTACGTACGT"
        broken = ("A" if primer[0] != "A" else "C") + primer[1:]
        hit = find(primer, PAD + broken + PAD)
        assert hit.mismatches == 1
        assert hit.three_prime_mismatches == 0, "the 5' end is not the business end"

    def test_the_window_boundary_is_where_it_says(self):
        primer = "ACGTACGTACGTACGTACGT"
        inside = len(primer) - THREE_PRIME_WINDOW
        for index, expected in ((inside, 1), (inside - 1, 0)):
            broken = list(primer)
            broken[index] = "A" if primer[index] != "A" else "C"
            hit = score_at(primer, "".join(broken), 0)
            assert hit.three_prime_mismatches == expected


class TestStrand:
    """A reverse primer appears in a forward-strand record as its reverse
    complement. Searching for it as written finds nothing and reports every
    assay as failed."""

    # Deliberately NOT a palindrome. The first version of this test used
    # ACGTTGCAACGTTGCAACGT, which is its own reverse complement, so it was
    # present on both strands and the strand logic could not be tested at all.
    PRIMER = "ACGTTGCAAGGTTGCAACTT"

    def test_the_fixture_is_not_its_own_reverse_complement(self):
        assert reverse_complement(self.PRIMER) != self.PRIMER

    def test_a_reverse_primer_is_found_via_its_reverse_complement(self):
        genome = PAD + reverse_complement(self.PRIMER) + PAD
        assert find_oligo(genome, self.PRIMER, "reverse").mismatches == 0

    def test_the_same_primer_read_forward_is_not_found(self):
        genome = PAD + reverse_complement(self.PRIMER) + PAD
        assert not find_oligo(genome, self.PRIMER, "forward").found

    def test_a_forward_primer_is_found_as_written(self):
        assert find_oligo(PAD + self.PRIMER + PAD, self.PRIMER, "forward").mismatches == 0


class TestEdges:
    def test_an_empty_oligo_is_not_found(self):
        assert not find("", PAD).found

    def test_a_target_shorter_than_the_oligo_is_not_found(self):
        assert not find("ACGTACGTACGTACGTACGT", "ACGT").found

    def test_score_at_past_the_end_returns_none(self):
        assert score_at("ACGTACGTACGTACGTACGT", "ACGT", 0) is None

    def test_a_primer_at_the_very_start_is_found(self):
        primer = "ACGTACGTACGTACGTACGT"
        hit = find(primer, primer + PAD)
        assert hit.found and hit.position == 0 and hit.mismatches == 0
