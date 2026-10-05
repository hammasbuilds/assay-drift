"""An indel in a primer site must read as an indel, not as a run of mismatches."""

from __future__ import annotations

import pytest

from assaydrift.match import find_oligo

PRIMER = "ACCAGGAACTAATCAGACAAG"
# Neither flank may contain the primer, or the matcher finds the clean copy and the
# mutated one is never scored - which is how an earlier version of this test "passed".
LEFT = "TTGCATTACGTTTGGTGGTCCCTCTGATTCTTCTGGCTGTTACC"
RIGHT = "GGTTTTTGGTTTCTGTTACTTTGGTTCTGATTTCTTTCTTTCTT"


def genome(middle: str) -> str:
    assert PRIMER not in LEFT + RIGHT
    return LEFT + middle + RIGHT


@pytest.mark.parametrize(
    ("label", "middle", "indels", "mismatches"),
    [
        ("1-base deletion", PRIMER[:10] + PRIMER[11:], 1, 0),
        ("2-base deletion", PRIMER[:10] + PRIMER[12:], 2, 0),
        ("3-base deletion", PRIMER[:9] + PRIMER[12:], 3, 0),
        ("3-base insertion", PRIMER[:10] + "GGG" + PRIMER[10:], 3, 0),
        ("deletion plus a substitution",
         PRIMER[:10] + PRIMER[11:15] + "T" + PRIMER[16:], 1, 1),
    ],
)
def test_an_indel_is_reported_as_an_indel(
    label: str, middle: str, indels: int, mismatches: int
) -> None:
    """One deleted base used to read as eight mismatches, three of them at the 3' end.

    That is the worst possible way to be wrong here: a 3' mismatch is the signal this
    tool reports as stopping the reaction, so a deletion in the MIDDLE of a primer site
    was reported as the primer's business end being destroyed.
    """
    hit = find_oligo(genome(middle), PRIMER, "forward")
    assert hit.indels == indels, label
    assert hit.mismatches == mismatches, label
    # The point of aligning before counting: a mid-primer indel must not manufacture
    # mismatches inside the 3' window.
    assert hit.three_prime_mismatches == 0, label


@pytest.mark.parametrize(
    ("label", "middle", "mismatches", "three_prime"),
    [
        ("exact", PRIMER, 0, 0),
        ("one substitution in the middle", PRIMER[:10] + "T" + PRIMER[11:], 1, 0),
        ("one substitution at the 3' end", PRIMER[:-1] + "T", 1, 1),
        ("two substitutions at the 3' end", PRIMER[:-2] + "TT", 2, 2),
    ],
)
def test_substitution_drift_is_not_refitted_as_gaps(
    label: str, middle: str, mismatches: int, three_prime: int
) -> None:
    """Substitution-only drift is the common case and must stay on the ungapped path.

    An aligner that opens a gap to save a single mismatch would turn ordinary drift
    into fictional deletions, and the 3'-end rule would be evaluated on the wrong
    positions. The affine gap cost is what keeps these ungapped.
    """
    hit = find_oligo(genome(middle), PRIMER, "forward")
    assert hit.indels == 0, label
    assert hit.mismatches == mismatches, label
    assert hit.three_prime_mismatches == three_prime, label
