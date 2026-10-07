"""A base in the SEQUENCE can be ambiguous too, not only `N`.

The comparison used to read `actual in IUPAC[base]`: it expanded the primer's code and
treated the subject's as literal unless it was exactly `N`. A consensus sequence carrying
`R` (A or G) where a primer wants `A` was therefore counted as a definite mismatch -
when R *could be* A. A definite mismatch is what moves an assay towards a drift verdict,
so the error ran in the direction of overstating drift, which is the mistake this
repository has already had to correct twice.

Nothing in the published result moves: the 2,765 sequences behind `results/drift.json`
contain only A, C, G, T and N, and no primer in `ASSAYS` contains an N. Re-running
`scripts/analyse.py` after the change produced a byte-identical `drift.json` apart from
its provenance block. These tests are for the data the tool meets elsewhere - GISAID
consensus sequences do carry R and Y.
"""

from __future__ import annotations

import pytest

from assaydrift.match import AMBIGUOUS, MATCH, MISMATCH, compare_base, find_oligo


@pytest.mark.parametrize(
    ("base", "actual", "expected"),
    [
        # Definite against definite.
        ("A", "A", MATCH),
        ("A", "G", MISMATCH),
        # The subject's possibilities are a subset of what the primer accepts, so it
        # matches whichever base it turns out to be.
        ("R", "A", MATCH),
        ("R", "G", MATCH),
        ("N", "C", MATCH),
        ("R", "R", MATCH),
        # Overlapping but not contained: it MIGHT be the base the primer wants. That is
        # not evidence of drift, and counting it as a mismatch is the bug.
        ("A", "R", AMBIGUOUS),
        ("G", "R", AMBIGUOUS),
        ("A", "N", AMBIGUOUS),
        ("R", "V", AMBIGUOUS),
        # Disjoint: it cannot be what the primer wants, whichever base it is.
        ("A", "Y", MISMATCH),
        ("R", "Y", MISMATCH),
        ("A", "C", MISMATCH),
        # Not a nucleotide code at all.
        ("A", "-", MISMATCH),
        ("A", "?", MISMATCH),
    ],
)
def test_the_three_way_comparison(base: str, actual: str, expected: str) -> None:
    assert compare_base(base, actual) == expected


def test_an_ambiguous_subject_base_is_not_counted_as_a_mismatch() -> None:
    """The case that would have read as drift.

    One `R` inside the oligo's footprint, at a position where the primer wants `A`.
    Before the fix this was one mismatch; it is one ambiguous base, and an ambiguous
    base is an absence of evidence rather than evidence of a changed target.
    """
    oligo = "ACGTACGTACGTACGTACGT"
    genome = "TTTT" + oligo + "TTTT"
    ambiguous_genome = genome.replace(oligo, "RCGTACGTACGTACGTACGT")

    clean = find_oligo(genome, oligo, "forward")
    assert clean is not None
    assert clean.mismatches == 0
    assert clean.ambiguous == 0

    hit = find_oligo(ambiguous_genome, oligo, "forward")
    assert hit is not None, "an ambiguous first base must not stop the oligo being found"
    assert hit.mismatches == 0, "R where the primer wants A is not a definite mismatch"
    assert hit.ambiguous == 1


def test_a_disjoint_ambiguous_base_is_still_a_mismatch() -> None:
    """Ambiguity is not a free pass.

    `Y` is C or T. Where the primer wants `A` it cannot be right whichever base it turns
    out to be, so it counts - otherwise the fix above would hide real drift instead.
    """
    oligo = "ACGTACGTACGTACGTACGT"
    genome = "TTTT" + oligo.replace("A", "Y", 1) + "TTTT"
    hit = find_oligo(genome, oligo, "forward")
    assert hit is not None
    assert hit.mismatches == 1
    assert hit.ambiguous == 0


def test_the_three_prime_window_still_sees_an_ambiguous_base_as_no_mismatch() -> None:
    """The 3' window is where a mismatch stops extension, so it must not be inflated."""
    oligo = "ACGTACGTACGTACGTACGT"
    # Last base, inside the 3' window, made ambiguous rather than wrong. It has to be a
    # code that OVERLAPS the primer's base: the oligo ends in T, and R is A-or-G, which
    # is disjoint from T and so a real mismatch. Y is C-or-T, which might be the T.
    genome = "TTTT" + oligo[:-1] + "Y" + "TTTT"
    hit = find_oligo(genome, oligo, "forward")
    assert hit is not None
    assert hit.three_prime_mismatches == 0, "an ambiguous 3' base is not a 3' mismatch"
    assert hit.ambiguous == 1
