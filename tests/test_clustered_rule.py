"""The clustered 3' rule, and the distance between it and the published one.

`likely_failing` calls an assay failing when **one** mismatch sits within five bases of a
primer's 3' end. That clause derives from Lefever et al. 2013, whose actual conclusion is
that single mismatches *more than* 5 bp from the 3' end can be tolerated - a statement
about the outside of the window, which the clause inverts into a failure prediction for
the inside of it.

Otwell et al. 2025 measured the inside case directly. A single reverse-primer mismatch
about one nucleotide from the 3' end cost roughly one cycle of Ct and produced no false
negative; two forward-primer mismatches in the same assay shifted Ct by 9-12 cycles and
did; three within the last four nucleotides failed below 5,000 copies.

`clustered_failing` requires two. These tests pin both directions - that one mismatch is
*not* enough and two *are* - because a threshold tested in only one direction passes
whether or not it is wired up, and the rate it produces on the real corpus is zero, which
is exactly the value a rule that never fires would also give.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.analyze import (  # noqa: E402
    CLUSTERED_THREE_PRIME,
    SEVERE_MISMATCH_LOAD,
    AssayResult,
)
from assaydrift.match import Hit  # noqa: E402


def hit(mismatches=0, three_prime=0, found=True, role="forward"):
    return Hit(
        position=10,
        mismatches=mismatches,
        ambiguous=0,
        three_prime_mismatches=three_prime,
        found=found,
        role=role,
    )


def result(**hits) -> AssayResult:
    r = AssayResult(accession="X", assay="test", period="2024-Q1")
    r.hits = hits
    return r


class TestTheThresholdIsTwo:
    def test_the_constant_is_what_the_evidence_supports(self):
        assert CLUSTERED_THREE_PRIME == 2

    def test_one_three_prime_mismatch_fires_the_published_rule(self):
        """The published rule's behaviour is unchanged - this is the baseline."""
        r = result(fwd=hit(mismatches=1, three_prime=1))
        assert r.likely_failing() is True

    def test_one_three_prime_mismatch_does_not_fire_the_clustered_rule(self):
        """The direction that matters: a single 3'-proximal mismatch is tolerated.

        Otwell et al. measured this configuration at about one cycle of Ct shift with no
        false negative, which is why the published rule overstates it.
        """
        r = result(fwd=hit(mismatches=1, three_prime=1))
        assert r.clustered_failing() is False

    def test_two_three_prime_mismatches_fire_the_clustered_rule(self):
        """The other direction. Without this, a rule that never fires would also pass."""
        r = result(fwd=hit(mismatches=2, three_prime=2))
        assert r.clustered_failing() is True

    def test_two_in_one_primer_not_split_across_two(self):
        """`three_prime_mismatches` sums over oligos, so one each in two primers counts.

        Recorded rather than defended. Otwell's two-mismatch false negative was two
        mismatches in the *same* primer, and the summed count cannot distinguish that
        from one in each of two primers. `Hit` keeps a per-oligo count, so a stricter
        rule is expressible; this one is not it, and the limitation belongs in the
        paper's threats rather than in a silent approximation.
        """
        r = result(
            fwd=hit(mismatches=1, three_prime=1),
            rev=hit(mismatches=1, three_prime=1, role="reverse"),
        )
        assert r.three_prime_mismatches == 2
        assert r.clustered_failing() is True


class TestTheOtherClausesAreUnchanged:
    def test_a_lost_oligo_still_fails_under_both(self):
        r = result(fwd=hit(found=False))
        assert r.likely_failing() is True
        assert r.clustered_failing() is True

    def test_a_heavy_load_still_fails_under_both(self):
        r = result(fwd=hit(mismatches=SEVERE_MISMATCH_LOAD))
        assert r.likely_failing() is True
        assert r.clustered_failing() is True

    def test_a_clean_assay_fails_under_neither(self):
        r = result(fwd=hit(), rev=hit(role="reverse"), probe=hit(role="probe"))
        assert r.likely_failing() is False
        assert r.clustered_failing() is False

    def test_the_baseline_is_honoured_by_both(self):
        """An assay with a known reference mismatch must not be charged for it twice."""
        r = result(fwd=hit(mismatches=SEVERE_MISMATCH_LOAD))
        assert r.likely_failing(baseline=1) is False
        assert r.clustered_failing(baseline=1) is False


class TestTheClusteredRuleIsNeverLouder:
    @pytest.mark.parametrize("three_prime", [0, 1, 2, 3])
    @pytest.mark.parametrize("mismatches", [0, 1, 2, 3, 4])
    def test_clustered_implies_published(self, three_prime, mismatches):
        """Anything the clustered rule flags, the published rule flags too.

        The clustered rule only ever narrows clause 2, so it is a subset. If this ever
        fails, the two rules have diverged somewhere other than the 3' threshold and the
        paper's framing of one as the conservative reading of the other is wrong.
        """
        r = result(fwd=hit(mismatches=max(mismatches, three_prime), three_prime=three_prime))
        if r.clustered_failing():
            assert r.likely_failing(), (
                f"clustered fired and published did not at "
                f"three_prime={three_prime}, mismatches={mismatches}"
            )


def test_the_corpus_distribution_is_what_the_paper_says():
    """Every 3'-proximal hit in the shipped corpus is a single mismatch, never a pair.

    This is the observation section 4.2 of the paper turns on, and it is the reason the
    clustered rate equals the no-3'-clause rate. If a future corpus contains a pair, the
    paper's "none is a pair" sentence is false and this fails.
    """
    data = ROOT / "data" / "sequences.jsonl"
    if not data.exists():
        pytest.skip("the corpus is not shipped with the tests")

    sys.path.insert(0, str(ROOT / "scripts"))
    from analyse import load  # noqa: PLC0415

    from assaydrift.match import find_oligo  # noqa: PLC0415
    from assaydrift.primers import ASSAYS  # noqa: PLC0415

    pairs = 0
    singles = 0
    for record in load(data):
        if not record.sequence:
            continue
        for assay in ASSAYS:
            r = AssayResult(accession=record.accession, assay=assay.name, period="")
            for oligo in assay.oligos:
                r.hits[oligo.name] = find_oligo(record.sequence, oligo.sequence, oligo.role)
            if not (r.usable or r.any_oligo_lost):
                continue
            if r.three_prime_mismatches >= 2:
                pairs += 1
            elif r.three_prime_mismatches == 1:
                singles += 1

    assert singles > 0, "no 3'-proximal mismatches at all - the check proves nothing"
    assert pairs == 0, (
        f"{pairs} assay-genome pair(s) carry two or more 3'-proximal mismatches. "
        f"The paper states there are none, and section 4.2's conclusion that the "
        f"clustered rate equals the no-3'-clause rate depends on it."
    )
