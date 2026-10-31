"""The two clauses `clustered_failing` can fire on, counted separately.

`failing_clustered` is a union: a lost oligo, two or more 3'-proximal primer mismatches,
or a severe total load. One number cannot say which clause carried it, and section 4.2's
claim is specifically about 3' pairs - so on the 291,827-genome corpus, where the
clustered count is not zero, the count alone could not tell a pair from a heavy load.
`three_prime_pair` and `severe_load` make the decomposition exact rather than an
inequality inferred from `three_prime` and `failing_clustered`.

These tests pin the wiring in both directions and then reconcile the counters against
`failing_clustered` on whichever shipped artefacts exist, because a counter that is
always zero passes any one-sided check.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.analyze import (  # noqa: E402
    CLUSTERED_THREE_PRIME,
    SEVERE_MISMATCH_LOAD,
    AssayResult,
    summarise,
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
    r = AssayResult(accession="X", assay="A", period="2024-Q1")
    r.hits = dict(hits)
    return r


class TestThreePrimePairCountsPairsOnly:
    def test_one_three_prime_mismatch_is_not_a_pair(self):
        s = summarise([result(f=hit(mismatches=1, three_prime=1))])["2024-Q1"]
        assert s.three_prime == 1
        assert s.three_prime_pair == 0

    def test_two_in_one_primer_is_a_pair(self):
        s = summarise([result(f=hit(mismatches=2, three_prime=2))])["2024-Q1"]
        assert s.three_prime_pair == 1

    def test_two_split_across_two_primers_is_a_pair(self):
        """`three_prime_mismatches` sums over oligos, so one each in two primers counts.

        That is the configuration Otwell et al. measured as a genuine false negative
        (two forward-primer mismatches); splitting them between primers is not evidence
        of safety and the rule does not treat it as such.
        """
        s = summarise(
            [
                result(
                    f=hit(mismatches=1, three_prime=1),
                    r=hit(mismatches=1, three_prime=1, role="reverse"),
                )
            ]
        )["2024-Q1"]
        assert s.three_prime_pair == 1

    def test_a_probe_three_prime_mismatch_never_counts(self):
        """A probe's 3' end carries the quencher and is never extended."""
        s = summarise([result(p=hit(mismatches=2, three_prime=2, role="probe"))])["2024-Q1"]
        assert s.three_prime == 0
        assert s.three_prime_pair == 0

    def test_the_threshold_is_the_published_constant(self):
        below = summarise(
            [
                result(
                    f=hit(
                        mismatches=CLUSTERED_THREE_PRIME - 1, three_prime=CLUSTERED_THREE_PRIME - 1
                    )
                )
            ]
        )["2024-Q1"]
        at = summarise(
            [result(f=hit(mismatches=CLUSTERED_THREE_PRIME, three_prime=CLUSTERED_THREE_PRIME))]
        )["2024-Q1"]
        assert below.three_prime_pair == 0
        assert at.three_prime_pair == 1


class TestSevereLoadCountsLoadOnly:
    def test_below_the_threshold_does_not_count(self):
        s = summarise([result(f=hit(mismatches=SEVERE_MISMATCH_LOAD - 1))])["2024-Q1"]
        assert s.severe_load == 0

    def test_at_the_threshold_counts(self):
        s = summarise([result(f=hit(mismatches=SEVERE_MISMATCH_LOAD))])["2024-Q1"]
        assert s.severe_load == 1

    def test_the_baseline_is_subtracted(self):
        """Known reference mismatches are not drift and must not push a record over."""
        r = [result(f=hit(mismatches=SEVERE_MISMATCH_LOAD))]
        assert summarise(r, baseline=1)["2024-Q1"].severe_load == 0

    def test_load_spread_across_oligos_still_counts(self):
        s = summarise(
            [
                result(
                    f=hit(mismatches=1),
                    r=hit(mismatches=1, role="reverse"),
                    p=hit(mismatches=1, role="probe"),
                )
            ]
        )["2024-Q1"]
        assert s.severe_load == (1 if SEVERE_MISMATCH_LOAD <= 3 else 0)


class TestUnusableRecordsAreNotCounted:
    def test_an_ambiguous_base_excludes_the_record_from_both(self):
        h = Hit(
            position=10,
            mismatches=3,
            ambiguous=1,
            three_prime_mismatches=2,
            found=True,
            role="forward",
        )
        s = summarise([result(f=h)])["2024-Q1"]
        assert s.usable == 0
        assert s.three_prime_pair == 0
        assert s.severe_load == 0

    def test_a_lost_oligo_excludes_the_record_from_both(self):
        """`lost_oligo` is counted, but the record is not usable, so it cannot also
        be counted as a pair or a load - which is why `failing_clustered` over usable
        records decomposes into exactly these two clauses."""
        s = summarise([result(f=hit(found=False, mismatches=9, three_prime=3))])["2024-Q1"]
        assert s.lost_oligo == 1
        assert s.usable == 0
        assert s.three_prime_pair == 0
        assert s.severe_load == 0


ARTEFACTS = ["results/drift.json", "results/drift_big.json"]


@pytest.mark.parametrize("name", ARTEFACTS)
def test_the_counters_reconcile_with_the_clustered_count(name):
    """`failing_clustered` is the union of the two clauses, so it is bounded by them.

    Over usable records a lost oligo cannot fire (a lost oligo makes the record
    unusable), leaving `three_prime_pair` and `severe_load`. The union is therefore at
    least the larger of the two and at most their sum. A clustered count outside that
    range means one of the three is measuring something else.
    """
    path = ROOT / name
    if not path.exists():
        pytest.skip(f"{name} is not present")
    data = json.loads(path.read_text(encoding="utf-8"))
    checked = 0
    for assay, body in data.items():
        if assay == "_run":
            continue
        for period, s in body["periods"].items():
            if "three_prime_pair" not in s:
                pytest.skip(f"{name} predates the clause counters - regenerate it")
            lo = max(s["three_prime_pair"], s["severe_load"])
            hi = s["three_prime_pair"] + s["severe_load"]
            assert lo <= s["failing_clustered"] <= hi, (
                f"{assay} {period}: failing_clustered={s['failing_clustered']} is "
                f"outside [{lo}, {hi}] implied by three_prime_pair="
                f"{s['three_prime_pair']} and severe_load={s['severe_load']}"
            )
            assert s["three_prime_pair"] <= s["three_prime"], (
                f"{assay} {period}: more pairs than records with any 3' mismatch"
            )
            checked += 1
    assert checked > 0, f"{name} carried no periods - nothing was checked"


def test_the_small_corpus_has_no_three_prime_pair_anywhere():
    """Section 4.2's claim, read straight off the artefact rather than re-derived.

    The slow version of this check re-scores every genome in
    `test_clustered_rule.py`; this one fails fast if a regenerated artefact ever
    disagrees with the sentence in the paper.
    """
    path = ROOT / "results" / "drift.json"
    if not path.exists():
        pytest.skip("results/drift.json is not present")
    data = json.loads(path.read_text(encoding="utf-8"))
    pairs = {}
    singles = 0
    for assay, body in data.items():
        if assay == "_run":
            continue
        for period, s in body["periods"].items():
            if "three_prime_pair" not in s:
                pytest.skip("results/drift.json predates the clause counters")
            if s["three_prime_pair"]:
                pairs[f"{assay} {period}"] = s["three_prime_pair"]
            singles += s["three_prime"]
    assert singles > 0, "no 3'-proximal mismatches at all - the check proves nothing"
    assert not pairs, (
        f"the 2,765-genome corpus now carries 3' pairs ({pairs}); section 4.2 says it carries none"
    )
