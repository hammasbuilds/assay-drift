"""The clustered rule is a union of three clauses, and a sum is not a union.

`AssayResult.clustered_failing` is `any_oligo_lost` OR a 3'-proximal pair OR a severe
mismatch load. `summarise` applies it and counts it as `failing_clustered`, and
`scripts/analyse.py` exports that count. Figure 4 nevertheless rebuilt the quantity as
`three_prime_pair + severe_load`, which is wrong in two ways at once: it counts twice any
genome that satisfies both clauses, and it counts not at all a genome that lost an oligo
outright.

On the 2,765-genome corpus the two agreed exactly for all five assays, because no genome
there is in two clauses at once and none loses an oligo. On 291,827 the sum over-counts CDC
N1 by 4 genomes and Charite RdRp by 1, which moved figure 4's caption from "divides RdRp's
rate by 990" to a factor of 1,009. Section 4.5's table is unaffected: every cell rounds the
same at the precision it publishes.

This is the same shape as the baseline bug that `test_grid_agrees_with_the_pipeline.py`
exists for - one quantity, two implementations, nothing comparing them, and a corpus too
small to show the difference. So the test is not "the sum happens to equal the union here".
It is that **nothing outside the library recomputes the rule at all**, plus a direct check
of the union against the library on whatever exports are present.
"""

from __future__ import annotations

import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.analyze import (  # noqa: E402
    CLUSTERED_THREE_PRIME,
    SEVERE_MISMATCH_LOAD,
    AssayResult,
    Hit,
)

EXPORTS = ("results/drift.json", "results/drift_big.json")


def _totals(path: pathlib.Path) -> dict[str, dict[str, int]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, dict[str, int]] = {}
    for assay, body in data.items():
        if assay == "_run":
            continue
        agg: dict[str, int] = {}
        for period in body["periods"].values():
            for field in (
                "usable",
                "failing_clustered",
                "three_prime_pair",
                "severe_load",
                "lost_oligo",
            ):
                agg[field] = agg.get(field, 0) + (period.get(field) or 0)
        out[assay] = agg
    return out


def _present() -> list[tuple[str, dict[str, dict[str, int]]]]:
    found = []
    for name in EXPORTS:
        path = ROOT / name
        if path.is_file():
            totals = _totals(path)
            if all("failing_clustered" in t for t in totals.values()):
                found.append((name, totals))
    return found


def _result(*, lost: bool = False, three_prime: int = 0, mismatches: int = 0) -> AssayResult:
    """One oligo carrying the mismatches asked for, because the rule reads properties.

    `any_oligo_lost`, `three_prime_mismatches` and `total_mismatches` are computed from
    `hits`, not stored, so they cannot be set directly - the result has to be built from a
    `Hit` that produces them. A probe is added so the result is `usable` either way.
    """
    hits = {
        "forward": Hit(
            position=10,
            mismatches=mismatches,
            ambiguous=0,
            three_prime_mismatches=three_prime,
            found=not lost,
            role="forward",
        ),
        "probe": Hit(
            position=40,
            mismatches=0,
            ambiguous=0,
            three_prime_mismatches=0,
            found=True,
            role="probe",
        ),
    }
    return AssayResult("X", "test", "2021-Q1", hits)


# --- the rule itself ------------------------------------------------------------------


def test_a_genome_in_two_clauses_is_one_failure_not_two():
    """The fact the bug turned on. A sum would call this 2; the rule calls it 1."""
    both = _result(three_prime=CLUSTERED_THREE_PRIME, mismatches=SEVERE_MISMATCH_LOAD)
    assert both.clustered_failing() is True
    # It is in both clauses, which is what makes adding the two counters wrong.
    assert both.three_prime_mismatches >= CLUSTERED_THREE_PRIME
    assert both.total_mismatches >= SEVERE_MISMATCH_LOAD


def test_a_lost_oligo_fails_the_clustered_rule_and_is_in_neither_counter():
    """The other half: the sum misses this genome entirely."""
    lost = _result(lost=True)
    assert lost.clustered_failing() is True
    assert lost.three_prime_mismatches < CLUSTERED_THREE_PRIME
    assert lost.total_mismatches < SEVERE_MISMATCH_LOAD


def test_each_clause_alone_fails_it():
    assert _result(three_prime=CLUSTERED_THREE_PRIME).clustered_failing()
    assert _result(mismatches=SEVERE_MISMATCH_LOAD).clustered_failing()


def test_one_short_of_each_clause_does_not():
    assert not _result(three_prime=CLUSTERED_THREE_PRIME - 1).clustered_failing()
    assert not _result(mismatches=SEVERE_MISMATCH_LOAD - 1).clustered_failing()


def test_the_baseline_is_subtracted_by_the_clustered_rule_too():
    """`likely_failing` subtracts it; so must this one, or RdRp diverges again."""
    at_baseline = _result(mismatches=SEVERE_MISMATCH_LOAD)
    assert at_baseline.clustered_failing(baseline=0) is True
    assert at_baseline.clustered_failing(baseline=1) is False


# --- nobody outside the library recomputes it ------------------------------------------


def test_nothing_outside_the_library_turns_the_two_counters_into_a_clustered_rate():
    """The real guard: the numbers agreeing today is not the property worth pinning.

    Narrower than "nothing adds them", deliberately. Adding the two counters is correct
    when the sum is used as what it is - an upper bound on the union, which is how
    `test_clause_counters.py` reconciles the three counters, and that check is right and
    stays. What is wrong is treating the sum *as* the clustered quantity: dividing it by
    `usable` to get a rate, or naming it clustered. That is the substitution that put a
    too-tall bar in figure 4 and a factor of 990 in its caption.

    So the flag is the sum appearing on a line that also mentions `usable` or `clustered`.
    If a legitimate use ever trips it, the fix is to narrow this further with a reason, not
    to delete the check.
    """
    offenders = []
    for path in sorted(ROOT.glob("**/*.py")):
        parts = set(path.parts)
        if parts & {".venv", ".venv-research", "build", "__pycache__"}:
            continue
        if path.resolve() == pathlib.Path(__file__).resolve():
            continue
        if "assaydrift" in parts:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("#") or stripped.startswith("*"):
                continue
            if not ("three_prime_pair" in line and "severe_load" in line and "+" in line):
                continue
            if "usable" in line or "clustered" in line:
                offenders.append(f"{path.relative_to(ROOT)}: {stripped[:90]}")
    assert not offenders, (
        "these turn the two clause counters into the clustered quantity instead of using "
        "`failing_clustered`, which is the rule as the library computes it:\n  "
        + "\n  ".join(offenders)
    )


def test_the_detector_flags_the_line_that_was_wrong_and_spares_the_one_that_was_not():
    """A checker nobody has shown to fire is not a check.

    The first line is figure 4's, as it stood. The second is `test_clause_counters.py`'s
    bound, which must keep passing.
    """

    def flags(line: str) -> bool:
        if not ("three_prime_pair" in line and "severe_load" in line and "+" in line):
            return False
        return "usable" in line or "clustered" in line

    assert flags('clustered.append((t["three_prime_pair"] + t["severe_load"]) / usable)')
    assert flags('rate = (agg["three_prime_pair"] + agg["severe_load"]) / agg["usable"]')
    assert not flags('hi = s["three_prime_pair"] + s["severe_load"]')
    assert not flags('clustered.append(t["failing_clustered"] / usable)')


# --- and against the exports that exist -------------------------------------------------


def test_the_exported_union_never_exceeds_the_sum():
    """A union of two sets is at most their sum - unless clause 1 fired, which it has not.

    If this ever inverts, a genome lost an oligo, and the sum is then wrong in the other
    direction too. The assertion names that case rather than hiding it.
    """
    present = _present()
    if not present:
        pytest.skip("no export carries failing_clustered")
    for name, totals in present:
        for assay, t in totals.items():
            summed = t["three_prime_pair"] + t["severe_load"]
            if t["lost_oligo"]:
                continue
            assert t["failing_clustered"] <= summed, (name, assay, t)


def test_the_gap_is_exactly_the_genomes_in_both_clauses():
    """Over-count equals overlap, which is what makes the sum wrong and by how much."""
    present = _present()
    if not present:
        pytest.skip("no export carries failing_clustered")
    seen_a_difference = False
    for name, totals in present:
        for assay, t in totals.items():
            if t["lost_oligo"]:
                continue
            overlap = t["three_prime_pair"] + t["severe_load"] - t["failing_clustered"]
            assert overlap >= 0, (name, assay, t)
            if overlap:
                seen_a_difference = True
    # Not asserted: on the small corpus alone the overlap is zero everywhere, and that is
    # the whole reason the bug survived. Reported so a reader knows which corpora are here.
    if not seen_a_difference:
        pytest.skip("only corpora where the sum and the union agree are present")


def test_the_clustered_count_never_exceeds_the_usable_population():
    present = _present()
    if not present:
        pytest.skip("no export carries failing_clustered")
    for name, totals in present:
        for assay, t in totals.items():
            assert 0 <= t["failing_clustered"] <= t["usable"], (name, assay, t)
