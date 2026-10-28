"""Tests for the grouping and rate logic.

Every bug this file guards against produces a plausible trend line rather than
an error, which is why they are worth pinning individually.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift.analyze import (  # noqa: E402
    MIN_PERIOD_SAMPLES,
    SEVERE_MISMATCH_LOAD,
    AssayResult,
    PeriodSummary,
    evaluate,
    period_of,
    summarise,
    trend,
)
from assaydrift.match import Hit  # noqa: E402
from assaydrift.ncbi import Record  # noqa: E402
from assaydrift.primers import by_name  # noqa: E402


def hit(mismatches=0, ambiguous=0, three_prime=0, found=True, role="forward"):
    return Hit(
        position=10,
        mismatches=mismatches,
        ambiguous=ambiguous,
        three_prime_mismatches=three_prime,
        found=found,
        role=role,
    )


class TestPeriodParsing:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("2021-03-14", "2021-Q1"),
            ("2021-04", "2021-Q2"),
            ("2020-12-31", "2020-Q4"),
            ("2022-07-01", "2022-Q3"),
            ("14-Mar-2021", "2021-Q1"),
        ],
    )
    def test_quarters(self, raw, expected):
        assert period_of(raw) == expected

    def test_year_only_has_no_quarter(self):
        # A bare year cannot be placed in a quarter, and guessing Q1 would pile
        # every undated-to-the-month record onto January.
        assert period_of("2021") is None
        assert period_of("2021", granularity="year") == "2021"

    @pytest.mark.parametrize("raw", ["", "unknown", "not a date", "missing"])
    def test_unusable_dates_are_dropped_not_defaulted(self, raw):
        """The rule that matters most.

        Dating an undated record to today moves old sequences into the current
        period - exactly where they do the most damage to a trend.
        """
        assert period_of(raw) is None

    def test_a_nonsense_month_does_not_produce_a_quarter(self):
        assert period_of("2021-19-01") is None

    def test_month_granularity(self):
        assert period_of("2021-03-14", granularity="month") == "2021-03"


class TestUsability:
    def test_a_result_with_unknown_bases_is_not_usable(self):
        r = AssayResult("X", "A", "2021-Q1", {"f": hit(), "p": hit(ambiguous=2)})
        assert not r.usable

    def test_a_clean_result_is_usable(self):
        r = AssayResult("X", "A", "2021-Q1", {"f": hit(), "p": hit()})
        assert r.usable

    def test_an_empty_result_is_not_usable(self):
        assert not AssayResult("X", "A", "2021-Q1", {}).usable


class TestPerfectMatch:
    """The sensitive measure: does every oligo still match exactly?"""

    def test_an_exact_match_is_perfect(self):
        assert AssayResult("X", "A", "p", {"f": hit()}).perfect()

    def test_any_mismatch_ends_a_perfect_match(self):
        assert not AssayResult("X", "A", "p", {"f": hit(mismatches=1)}).perfect()

    def test_a_missing_binding_site_is_not_a_perfect_match(self):
        r = AssayResult("X", "A", "p", {"f": hit(found=False)})
        assert r.any_oligo_lost
        assert not r.perfect()

    def test_an_assay_with_a_known_baseline_mismatch_is_not_drifted_on_day_one(self):
        """Charite RdRp shipped with one mismatch to SARS-CoV-2.

        Scoring it against a baseline of zero would report it as broken for
        every sequence ever collected, including the 2019 reference.
        """
        r = AssayResult("X", "Charite RdRp", "p", {"r": hit(mismatches=1)})
        assert not r.perfect(baseline=0)
        assert r.perfect(baseline=1)


class TestLikelyFailing:
    """The consequential measure, and a much higher bar than perfect()."""

    def test_one_mismatch_away_from_the_three_prime_end_is_not_failure(self):
        """The correction that mattered most.

        Requiring a perfect match to call an assay working reported CDC N1 and
        Charite E at 0% through 2022-2024, while laboratories were still
        running both successfully. A single mismatch mid-oligo barely moves the
        melting temperature.
        """
        r = AssayResult("X", "A", "p", {"f": hit(mismatches=1, three_prime=0)})
        assert not r.perfect()
        assert not r.likely_failing()

    def test_a_three_prime_mismatch_on_a_primer_is_failure(self):
        r = AssayResult("X", "A", "p", {"f": hit(mismatches=1, three_prime=1, role="forward")})
        assert r.likely_failing()

    def test_a_three_prime_mismatch_on_a_PROBE_is_not_failure(self):
        """A hydrolysis probe is never extended.

        Its 3' end carries the quencher and is chemically blocked, so the
        extension-blocking logic does not apply. Counting probe 3' mismatches as
        fatal is what first reported CDC N1 as dead in 2022 - the Omicron
        mutation under its probe happens to sit two bases from the probe's 3'
        end, which is meaningless for a probe.
        """
        r = AssayResult("X", "A", "p", {"p": hit(mismatches=1, three_prime=1, role="probe")})
        assert not r.likely_failing()
        assert r.three_prime_mismatches == 0

    def test_a_missing_binding_site_is_failure(self):
        assert AssayResult("X", "A", "p", {"f": hit(found=False)}).likely_failing()

    def test_a_heavy_mismatch_load_is_failure(self):
        r = AssayResult("X", "A", "p", {"f": hit(mismatches=SEVERE_MISMATCH_LOAD)})
        assert r.likely_failing()

    def test_the_baseline_is_subtracted_before_the_load_threshold(self):
        r = AssayResult("X", "A", "p", {"f": hit(mismatches=SEVERE_MISMATCH_LOAD)})
        assert not r.likely_failing(baseline=1)


class TestSummaries:
    def test_unknown_bases_are_excluded_rather_than_counted_as_failures(self):
        results = [
            AssayResult("a", "A", "2021-Q1", {"f": hit()}),
            AssayResult("b", "A", "2021-Q1", {"f": hit(ambiguous=3)}),
        ]
        s = summarise(results)["2021-Q1"]
        assert s.total == 2
        assert s.usable == 1
        assert s.perfect == 1
        assert s.perfect_rate == 1.0, "a gappy genome must not read as a failed assay"
        assert s.excluded_rate == 0.5

    def test_a_period_with_nothing_usable_reports_none_not_zero(self):
        """0.0 would be plotted as total assay failure; None says 'no data'."""
        results = [AssayResult("a", "A", "2021-Q1", {"f": hit(ambiguous=1)})]
        s = summarise(results)["2021-Q1"]
        assert s.usable == 0
        assert s.perfect_rate is None
        assert s.failing_rate is None

    def test_periods_come_back_in_order(self):
        results = [
            AssayResult("a", "A", "2022-Q1", {"f": hit()}),
            AssayResult("b", "A", "2020-Q3", {"f": hit()}),
            AssayResult("c", "A", "2021-Q2", {"f": hit()}),
        ]
        assert list(summarise(results)) == ["2020-Q3", "2021-Q2", "2022-Q1"]

    def test_a_small_period_is_marked_unreliable(self):
        small = PeriodSummary("2021-Q1", total=3, usable=3, perfect=3)
        big = PeriodSummary("2021-Q2", total=999, usable=MIN_PERIOD_SAMPLES, perfect=1)
        assert not small.reliable
        assert big.reliable


class TestTrend:
    def _summaries(self, pairs):
        return {
            period: PeriodSummary(period, total=n, usable=n, perfect=int(round(rate * n)))
            for period, rate, n in pairs
        }

    def test_a_falling_rate_is_reported_as_a_negative_change(self):
        t = trend(self._summaries([("2020-Q1", 1.0, 100), ("2024-Q1", 0.5, 100)]))
        assert t["enough_to_say"]
        assert t["change"] == pytest.approx(-0.5)

    def test_one_period_cannot_show_a_slope(self):
        """A single point is not a trend, however many sequences are behind it."""
        t = trend(self._summaries([("2020-Q1", 1.0, 5000)]))
        assert not t["enough_to_say"]

    def test_small_periods_do_not_count_towards_a_trend(self):
        t = trend(self._summaries([("2020-Q1", 1.0, 3), ("2024-Q1", 0.2, 4)]))
        assert not t["enough_to_say"], "two noisy points are not evidence of drift"

    def test_the_worst_period_is_reported_separately_from_the_last(self):
        """An assay that dipped and recovered is a different story from one
        that is failing now, and a first-to-last delta hides the difference."""
        t = trend(
            self._summaries(
                [
                    ("2020-Q1", 1.0, 100),
                    ("2022-Q1", 0.3, 100),
                    ("2024-Q1", 0.95, 100),
                ]
            )
        )
        assert t["worst_period"] == "2022-Q1"
        assert t["worst_rate"] == pytest.approx(0.3)
        assert t["change"] == pytest.approx(-0.05, abs=0.01)


class TestEvaluate:
    def test_a_record_without_a_date_is_skipped(self):
        record = Record(accession="X", organism="o", sequence="ACGT" * 50, collected="")
        assert evaluate(record, by_name("CDC N1")) is None

    def test_a_dated_record_produces_one_hit_per_oligo(self):
        record = Record(accession="X", organism="o", sequence="ACGT" * 50, collected="2021-05-02")
        result = evaluate(record, by_name("CDC N1"))
        assert result is not None
        assert result.period == "2021-Q2"
        assert set(result.hits) == {"N1-F", "N1-R", "N1-P"}


class TestTrendEndpoints:
    """An endpoint is a stronger claim than a table row, and takes a stronger floor.

    `reliable` counts sequences on purpose - gating the per-period table on independence
    would leave 2 of 24 quarters and delete the series rather than qualify it - and
    `effective_n` is published on every row for the reader to weigh. `trend` then read the
    first and last reliable period whatever their independence, and on the
    291,827-sequence corpus the last quarter is partial: 67-77 usable sequences from
    submitter groups so lopsided that Kish's effective n is 1.5.

    Three of the five published headlines were that quarter rather than the assay. CDC N2
    reads 56.6% perfect in 2026-Q1, 47.3% in Q2, 80.1% in Q3 and then 98.5% on an
    effective n of 1.6 - and 98.5% is what the summary quoted, so an assay that is down
    19 points was published as down 0.9.
    """

    def _period(self, name, rate, n, groups, lopsided=False):
        """A period of `n` usable sequences across `groups` submitters.

        Evenly by default, which makes Kish's effective n equal to `groups` - the shape
        of the corpus's large quarters, whose effective n runs 35 to 245. With
        `lopsided`, one submission holds all but one sequence per other group, which is
        the shape of the partial final quarter: 68 sequences, 7 groups, effective n 1.6.

        The first version of this fixture was lopsided for every period, so effective n
        came out near 1 everywhere, every period failed the floor, and the fallback made
        the two tests pass for the wrong reason - a fixture that cannot express the
        condition under test is a test of the fallback.
        """
        summary = PeriodSummary(name, total=n, usable=n, perfect=int(round(rate * n)))
        if groups == 1:
            summary.submitters = {"one": n}
        elif lopsided:
            summary.submitters = {f"s{i}": 1 for i in range(groups - 1)}
            summary.submitters["big"] = n - (groups - 1)
        else:
            share, extra = divmod(n, groups)
            summary.submitters = {
                f"s{i}": share + (1 if i < extra else 0) for i in range(groups)
            }
        assert sum(summary.submitters.values()) == n
        return summary

    def _summaries(self, rows):
        return {r[0]: self._period(*r) for r in rows}

    def test_the_fixture_produces_the_independence_it_claims(self):
        """Prove the fixture can express both conditions before trusting a verdict."""
        even = self._period("2026-Q3", 0.8, 1860, 60)
        assert even.effective_n == pytest.approx(60, abs=1)
        lopsided = self._period("2026-Q4", 0.985, 68, 7, lopsided=True)
        assert lopsided.effective_n == pytest.approx(1.2, abs=0.5)
        assert lopsided.reliable, "it passes the sequence floor; that is the point"

    def test_a_lopsided_final_period_is_not_the_endpoint(self):
        t = trend(self._summaries([
            ("2020-Q1", 0.994, 9000, 200),
            ("2026-Q3", 0.801, 1860, 60),
            ("2026-Q4", 0.985, 67, 7, True),
        ]))
        assert t["last_period"] == "2026-Q3"
        assert t["change"] == pytest.approx(-0.193, abs=0.002), (
            "the headline must be the assay, not a quarter worth 1.6 observations"
        )
        assert t["endpoints_meet_the_independence_floor"]

    def test_the_skipped_periods_are_named_with_their_rates(self):
        """Withdrawn silently is worse than qualified out loud."""
        t = trend(self._summaries([
            ("2020-Q1", 0.994, 9000, 200),
            ("2026-Q3", 0.801, 1860, 60),
            ("2026-Q4", 0.985, 67, 7, True),
        ]))
        skipped = t["periods_after_the_endpoint"]
        assert [s["period"] for s in skipped] == ["2026-Q4"]
        assert skipped[0]["perfect_rate"] == pytest.approx(0.985, abs=0.01)
        assert skipped[0]["effective_n"] < t["independence_floor"]

    def test_an_independent_final_period_is_kept(self):
        """The floor must not simply discard the end of every series."""
        t = trend(self._summaries([
            ("2020-Q1", 1.0, 9000, 200),
            ("2026-Q4", 0.5, 2000, 400),
        ]))
        assert t["last_period"] == "2026-Q4"
        assert t["periods_after_the_endpoint"] == []
        assert t["change"] == pytest.approx(-0.5)

    def test_a_corpus_where_no_period_is_independent_still_gets_a_trend(self):
        """Qualified, not refused: such a corpus still has a story to tell."""
        t = trend(self._summaries([
            ("2020-Q1", 1.0, 100, 1),
            ("2024-Q1", 0.4, 100, 1),
        ]))
        assert t["enough_to_say"]
        assert t["change"] == pytest.approx(-0.6)
        assert not t["endpoints_meet_the_independence_floor"], (
            "the caveat is the whole point of the fallback"
        )

    def test_the_worst_period_may_be_one_that_cannot_be_an_endpoint(self):
        """"The worst it ever got" is a different question from "where did it end".

        A quarter too lopsided to anchor a trend is still evidence that a rate that low
        occurred, so the worst-period figures stay over every reliable period.
        """
        t = trend(self._summaries([
            ("2020-Q1", 1.0, 9000, 200),
            ("2023-Q3", 0.051, 469, 5, True),
            ("2026-Q3", 0.868, 1891, 60),
        ]))
        assert t["worst_period"] == "2023-Q3"
        assert t["worst_rate"] == pytest.approx(0.051, abs=0.01)
        assert t["last_period"] == "2026-Q3"
