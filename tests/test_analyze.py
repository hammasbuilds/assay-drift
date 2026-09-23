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


def hit(mismatches=0, ambiguous=0, three_prime=0, found=True):
    return Hit(
        position=10,
        mismatches=mismatches,
        ambiguous=ambiguous,
        three_prime_mismatches=three_prime,
        found=found,
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


class TestIntactness:
    def test_a_perfect_match_is_intact(self):
        assert AssayResult("X", "A", "p", {"f": hit()}).intact()

    def test_a_three_prime_mismatch_breaks_the_assay(self):
        """One mismatch under the polymerase's start point is not like one in
        the middle, so it is not allowed to average away."""
        r = AssayResult("X", "A", "p", {"f": hit(mismatches=1, three_prime=1)})
        assert not r.intact()

    def test_a_missing_binding_site_breaks_the_assay(self):
        r = AssayResult("X", "A", "p", {"f": hit(found=False)})
        assert r.any_oligo_lost
        assert not r.intact()

    def test_an_assay_with_a_known_baseline_mismatch_is_not_drifted_on_day_one(self):
        """Charite RdRp shipped with one mismatch to SARS-CoV-2.

        Scoring it against a baseline of zero would report it as broken for
        every sequence ever collected, including the 2019 reference.
        """
        r = AssayResult("X", "Charite RdRp", "p", {"r": hit(mismatches=1)})
        assert not r.intact(baseline=0)
        assert r.intact(baseline=1)


class TestSummaries:
    def test_unknown_bases_are_excluded_rather_than_counted_as_failures(self):
        results = [
            AssayResult("a", "A", "2021-Q1", {"f": hit()}),
            AssayResult("b", "A", "2021-Q1", {"f": hit(ambiguous=3)}),
        ]
        s = summarise(results)["2021-Q1"]
        assert s.total == 2
        assert s.usable == 1
        assert s.intact == 1
        assert s.intact_rate == 1.0, "a gappy genome must not read as a failed assay"
        assert s.excluded_rate == 0.5

    def test_a_period_with_nothing_usable_reports_none_not_zero(self):
        """0.0 would be plotted as total assay failure; None says 'no data'."""
        results = [AssayResult("a", "A", "2021-Q1", {"f": hit(ambiguous=1)})]
        s = summarise(results)["2021-Q1"]
        assert s.usable == 0
        assert s.intact_rate is None

    def test_periods_come_back_in_order(self):
        results = [
            AssayResult("a", "A", "2022-Q1", {"f": hit()}),
            AssayResult("b", "A", "2020-Q3", {"f": hit()}),
            AssayResult("c", "A", "2021-Q2", {"f": hit()}),
        ]
        assert list(summarise(results)) == ["2020-Q3", "2021-Q2", "2022-Q1"]

    def test_a_small_period_is_marked_unreliable(self):
        small = PeriodSummary("2021-Q1", total=3, usable=3, intact=3)
        big = PeriodSummary("2021-Q2", total=999, usable=MIN_PERIOD_SAMPLES, intact=1)
        assert not small.reliable
        assert big.reliable


class TestTrend:
    def _summaries(self, pairs):
        return {
            period: PeriodSummary(period, total=n, usable=n, intact=int(round(rate * n)))
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
