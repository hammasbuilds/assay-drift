"""Per-period assay performance: does this test still detect what is circulating?

A rate on its own answers the wrong question. "97% of sequences match" tells you
nothing about whether an assay is fine or failing, because an assay that was
100% in 2020 and 80% today averages to something reassuring. Drift is a *slope*,
so everything here is grouped by when the sample was collected.

Three rules the grouping follows, each of which exists because breaking it
produces a confident wrong trend:

**A sequence with no collection date is dropped, never dated to today.**
Defaulting an undated record to the present moment moves old sequences into the
current period, which is precisely where they would corrupt the answer.

**Periods below a minimum sample size are reported but not plotted as trend.**
A quarter with four sequences can show any rate at all. The count is printed
beside every rate so a reader can discount it themselves.

**An oligo whose binding site contains `N` is excluded from the rate, not
counted as a failure.** See `match.py` - this is the difference between
measuring viral drift and measuring sequencing quality.
"""

from __future__ import annotations

import re
import statistics
from collections import defaultdict
from dataclasses import dataclass, field

from .match import Hit, find_oligo
from .ncbi import Record
from .primers import Assay

# Below this, a period's rate is too noisy to read as a trend point.
MIN_PERIOD_SAMPLES = 25

# Total mismatches across an assay's three oligos, above the assay's own
# baseline, at which binding is assumed to fail even without a 3' hit. Three is
# a judgement call rather than a measured threshold: it is the point at which a
# 20-25mer has lost enough complementarity that the melting temperature falls
# below a typical 60C annealing step. It is stated here so it can be argued
# with rather than buried.
SEVERE_MISMATCH_LOAD = 3


def submitter_block(accession: str, country: str) -> str:
    """Group records that plausibly arrived in one submission.

    GenBank hands out consecutive accessions to a single submission, so the accession
    prefix plus the thousand-block, together with the reported country, approximates
    "same laboratory, same batch". It is a proxy - the records carry no submitter field -
    and it is deliberately coarse, because the point is not to identify laboratories but
    to stop a quarter of 159 sequences from 4 groups being read as 159 observations.
    """
    match = re.match(r"([A-Za-z]+)(\d+)", accession or "")
    stem = f"{match.group(1)}{int(match.group(2)) // 1000}" if match else (accession or "?")
    return f"{country or '?'}|{stem}"


def period_of(collected: str, granularity: str = "quarter") -> str | None:
    """Bucket a GenBank collection_date. Returns None when it cannot be read.

    GenBank dates are free text and arrive in several shapes: `2021-03-14`,
    `2021-03`, `2021`, `14-Mar-2021`, and occasionally a range. Only the year
    and month are needed here, and anything that yields neither is dropped.
    """
    if not collected:
        return None
    year = re.search(r"(19|20)\d{2}", collected)
    if not year:
        return None
    y = year.group(0)
    if granularity == "year":
        return y

    month = None
    iso = re.search(rf"{y}-(\d{{2}})", collected)
    if iso:
        month = int(iso.group(1))
    else:
        names = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
        found = re.search(r"([A-Za-z]{3})", collected)
        if found and found.group(1).lower() in names:
            month = names.index(found.group(1).lower()) + 1
    if month is None or not 1 <= month <= 12:
        return y if granularity == "year" else None
    if granularity == "month":
        return f"{y}-{month:02d}"
    return f"{y}-Q{(month - 1) // 3 + 1}"


@dataclass
class AssayResult:
    """How one assay performed on one sequence."""

    accession: str
    assay: str
    period: str
    hits: dict[str, Hit] = field(default_factory=dict)
    submitter: str = ""
    """A proxy for "same submitting laboratory, same batch".

    GenBank records carry no submitter field, but consecutive accessions are handed out to
    one submission, so the accession prefix plus thousand-block together with the reported
    country groups records that arrived together. Used only to count how many independent
    groups a period's rate rests on - see PeriodSummary.effective_n.
    """

    @property
    def usable(self) -> bool:
        """Every oligo located, with no unknown bases under any of them."""
        return bool(self.hits) and all(h.usable for h in self.hits.values())

    @property
    def total_mismatches(self) -> int:
        return sum(h.mismatches for h in self.hits.values())

    @property
    def three_prime_mismatches(self) -> int:
        """3' mismatches that can stop extension - primers only.

        A hydrolysis probe's 3' end is blocked by its quencher and is never
        extended, so a mismatch there is an ordinary binding penalty rather than
        a stop signal. Counting probe 3' mismatches here reported CDC N1 as
        having failed outright in 2022, because the Omicron mutation under its
        probe happens to sit two bases from that probe's 3' end.
        """
        return sum(h.blocks_extension for h in self.hits.values())

    @property
    def any_oligo_lost(self) -> bool:
        """An oligo that cannot be found at all - the assay has no binding site."""
        return any(not h.found for h in self.hits.values())

    def perfect(self, baseline: int = 0) -> bool:
        """Every oligo matches exactly.

        This is the clean drift signal: it falls as soon as the target changes
        under any oligo, whether or not the change is enough to break the test.

        `baseline` is the assay's own mismatch count against the reference
        genome, so an assay that shipped with a known mismatch is not scored as
        drifted on day one. See primers.KNOWN_REFERENCE_MISMATCHES.
        """
        return not self.any_oligo_lost and self.total_mismatches <= baseline

    def likely_failing(self, baseline: int = 0) -> bool:
        """Whether the assay would plausibly stop detecting the sample.

        Deliberately a much higher bar than `perfect`. A single mismatch near
        the 5' end of a 26-base primer barely moves the melting temperature, and
        an assay carrying one keeps working - the Charite E assay stayed in
        clinical use through Omicron with exactly that. Treating any mismatch as
        failure reported these assays at 0% while laboratories were still
        running them successfully, which is an overclaim in the alarming
        direction.

        What actually stops a reaction:
          - an oligo with no binding site left at all
          - a mismatch in the last five bases of a PRIMER, where extension starts
          - a heavy mismatch load across the assay, which drops binding outright
        """
        if self.any_oligo_lost:
            return True
        if self.three_prime_mismatches > 0:
            return True
        return self.total_mismatches - baseline >= SEVERE_MISMATCH_LOAD


def evaluate(record: Record, assay: Assay, granularity: str = "quarter") -> AssayResult | None:
    """Run one assay against one sequence. None when the record has no date."""
    period = period_of(record.collected, granularity)
    if period is None:
        return None
    result = AssayResult(
        accession=record.accession,
        assay=assay.name,
        period=period,
        submitter=submitter_block(record.accession, getattr(record, "country", "")),
    )
    for oligo in assay.oligos:
        result.hits[oligo.name] = find_oligo(record.sequence, oligo.sequence, oligo.role)
    return result


@dataclass
class PeriodSummary:
    period: str
    total: int = 0  # sequences with this assay evaluated
    usable: int = 0  # ... and with no N under any oligo
    perfect: int = 0  # ... and matching every oligo exactly
    failing: int = 0  # ... and plausibly no longer amplifying
    lost_oligo: int = 0  # an oligo with no binding site at all
    three_prime: int = 0  # a mismatch in the last 5 bases of a PRIMER
    mismatch_counts: list[int] = field(default_factory=list)
    submitters: dict[str, int] = field(default_factory=dict)
    """Usable records per submitter group, for effective_n.

    `usable` counts sequences, and sequences from one submission batch are not
    independent evidence about a population: a 159-sequence quarter here can be 97%
    one batch.
    """

    @property
    def perfect_rate(self) -> float | None:
        """Share of usable sequences the assay matches exactly.

        The sensitive measure: it moves as soon as the target changes under any
        oligo, whether or not the change is enough to break the test.

        None rather than 0.0 when nothing was usable: an assay that could not be
        measured is not an assay that failed, and 0.0 would be plotted as
        catastrophic drift.
        """
        return self.perfect / self.usable if self.usable else None

    @property
    def failing_rate(self) -> float | None:
        """Share that would plausibly no longer be detected.

        The consequential measure, and a much higher bar - see
        AssayResult.likely_failing. Reported alongside perfect_rate because the
        gap between them is the interesting part: an assay can drift a long way
        before it stops working.
        """
        return self.failing / self.usable if self.usable else None

    @property
    def excluded_rate(self) -> float | None:
        """Share dropped for unknown bases. High values mean poor sequencing,
        not drift - and a period where most data is excluded cannot support a
        claim in either direction."""
        return (self.total - self.usable) / self.total if self.total else None

    @property
    def median_mismatches(self) -> float | None:
        return statistics.median(self.mismatch_counts) if self.mismatch_counts else None

    @property
    def clusters(self) -> int:
        """How many distinct submitter groups the usable records came from."""
        return len(self.submitters)

    @property
    def effective_n(self) -> float | None:
        """Kish effective sample size over the submitter groups.

        n^2 / sum(group sizes squared): the number of INDEPENDENT observations this
        period's rate is worth. A period of 159 sequences that is 97% one submission is
        worth about one. Reported beside every rate because `usable` counts sequences,
        and sequences from one batch are not independent evidence about a population.
        """
        sizes = list(self.submitters.values())
        total = sum(sizes)
        if not total:
            return None
        return total * total / sum(s * s for s in sizes)

    @property
    def reliable(self) -> bool:
        # Counts sequences, deliberately unchanged: gating this on effective_n would
        # leave 2 of 24 quarters and delete the series rather than qualify it. The
        # independence figure is published alongside instead, so a reader can see that a
        # "reliable" quarter can rest on one submission.
        return self.usable >= MIN_PERIOD_SAMPLES


def summarise(results: list[AssayResult], baseline: int = 0) -> dict[str, PeriodSummary]:
    """Collapse per-sequence results into per-period rates."""
    out: dict[str, PeriodSummary] = defaultdict(lambda: PeriodSummary(period=""))
    for result in results:
        summary = out[result.period]
        summary.period = result.period
        summary.total += 1
        if result.any_oligo_lost:
            summary.lost_oligo += 1
        if not result.usable:
            continue
        summary.usable += 1
        summary.submitters[result.submitter or result.accession] = (
            summary.submitters.get(result.submitter or result.accession, 0) + 1
        )
        summary.mismatch_counts.append(result.total_mismatches)
        if result.three_prime_mismatches > 0:
            summary.three_prime += 1
        if result.perfect(baseline):
            summary.perfect += 1
        if result.likely_failing(baseline):
            summary.failing += 1
    return dict(sorted(out.items()))


def trend(summaries: dict[str, PeriodSummary]) -> dict:
    """First and last reliable periods, and the change between them.

    Deliberately not a regression line. With a handful of periods, a fitted
    slope invites more confidence than the data supports; the honest summary is
    where it started, where it ended, and how many sequences are behind each.
    """
    reliable = [s for s in summaries.values() if s.reliable and s.perfect_rate is not None]
    if len(reliable) < 2:
        return {"periods": len(reliable), "enough_to_say": False}
    first, last = reliable[0], reliable[-1]
    return {
        "periods": len(reliable),
        "enough_to_say": True,
        "first_period": first.period,
        "first_rate": first.perfect_rate,
        "first_n": first.usable,
        "last_period": last.period,
        "last_rate": last.perfect_rate,
        "last_n": last.usable,
        "change": last.perfect_rate - first.perfect_rate,
        "worst_period": min(reliable, key=lambda s: s.perfect_rate).period,
        "worst_rate": min(s.perfect_rate for s in reliable),
        # The consequential measure, alongside the sensitive one. The gap
        # between them is the point: an assay drifts long before it breaks.
        "first_failing_rate": first.failing_rate,
        "last_failing_rate": last.failing_rate,
        "worst_failing_rate": max(s.failing_rate for s in reliable),
        "worst_failing_period": max(reliable, key=lambda s: s.failing_rate).period,
    }
