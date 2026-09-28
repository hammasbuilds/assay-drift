"""Tests for the named-mutation report.

Every bug this file guards against produces a plausible mutation name rather
than an error, which is why they are worth pinning individually.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift.mutations import observe, peak, rank, reference_sites  # noqa: E402
from assaydrift.ncbi import Record  # noqa: E402
from assaydrift.primers import by_name  # noqa: E402


def _record(accession: str, sequence: str, collected: str) -> Record:
    return Record(
        accession=accession,
        organism="Severe acute respiratory syndrome coronavirus 2",
        sequence=sequence,
        collected=collected,
        country="",
        length=len(sequence),
    )


def _genome_with(position: int, base: str) -> str:
    """The reference genome with one base swapped, 1-based `position`."""
    from assaydrift.reference import genome

    ref = genome()
    return ref[: position - 1] + base + ref[position:]


class TestReferenceSites:
    def test_forward_oligo_numbers_5_to_3(self):
        assay = by_name("CDC N1")
        sites = reference_sites([assay])
        forward = sites["N1-F"]
        assert [s.base_in_oligo for s in forward] == list(range(1, len(forward) + 1))
        assert forward[0].bases_from_3_prime_end == len(forward) - 1
        assert forward[-1].bases_from_3_prime_end == 0

    def test_reverse_oligo_numbers_are_flipped_from_the_search_direction(self):
        assay = by_name("CDC N1")
        sites = reference_sites([assay])
        reverse = sites["N1-R"]
        # Site j=0 walks the search string (the reverse complement, 5' of the
        # genome plus strand); it is the oligo's own 3' end.
        assert reverse[0].bases_from_3_prime_end == 0
        assert reverse[0].base_in_oligo == len(reverse)
        assert reverse[-1].bases_from_3_prime_end == len(reverse) - 1
        assert reverse[-1].base_in_oligo == 1


class TestObserve:
    def test_a_sample_matching_the_reference_is_not_a_mutation(self):
        assay = by_name("CDC N1")
        ref = _genome_with(0, "A")  # no-op copy of the reference
        record = _record("REF", ref, "2021")
        observations, denominators = observe([record], [assay], "year")
        assert observations == {}
        assert denominators[f"{assay.name}|2021"] == 1

    def test_a_real_substitution_is_named_reference_position_alt(self):
        assay = by_name("CDC N1")
        sites = reference_sites([assay])
        # N1-F's first base, forward strand: genome_position is 1-based.
        site = sites["N1-F"][0]
        mutated = _genome_with(site.genome_position, "A" if site.reference_base != "A" else "C")
        record = _record("MUT", mutated, "2022")
        observations, denominators = observe([record], [assay], "year")
        assert len(observations) == 1
        (obs,) = observations.values()
        assert obs.name == f"{site.reference_base}{site.genome_position}{obs.alt_base}"
        assert obs.by_period["2022"] == 1

    def test_the_charite_rdrp_designed_mismatch_is_never_reported_as_a_mutation(self):
        """RdRp-R was designed with one mismatch against the 2019 reference
        (documented in primers.KNOWN_REFERENCE_MISMATCHES). A sample that
        still reads as the reference at that position has not drifted - it is
        the assay's own primer that departs from the reference, not the virus.
        Counting it as a "mutation" would report it as carried by ~100% of
        every sample, in every period, forever.
        """
        assay = by_name("Charite RdRp")
        ref = _genome_with(0, "A")  # unmodified reference sequence
        records = [_record(f"R{i}", ref, "2020") for i in range(5)]
        observations, _denominators = observe(records, [assay], "year")
        assert observations == {}, [o.name for o in observations.values()]

    def test_an_unknown_base_is_not_counted_as_a_substitution(self):
        assay = by_name("CDC N1")
        sites = reference_sites([assay])
        site = sites["N1-F"][0]
        genome = _genome_with(site.genome_position, "N")
        record = _record("AMBIG", genome, "2023")
        observations, denominators = observe([record], [assay], "year")
        assert observations == {}
        # the record is still usable overall (only this one base is N)
        assert denominators.get(f"{assay.name}|2023", 0) in (0, 1)


class TestRankAndPeak:
    def test_rank_orders_most_frequent_first(self):
        assay = by_name("CDC N1")
        sites = reference_sites([assay])
        site_a = sites["N1-F"][0]
        site_b = sites["N1-F"][1]
        common = _genome_with(site_a.genome_position, "A" if site_a.reference_base != "A" else "C")
        rare = _genome_with(site_b.genome_position, "A" if site_b.reference_base != "A" else "C")
        records = (
            [_record(f"C{i}", common, "2021") for i in range(8)]
            + [_record(f"R{i}", rare, "2021") for i in range(2)]
            + [_record(f"P{i}", common[:], "2021") for i in range(0)]
        )
        observations, denominators = observe(records, [assay], "year")
        ranked = rank(observations, denominators)
        assert len(ranked) == 2
        (top_obs, top_rate), (second_obs, second_rate) = ranked
        assert top_rate >= second_rate

    def test_peak_reports_the_period_with_the_highest_rate_not_the_most_recent(self):
        assay = by_name("CDC N1")
        sites = reference_sites([assay])
        site = sites["N1-F"][0]
        mutated = _genome_with(site.genome_position, "A" if site.reference_base != "A" else "C")
        ref = _genome_with(0, "A")
        records = (
            [_record(f"E{i}", mutated, "2020") for i in range(9)]
            + [_record(f"F{i}", ref, "2020") for i in range(1)]
            + [_record(f"L{i}", mutated, "2023") for i in range(1)]
            + [_record(f"M{i}", ref, "2023") for i in range(9)]
        )
        observations, denominators = observe(records, [assay], "year")
        (obs,) = observations.values()
        period, rate = peak(obs, denominators)
        assert period == "2020"
        assert rate == 0.9
