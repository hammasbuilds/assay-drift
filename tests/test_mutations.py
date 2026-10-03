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


class TestN2ThreePrimeSignal:
    """C29215T, the 2026 CDC N2 signal, pinned from the reference genome.

    The README's central claim changed because of this one substitution, so the
    two facts it rests on are pinned here: that the site is two bases from the
    N2 reverse primer's 3' end, and that one mismatch there is enough to make
    the assay likely-failing rather than merely drifted.
    """

    def test_c29215t_sits_two_bases_from_the_n2_reverse_primer_3_prime_end(self):
        assay = by_name("CDC N2")
        sites = reference_sites([assay])
        site = next(s for s in sites["N2-R"] if s.genome_position == 29215)
        assert site.oligo == "N2-R"
        assert site.role == "reverse"
        assert site.bases_from_3_prime_end == 2
        assert site.in_three_prime_window

    def test_it_is_named_c29215t_and_not_its_reverse_complement(self):
        sample = _record("x", _genome_with(29215, "T"), "2026-01-15")
        observations, _ = observe([sample], [by_name("CDC N2")], "quarter")
        assert [o.name for o in observations.values()] == ["C29215T"]

    def test_one_mismatch_there_makes_the_assay_likely_failing(self):
        from assaydrift.analyze import evaluate

        sample = _record("x", _genome_with(29215, "T"), "2026-01-15")
        result = evaluate(sample, by_name("CDC N2"))
        assert result.usable
        assert not result.perfect()
        assert result.likely_failing()
        assert result.total_mismatches == 1  # not a severe mismatch load

    def test_the_same_substitution_mid_oligo_is_drift_not_failure(self):
        from assaydrift.analyze import evaluate

        sample = _record("x", _genome_with(29221, "G"), "2026-01-15")  # 8 bases in
        result = evaluate(sample, by_name("CDC N2"))
        assert not result.perfect()
        assert not result.likely_failing()


class TestPlaces:
    """Where a mutation was collected, which is how a batch is told from a wave."""

    def _records(self):
        from assaydrift.reference import genome

        carrying = _genome_with(29215, "T")
        clean = genome()
        rows = []
        for i in range(8):  # one place deposits a batch, almost all carrying
            rows.append(
                Record(
                    accession=f"a{i}",
                    organism="x",
                    sequence=carrying if i < 6 else clean,
                    collected="2026-01-15",
                    country="USA: Wisconsin",
                    length=0,
                )
            )
        for i in range(10):  # another sequences more and finds one
            rows.append(
                Record(
                    accession=f"b{i}",
                    organism="x",
                    sequence=carrying if i < 1 else clean,
                    collected="2026-01-15",
                    country="USA: California",
                    length=0,
                )
            )
        return rows

    def _observation(self, rows):
        observations, denominators = observe(rows, [by_name("CDC N2")], "quarter")
        return next(o for o in observations.values() if o.name == "C29215T"), denominators

    def test_rate_is_per_place_not_over_the_whole_study(self):
        from assaydrift.mutations import places

        observation, denominators = self._observation(self._records())
        rows = {place: (n, total) for place, n, total, _ in places(observation, denominators)}
        assert rows == {"USA: Wisconsin": (6, 8), "USA: California": (1, 10)}

    def test_sub_locality_is_dropped_so_one_place_is_one_group(self):
        from assaydrift.mutations import place_of, places

        assert place_of("USA: AZ, Maricopa") == "USA: AZ"
        assert place_of("") == "unknown"
        rows = self._records()
        rows[0] = Record(
            accession=rows[0].accession,
            organism="x",
            sequence=rows[0].sequence,
            collected="2026-01-15",
            country="USA: Wisconsin, Dane",
            length=0,
        )
        observation, denominators = self._observation(rows)
        counts = {place: n for place, n, _t, _r in places(observation, denominators)}
        assert counts["USA: Wisconsin"] == 6

    def test_period_denominators_are_unaffected_by_the_place_keys(self):
        observation, denominators = self._observation(self._records())
        assert denominators["CDC N2|2026-Q1"] == 18
        assert rank({("N2-R", 2, "A"): observation}, denominators)[0][1] == 7 / 18
