"""The seven gaps suite-auditor proved in this package's own suite.

Each one is a mutant the suite did not notice, with the call that separates it from the
original - all seven on inputs the tests really use, so none is an artefact of a
made-up argument. They are grouped by the function they were found in, with the mutation
each one defeats named, because a test whose purpose is not written down is a test that
gets deleted in the next refactor for looking redundant.

Found by running `suite-auditor audit .` over this repository: 152 mutants scored, 79.6%
kill rate, 81.2% of functions reached by a test.
"""

from __future__ import annotations

import pytest

from assaydrift.analyze import period_of, submitter_block
from assaydrift.match import MAX_SEED, MIN_SEED, SEED_GUARANTEE, _seed_length, _seeds


class TestSubmitterBlock:
    """Two mutants survived here, both on `submitter_block("X", "")`.

    An accession that does not match `([A-Za-z]+)(\\d+)` and an empty country together
    exercise both fallbacks at once, and the suite only ever passed well-formed pairs.
    """

    def test_an_unparsable_accession_falls_back_to_the_accession_itself(self) -> None:
        # Defeats `drop_return`: the function returned None instead of the group key, so
        # every record with an odd accession would have shared one None cluster.
        assert submitter_block("X", "") == "?|X"

    def test_an_empty_country_becomes_a_question_mark_not_the_accession(self) -> None:
        # Defeats `boolop`: `country or '?'` mutated to put '?' on both sides, which
        # merges every unparsable accession from every country into one cluster - and
        # `effective_n` is computed from cluster counts, so that inflates independence.
        assert submitter_block("X", "") == "?|X"
        assert submitter_block("X", "USA") == "USA|X"

    def test_a_well_formed_accession_is_grouped_by_its_thousand_block(self) -> None:
        assert submitter_block("MW642250", "USA") == "USA|MW642"
        assert submitter_block("MW642999", "USA") == "USA|MW642"
        assert submitter_block("MW643000", "USA") == "USA|MW643"

    def test_two_countries_never_share_a_block(self) -> None:
        """The whole point of including the country: a block is a laboratory proxy."""
        assert submitter_block("MW642250", "USA") != submitter_block("MW642250", "India")


class TestPeriodOf:
    """One mutant survived on `period_of("14-Mar-2021", "month")`, giving 2021-01.

    The suite covered quarter granularity and ISO dates. Month granularity on a
    day-first, month-name date reached the `names.index(...) + 1` arithmetic, and that
    `+ 1` was never checked - a mutation to it shifts every month-name date by one month,
    which moves sequences between quarters at a quarter boundary.
    """

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("14-Mar-2021", "2021-03"),
            ("01-Jan-2021", "2021-01"),
            ("31-Dec-2021", "2021-12"),
            ("2021-03-14", "2021-03"),
            ("2021-03", "2021-03"),
        ],
    )
    def test_month_granularity_reads_the_month(self, text: str, expected: str) -> None:
        assert period_of(text, "month") == expected

    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("01-Jan-2021", "2021-Q1"),
            ("31-Mar-2021", "2021-Q1"),
            ("01-Apr-2021", "2021-Q2"),
            ("30-Jun-2021", "2021-Q2"),
            ("01-Jul-2021", "2021-Q3"),
            ("01-Oct-2021", "2021-Q4"),
            ("31-Dec-2021", "2021-Q4"),
        ],
    )
    def test_every_quarter_boundary_from_a_month_name(self, text: str, expected: str) -> None:
        """A one-month shift is invisible except at a boundary, so all four are here."""
        assert period_of(text) == expected

    def test_a_month_name_that_is_not_a_month_is_not_a_month(self) -> None:
        assert period_of("14-Foo-2021", "month") is None
        assert period_of("14-Foo-2021", "year") == "2021"

    def test_a_year_with_no_readable_month_is_dropped_except_at_year_granularity(self) -> None:
        assert period_of("2021", "month") is None
        assert period_of("2021") is None
        assert period_of("2021", "year") == "2021"


class TestSeedLength:
    """Three mutants survived on `_seed_length(20)`, giving 7, 5 and 8 instead of 6.

    `(oligo_length - SEED_GUARANTEE) // (SEED_GUARANTEE + 1)` is the whole reason the
    seed search cannot be defeated by two mismatches, and nothing checked the arithmetic.
    A seed one base too long stops tiling the oligo, which is exactly the defect the
    docstring says cost 31.4% of two-mismatch sites in an 18-mer.
    """

    @pytest.mark.parametrize(
        ("oligo_length", "expected"),
        [(18, 5), (20, 6), (22, 6), (24, 7), (26, 8), (30, 8), (12, 5), (5, 5)],
    )
    def test_the_seed_length_for_each_real_oligo_size(
        self, oligo_length: int, expected: int
    ) -> None:
        assert _seed_length(oligo_length) == expected

    def test_the_result_stays_inside_its_bounds(self) -> None:
        for n in range(1, 60):
            assert MIN_SEED <= _seed_length(n) <= MAX_SEED

    def test_the_guarantee_the_formula_exists_for(self) -> None:
        """SEED_GUARANTEE mismatches cannot break every seed of a tiling.

        With k = (n - g) // (g + 1), the oligo holds at least g + 1 non-overlapping
        seeds, so g mismatches cannot touch all of them. This is the property the three
        surviving mutants each broke - 7, 5 and 8 for a 20-mer - and the one worth
        asserting rather than the number.
        """
        for n in range(MIN_SEED * (SEED_GUARANTEE + 1), 40):
            k = _seed_length(n)
            disjoint = n // k
            assert disjoint >= SEED_GUARANTEE + 1, (
                f"a {n}-mer with {k}-base seeds holds only {disjoint} disjoint seeds, "
                f"so {SEED_GUARANTEE} mismatches could break all of them"
            )


class TestSeeds:
    """A `slice_upper` mutant survived on `_seeds("ACGTACGTACGTACGTACGT", None)`.

    It shortened every seed by one base. Each seed was still a valid k-mer and the list
    was still the right length, so only the seed CONTENTS separate the two - which is
    what the suite never looked at.
    """

    def test_the_seeds_are_the_oligos_own_k_mers_at_stride_one(self) -> None:
        oligo = "ACGTACGTACGTACGTACGT"
        k = _seed_length(len(oligo))
        seeds = _seeds(oligo)
        assert seeds[0] == (0, oligo[:k])
        assert all(piece == oligo[start : start + k] for start, piece in seeds)
        assert all(len(piece) == k for _, piece in seeds)
        assert len(seeds) == len(oligo) - k + 1

    def test_the_seeds_tile_the_oligo_with_no_base_left_out(self) -> None:
        """The property the docstring is about: no base inside no seed."""
        oligo = "ACGTACGTACGTACGTAC"
        covered: set[int] = set()
        for start, piece in _seeds(oligo):
            covered.update(range(start, start + len(piece)))
        assert covered == set(range(len(oligo)))

    def test_an_explicit_k_overrides_the_computed_one(self) -> None:
        oligo = "ACGTACGTACGTACGTACGT"
        assert all(len(piece) == 4 for _, piece in _seeds(oligo, 4))

    def test_a_seed_containing_an_ambiguity_code_is_dropped(self) -> None:
        """Only unambiguous k-mers can be searched with `str.find`."""
        oligo = "ACGTACGTRCGTACGTACGT"
        for start, piece in _seeds(oligo):
            assert "R" not in piece
            assert piece == oligo[start : start + len(piece)]
