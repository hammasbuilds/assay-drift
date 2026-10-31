"""A fragment that never reached the target region is not an assay with no binding site.

The 291,827-genome corpus contains LR877184: 2,859 bases, collected 2020-04-26 in
Stockholm, with no unknown bases at all. SARS-CoV-2 is about 29,900 bases, so this
record is a partial deposit. Every oligo of all five assays is absent from it, which
made `any_oligo_lost` true five times and put one record into clause 1 of
`likely_failing` - the only clause with no judgement call in it, the one the paper calls
"the only clause a bench would not dispute". The strongest claim the tool can make was
being produced by a short record.

`region_absent` separates the two facts, and the length condition is half of it on
purpose: a full-length genome that has genuinely lost an entire target region must still
be reported as a lost oligo, because that is the finding this rule must not be able to
hide. These tests pin both halves and both directions.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.analyze import (  # noqa: E402
    MIN_WHOLE_GENOME,
    AssayResult,
    summarise,
)
from assaydrift.match import Hit  # noqa: E402

LR877184_LENGTH = 2859  # the record that prompted this rule
WHOLE_GENOME = 29903  # Wuhan-Hu-1


def hit(found=True, mismatches=0, three_prime=0, role="forward"):
    return Hit(
        position=10 if found else -1,
        mismatches=mismatches,
        ambiguous=0,
        three_prime_mismatches=three_prime,
        found=found,
        role=role,
    )


def result(length, **hits) -> AssayResult:
    r = AssayResult(accession="X", assay="A", period="2020-Q2", sequence_length=length)
    r.hits = dict(hits)
    return r


ALL_ABSENT = {
    "f": hit(found=False),
    "r": hit(found=False, role="reverse"),
    "p": hit(found=False, role="probe"),
}


class TestTheFragmentCase:
    def test_a_fragment_with_every_oligo_absent_is_region_absent(self):
        assert result(LR877184_LENGTH, **ALL_ABSENT).region_absent

    def test_a_whole_genome_with_every_oligo_absent_is_not(self):
        """A real total dropout stays a finding, not an exclusion."""
        r = result(WHOLE_GENOME, **ALL_ABSENT)
        assert not r.region_absent
        assert r.any_oligo_lost

    def test_the_boundary_is_the_published_constant(self):
        assert result(MIN_WHOLE_GENOME - 1, **ALL_ABSENT).region_absent
        assert not result(MIN_WHOLE_GENOME, **ALL_ABSENT).region_absent

    def test_one_oligo_absent_in_a_fragment_is_still_a_loss(self):
        """Partial absence is interpretable and must not be excused by length.

        A fragment that aligns two oligos did reach the region, so a third one missing
        is evidence about that oligo's site rather than about the record.
        """
        r = result(LR877184_LENGTH, f=hit(found=False), r=hit(role="reverse"), p=hit(role="probe"))
        assert not r.region_absent
        assert r.any_oligo_lost

    def test_an_unrecorded_length_is_not_treated_as_a_fragment(self):
        """Zero means "not recorded". Guessing "fragment" there would silently delete
        every lost-oligo observation from any caller that does not set the field."""
        assert not result(0, **ALL_ABSENT).region_absent

    def test_no_hits_at_all_is_not_region_absent(self):
        assert not result(LR877184_LENGTH).region_absent


class TestSummariseDropsThemWithoutHidingThem:
    def test_the_fragment_is_counted_and_excluded(self):
        s = summarise([result(LR877184_LENGTH, **ALL_ABSENT)])["2020-Q2"]
        assert s.region_absent == 1
        assert s.lost_oligo == 0, "a fragment must not reach clause 1"
        assert s.usable == 0
        assert s.total == 1, "dropping it from `total` would hide it from excluded_rate"

    def test_a_whole_genome_dropout_still_reaches_clause_one(self):
        s = summarise([result(WHOLE_GENOME, **ALL_ABSENT)])["2020-Q2"]
        assert s.region_absent == 0
        assert s.lost_oligo == 1

    def test_ordinary_records_are_untouched(self):
        s = summarise([result(WHOLE_GENOME, f=hit(), r=hit(role="reverse"), p=hit(role="probe"))])[
            "2020-Q2"
        ]
        assert s.region_absent == 0
        assert s.usable == 1
        assert s.perfect == 1


class TestTheCliUsesTheSamePredicate:
    """cli.py had its own length guard (`too_short`) and this record walked past it.

    2,859 bases is longer than every oligo, so `too_short` returns 0 and the CLI
    reported the fragment as `likely_failing` with three oligos lost - the corpus-path
    bug on the single-sequence path. The CLI now calls `region_absent` rather than
    restating the rule.
    """

    def test_a_fragment_is_unknown_not_likely_failing(self):
        from assaydrift.cli import check  # noqa: PLC0415
        from assaydrift.primers import ASSAYS  # noqa: PLC0415

        assay = ASSAYS[0]
        fragment = "ACGT" * 700  # 2,800 bases of nothing the assay binds
        rows = check(assay, [("fragment", fragment)])
        assert rows[0]["verdict"] == "unknown", rows[0]
        assert all(not o["found"] for o in rows[0]["oligos"].values())

    def test_too_short_still_covers_what_it_covered(self):
        from assaydrift.cli import too_short  # noqa: PLC0415
        from assaydrift.primers import ASSAYS  # noqa: PLC0415

        assay = ASSAYS[0]
        assert too_short("ACGT", assay) == len(assay.oligos)
        assert too_short("ACGT" * 700, assay) == 0


def test_clause_one_is_zero_in_the_published_results():
    """The paper's sentence about clause 1, read off whichever artefacts are present.

    It was true of the 2,765-genome corpus and false of the 291,827-genome one by
    exactly one record per assay - the same record, which is the shape of a bad deposit
    rather than of five findings.
    """
    checked = []
    for name in ("results/drift.json", "results/drift_big.json"):
        path = ROOT / name
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        if "region_absent" not in next(
            iter(next(v for k, v in data.items() if k != "_run")["periods"].values())
        ):
            continue  # predates the rule; regenerate it
        checked.append(name)
        for assay, body in data.items():
            if assay == "_run":
                continue
            lost = {p: s["lost_oligo"] for p, s in body["periods"].items() if s["lost_oligo"]}
            assert not lost, (
                f"{name}: {assay} reports a lost binding site in {lost}. The paper says "
                f"clause 1 fires for no sequence; check whether these are whole genomes "
                f"or fragments that `region_absent` should have excluded."
            )
    if not checked:
        pytest.skip("no regenerated results carry region_absent yet")
