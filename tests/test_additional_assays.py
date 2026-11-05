"""The second cohort's sequences, against the reference and against their own geometry.

Only one of these eight assays comes from the article that published it. The other seven
were transcribed from a compilation, and a compilation can carry transcription errors - so
each is admitted on three independent structural checks rather than on trust:

  1. every oligo binds NC_045512.2 with zero mismatches, since all of them were designed
     against that genome;
  2. each binds in the region its own name claims;
  3. the three oligos form an amplicon - forward before reverse, probe strictly inside.

The checks are not decoration. The Japan NIID N2 set failed the first and is absent: its
reverse primer as compiled reads TGGCAGCTGTGTAGGTCAAC, one mismatch from the reference's
TGGCACCTGTGTAGGTCAAC at base 6. Editing that base would have produced a sequence passing
every check here which might not be the assay anyone ran, so neither version was admitted.
`test_the_niid_set_is_absent_and_stays_absent` pins that decision, because the tempting
fix is exactly the one that would make it silently wrong.
"""

from __future__ import annotations

import gzip
import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.match import find_oligo  # noqa: E402
from assaydrift.primers import ADDITIONAL_ASSAYS, ASSAYS, all_oligos  # noqa: E402

REFERENCE_FILE = ROOT / "src" / "assaydrift" / "data" / "NC_045512.2.json.gz"

# assay -> the genomic window its name or target implies, as (low, high).
EXPECTED_REGION = {
    "Charite N": (28_200, 29_600),
    "HKU ORF1b": (18_000, 19_500),
    "HKU N": (28_200, 29_600),
    "China CDC ORF1ab": (13_000, 14_000),
    "China CDC N": (28_200, 29_600),
    "Pasteur IP2": (12_400, 13_000),
    "Pasteur IP4": (13_900, 14_400),
    "Thailand N": (28_200, 29_600),
}


@pytest.fixture(scope="module")
def reference() -> str:
    if not REFERENCE_FILE.is_file():
        pytest.skip("the vendored reference genome is not present")
    payload = json.loads(gzip.decompress(REFERENCE_FILE.read_bytes()).decode())
    sequence = (
        payload
        if isinstance(payload, str)
        else payload.get("sequence") or next(v for v in payload.values() if isinstance(v, str))
    )
    assert len(sequence) == 29_903, len(sequence)
    return sequence


def _hits(reference: str, assay) -> dict[str, tuple[int, int]]:
    out = {}
    for oligo in assay.oligos:
        hit = find_oligo(reference, oligo.sequence, oligo.role)
        assert hit.found, f"{oligo.name} ({assay.name}) does not bind the reference"
        out[oligo.role] = (hit.position, hit.position + len(oligo.sequence))
    return out


@pytest.mark.parametrize("assay", ADDITIONAL_ASSAYS, ids=lambda a: a.name)
def test_every_oligo_binds_the_reference_with_no_mismatches(assay, reference):
    """Check one. A compiled sequence that does not match is a transcription error."""
    for oligo in assay.oligos:
        hit = find_oligo(reference, oligo.sequence, oligo.role)
        assert hit.found, f"{oligo.name} has no binding site in the reference"
        assert hit.mismatches == 0, (
            f"{oligo.name} ({assay.name}) has {hit.mismatches} mismatches against "
            f"NC_045512.2; a compiled sequence that does not match the genome it was "
            f"designed against is a transcription error, not a finding"
        )
        assert hit.ambiguous == 0, f"{oligo.name} has {hit.ambiguous} ambiguous positions"


@pytest.mark.parametrize("assay", ADDITIONAL_ASSAYS, ids=lambda a: a.name)
def test_each_assay_binds_where_its_name_says(assay, reference):
    """Check two. A correct sequence in the wrong place is a different assay."""
    low, high = EXPECTED_REGION[assay.name]
    for role, (start, end) in _hits(reference, assay).items():
        assert low <= start and end <= high, (assay.name, role, start, end, (low, high))


@pytest.mark.parametrize("assay", ADDITIONAL_ASSAYS, ids=lambda a: a.name)
def test_the_three_oligos_form_an_amplicon(assay, reference):
    """Check three. Forward before reverse, probe strictly between them."""
    hits = _hits(reference, assay)
    assert set(hits) == {"forward", "reverse", "probe"}, (assay.name, sorted(hits))
    (f0, f1), (r0, r1), (p0, p1) = hits["forward"], hits["reverse"], hits["probe"]
    assert f0 < r0, f"{assay.name}: the forward primer does not precede the reverse"
    amplicon = r1 - f0
    assert 50 <= amplicon <= 250, f"{assay.name}: implausible {amplicon}-base amplicon"
    assert f1 <= p0 and p1 <= r0, f"{assay.name}: the probe is not inside the amplicon"


def test_the_niid_set_is_absent_and_stays_absent():
    """The one candidate that failed, and the edit that would hide the failure.

    Neither the compiled sequence nor the reference-matching one could be attributed to
    NIID, so the assay is not here. The failure mode this guards is subtle: changing the
    one base makes the sequence pass every other test in this file.
    """
    names = {assay.name for assay in ASSAYS + ADDITIONAL_ASSAYS}
    assert not any("NIID" in name or "Japan" in name for name in names), names
    sequences = {oligo.sequence for _assay, oligo in all_oligos(ASSAYS + ADDITIONAL_ASSAYS)}
    assert "TGGCAGCTGTGTAGGTCAAC" not in sequences, "the compiled NIID reverse primer"
    assert "TGGCACCTGTGTAGGTCAAC" not in sequences, (
        "the reference-matching NIID reverse primer - admitting this would be inventing "
        "a published sequence from the genome it is meant to be tested against"
    )


def test_the_reference_version_of_the_niid_primer_really_would_have_passed(reference):
    """A guard is only worth having if the thing it forbids would otherwise slip through."""
    compiled = find_oligo(reference, "TGGCAGCTGTGTAGGTCAAC", "reverse")
    corrected = find_oligo(reference, "TGGCACCTGTGTAGGTCAAC", "reverse")
    assert compiled.found and compiled.mismatches == 1, compiled.mismatches
    assert corrected.found and corrected.mismatches == 0, corrected.mismatches


def test_the_cohorts_are_separate_and_the_primary_one_is_unchanged():
    """Every published per-assay number is about `ASSAYS`; the cohort must not touch it."""
    assert [a.name for a in ASSAYS] == [
        "CDC N1",
        "CDC N2",
        "CDC N3",
        "Charite E",
        "Charite RdRp",
    ]
    assert len(all_oligos()) == 15, "the primary cohort's oligo count moved"
    assert len(ADDITIONAL_ASSAYS) == 8
    assert len(all_oligos(ADDITIONAL_ASSAYS)) == 24
    overlap = {a.name for a in ASSAYS} & {a.name for a in ADDITIONAL_ASSAYS}
    assert not overlap, overlap


def test_every_additional_assay_records_where_its_sequences_came_from():
    """Seven of eight are second-hand, and the file has to keep saying which."""
    for assay in ADDITIONAL_ASSAYS:
        assert assay.source, assay.name
    primary = [a for a in ADDITIONAL_ASSAYS if "read in full" in a.source]
    compiled = [a for a in ADDITIONAL_ASSAYS if "compilation" in a.source]
    assert len(primary) == 1 and primary[0].name == "Charite N", [a.name for a in primary]
    assert len(compiled) == 7, [a.name for a in compiled]
