"""Every published primer must match the genome it was designed against.

This is the test that decides whether any number this project reports means
anything. A single mistyped base in `primers.py` would show up as drift that is
not there - the assay would appear to mismatch every sequence ever collected,
including the reference, and the write-up would read as a discovery.

So: each oligo is aligned against NC_045512.2, the first SARS-CoV-2 genome
deposited, and must match perfectly. The one real exception is listed by name in
`primers.KNOWN_REFERENCE_MISMATCHES` with the reason, which is the difference
between documenting a property of an assay and quietly loosening a threshold.

The reference is vendored at `tests/data/NC_045512.2.json.gz` (9 KB) so the
suite never touches the network. A test that needs NCBI to be reachable is a
test that fails on a plane, in CI behind a proxy, and on the day NCBI has an
outage - and the reference genome has not changed since 2020.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift import primers  # noqa: E402
from assaydrift.match import IUPAC, find_oligo  # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "data" / "NC_045512.2.json.gz"

# Pinned so a corrupted or swapped fixture fails loudly rather than silently
# changing what "matches the reference" means.
REFERENCE_SHA256 = "7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"


@pytest.fixture(scope="module")
def reference() -> str:
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        record = json.load(handle)
    sequence = record["sequence"]
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    assert digest == REFERENCE_SHA256, "the vendored reference genome is not the expected one"
    assert len(sequence) == 29903, "SARS-CoV-2 reference is 29,903 bases"
    return sequence


class TestOligoDefinitions:
    def test_every_oligo_is_valid_iupac(self):
        for _, oligo in primers.all_oligos():
            assert set(oligo.sequence) <= set(IUPAC), oligo.name

    def test_rejects_a_non_nucleotide_sequence(self):
        with pytest.raises(ValueError, match="IUPAC"):
            primers.Oligo("bad", "ACGTXQ", "forward")

    def test_every_assay_has_a_forward_a_reverse_and_a_probe(self):
        for assay in primers.ASSAYS:
            roles = sorted(o.role for o in assay.oligos)
            assert roles == ["forward", "probe", "reverse"], assay.name

    def test_oligo_lengths_are_plausible(self):
        # Real PCR primers are roughly 18-30 bases. Anything far outside that
        # is a transcription error, not an assay design.
        for _, oligo in primers.all_oligos():
            assert 15 <= len(oligo) <= 35, f"{oligo.name} is {len(oligo)} bases"

    def test_assay_names_are_unique(self):
        names = [a.name for a in primers.ASSAYS]
        assert len(names) == len(set(names))

    def test_by_name_raises_clearly_for_an_unknown_assay(self):
        with pytest.raises(KeyError, match="no assay named"):
            primers.by_name("not a real assay")


class TestAgainstTheReference:
    """The load-bearing test."""

    def test_every_oligo_binds_the_reference_genome(self, reference):
        for assay, oligo in primers.all_oligos():
            hit = find_oligo(reference, oligo.sequence, oligo.role)
            assert hit.found, (
                f"{oligo.name} ({assay.name}) has no binding site in the reference "
                f"genome at all - almost certainly a typo in primers.py"
            )

    def test_every_oligo_matches_the_reference_exactly(self, reference):
        """...apart from the one documented exception."""
        for assay, oligo in primers.all_oligos():
            hit = find_oligo(reference, oligo.sequence, oligo.role)
            expected = primers.expected_reference_mismatches(oligo.name)
            assert hit.mismatches == expected, (
                f"{oligo.name} ({assay.name}) has {hit.mismatches} mismatches against "
                f"the reference, expected {expected}. Either the sequence in "
                f"primers.py is wrong, or this is a real property that belongs in "
                f"KNOWN_REFERENCE_MISMATCHES with a reason."
            )

    def test_no_unknown_bases_in_the_reference_binding_sites(self, reference):
        for _, oligo in primers.all_oligos():
            assert find_oligo(reference, oligo.sequence, oligo.role).ambiguous == 0

    def test_the_known_exception_is_the_charite_rdrp_reverse_primer(self, reference):
        """Pinned in detail, because it is the one place the rule is relaxed.

        RdRp_SARSr-R was designed for SARS-related coronaviruses generally, from
        sequences that predate SARS-CoV-2. It carries S (G or C) where SARS-CoV-2
        has T, so this assay has never matched its target perfectly.
        """
        assert list(primers.KNOWN_REFERENCE_MISMATCHES) == ["RdRp-R"]
        known = primers.KNOWN_REFERENCE_MISMATCHES["RdRp-R"]
        oligo = next(o for _, o in primers.all_oligos() if o.name == "RdRp-R")
        hit = find_oligo(reference, oligo.sequence, oligo.role)
        assert hit.mismatches == known["mismatches"] == 1
        # Mid-primer, not at the 3' terminus, which is why the assay worked.
        assert hit.three_prime_mismatches == 0
        assert known["bases_from_3_prime_end"] == 14

    def test_binding_sites_land_in_the_right_genes(self, reference):
        """A sanity check on the matcher, not only on the sequences.

        SARS-CoV-2 gene coordinates in NC_045512.2 are fixed and published. If
        the N assays were landing in ORF1ab, the matcher would be finding
        coincidental alignments in a 30,000 base genome.
        """
        genes = {
            "N (nucleocapsid)": (28274, 29533),
            "E (envelope)": (26245, 26472),
            "RdRp (ORF1ab)": (266, 21555),
        }
        for assay, oligo in primers.all_oligos():
            low, high = genes[assay.target_gene]
            hit = find_oligo(reference, oligo.sequence, oligo.role)
            assert low <= hit.position <= high, (
                f"{oligo.name} binds at {hit.position}, outside {assay.target_gene} ({low}-{high})"
            )
