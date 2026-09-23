"""Published PCR diagnostic assays, as oligonucleotide sequences.

These are the actual primers and probes that clinical labs ran on millions of
samples. Each one is a short stretch of DNA that must match the target genome
for the test to detect it; when the organism mutates underneath an assay, the
test starts returning false negatives while reporting nothing unusual.

**Every sequence here is checked against the reference genome it was designed
from.** A primer that does not match the 2019 reference with zero mismatches is
a typo in this file, not a discovery about drift, and `tests/test_primers.py`
refuses to let one through. Getting a base wrong would manufacture drift that
does not exist - the single most likely way for this project to publish a
confident wrong answer.

Sources
-------
CDC 2019-nCoV Real-Time RT-PCR Diagnostic Panel, Instructions for Use
    CDC-006-00019, Revision 06 (2020). Primers N1, N2 and the retired N3.
Corman VM et al. "Detection of 2019 novel coronavirus (2019-nCoV) by real-time
    RT-PCR." Euro Surveill. 2020;25(3):2000045. Assays E_Sarbeco and RdRp_SARSr.

The N3 assay is included **because it was withdrawn**. CDC removed it from the
panel in 2020, so it is the one assay here already known to have had a problem,
which makes it the closest thing available to a positive control.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# The genome every assay below was designed against: Wuhan-Hu-1, the first
# SARS-CoV-2 sequence deposited. Drift is measured as distance from this.
REFERENCE = "NC_045512.2"
REFERENCE_ORGANISM = "Severe acute respiratory syndrome coronavirus 2"


@dataclass(frozen=True)
class Oligo:
    """One primer or probe."""

    name: str
    sequence: str
    role: str  # "forward", "reverse" or "probe"

    def __post_init__(self) -> None:
        bad = set(self.sequence.upper()) - set("ACGTRYSWKMBDHVN")
        if bad:
            raise ValueError(f"{self.name}: not IUPAC nucleotide codes: {sorted(bad)}")

    def __len__(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class Assay:
    """A complete diagnostic assay: two primers and a probe."""

    name: str
    target_gene: str
    source: str
    oligos: list[Oligo] = field(default_factory=list)
    retired: bool = False
    note: str = ""

    def __iter__(self):
        return iter(self.oligos)


ASSAYS: list[Assay] = [
    Assay(
        name="CDC N1",
        target_gene="N (nucleocapsid)",
        source="CDC-006-00019 rev.06",
        oligos=[
            Oligo("N1-F", "GACCCCAAAATCAGCGAAAT", "forward"),
            Oligo("N1-R", "TCTGGTTACTGCCAGTTGAATCTG", "reverse"),
            Oligo("N1-P", "ACCCCGCATTACGTTTGGTGGACC", "probe"),
        ],
    ),
    Assay(
        name="CDC N2",
        target_gene="N (nucleocapsid)",
        source="CDC-006-00019 rev.06",
        oligos=[
            Oligo("N2-F", "TTACAAACATTGGCCGCAAA", "forward"),
            Oligo("N2-R", "GCGCGACATTCCGAAGAA", "reverse"),
            Oligo("N2-P", "ACAATTTGCCCCCAGCGCTTCAG", "probe"),
        ],
    ),
    Assay(
        name="CDC N3",
        target_gene="N (nucleocapsid)",
        source="CDC-006-00019 rev.03",
        retired=True,
        note=(
            "Withdrawn by CDC in 2020. Included as the one assay already known "
            "to have had a problem - the closest thing to a positive control "
            "this study can have."
        ),
        oligos=[
            Oligo("N3-F", "GGGAGCCTTGAATACACCAAAA", "forward"),
            Oligo("N3-R", "TGTAGCACGATTGCAGCATTG", "reverse"),
            Oligo("N3-P", "AYCACATTGGCACCCGCAATCCTG", "probe"),
        ],
    ),
    Assay(
        name="Charite E",
        target_gene="E (envelope)",
        source="Corman et al. 2020, Euro Surveill 25(3)",
        oligos=[
            Oligo("E-F", "ACAGGTACGTTAATAGTTAATAGCGT", "forward"),
            Oligo("E-R", "ATATTGCAGCAGTACGCACACA", "reverse"),
            Oligo("E-P", "ACACTAGCCATCCTTACTGCGCTTCG", "probe"),
        ],
    ),
    Assay(
        name="Charite RdRp",
        target_gene="RdRp (ORF1ab)",
        source="Corman et al. 2020, Euro Surveill 25(3)",
        note=(
            "Carries IUPAC ambiguity codes (R, S) by design, because it was "
            "written to detect SARS-related coronaviruses broadly rather than "
            "SARS-CoV-2 alone. An ambiguous base matches several nucleotides, "
            "so this assay should tolerate drift better than the others. "
            "It is also the one assay that never matched SARS-CoV-2 perfectly: "
            "see KNOWN_REFERENCE_MISMATCHES."
        ),
        oligos=[
            Oligo("RdRp-F", "GTGARATGGTCATGTGTGGCGG", "forward"),
            Oligo("RdRp-R", "CARATGTTAAASACACTATTAGCATA", "reverse"),
            Oligo("RdRp-P", "CAGGTGGAACCTCATCAGGAGATGC", "probe"),
        ],
    ),
]


# Oligos that do NOT match the 2019 reference perfectly, with the reason.
#
# The rule in this file is that a primer failing against NC_045512.2 is a typo.
# There is exactly one real exception, and listing it here is the difference
# between a documented property of an assay and a loosened test.
#
# RdRp_SARSr-R was designed as a pan-Sarbecovirus primer from SARS-CoV (2003)
# and bat SARS-related coronaviruses, before SARS-CoV-2 existed. It carries `S`
# (G or C) at a position where SARS-CoV-2 has `T`, so the Charite RdRp assay has
# had one mismatch to its target since the day it was published.
#
# It sits 14 bases from the 3' end - the middle of the primer - which is why the
# assay worked anyway. The same mismatch within five bases of the 3' terminus
# would likely have stopped extension. This is the clearest illustration in the
# whole project of why position matters more than count.
KNOWN_REFERENCE_MISMATCHES: dict[str, dict] = {
    "RdRp-R": {
        "mismatches": 1,
        "oligo_base": "S",
        "reference_base": "T",
        "bases_from_3_prime_end": 14,
        "why": (
            "designed for SARS-related coronaviruses generally, from sequences "
            "that predate SARS-CoV-2"
        ),
    },
}


def expected_reference_mismatches(oligo_name: str) -> int:
    """How many mismatches against NC_045512.2 are known and accounted for."""
    return KNOWN_REFERENCE_MISMATCHES.get(oligo_name, {}).get("mismatches", 0)


def by_name(name: str) -> Assay:
    for assay in ASSAYS:
        if assay.name == name:
            return assay
    raise KeyError(f"no assay named {name!r}; have {[a.name for a in ASSAYS]}")


def all_oligos() -> list[tuple[Assay, Oligo]]:
    return [(assay, oligo) for assay in ASSAYS for oligo in assay.oligos]
