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

Two cohorts
-----------
`ASSAYS` is those five, and every published per-assay number in this repository
is about them. `ADDITIONAL_ASSAYS` is eight further deployed protocols - Charite
N, HKU's ORF1b and N, China CDC's ORF1ab and N, Institut Pasteur's IP2 and IP4,
and Thailand NIH's N - added so that section 4's claim about a *family* of rules
can be tested somewhere other than the five assays it was derived from. They are
a separate list precisely so that adding them moves nothing already published;
`research/cohort_generality.py` scores both with the same rule.

Only Charite N was read from the article that published it. The other seven come
from a compilation and are admitted on three structural checks - exact binding to
NC_045512.2, binding in the region their names claim, and a valid amplicon with
the probe inside - which one candidate failed. See the comment above
`ADDITIONAL_ASSAYS` for that case and why it was dropped rather than repaired.
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


# A second cohort, for asking whether section 4's findings are about these five assays or
# about the rule. The five above are the ones this study was built on: the CDC panel, the
# two Charite assays, and N3 as a positive control. Eight assays is a small base for a
# claim about how a *family* of rules behaves, so these are the other widely deployed
# protocols, kept in their own list so that every published per-assay number above is
# unaffected by their arrival.
#
# **Provenance, and why it is stated separately for each.** Only Charite N comes from the
# article that published it: Corman et al., Euro Surveill 2020;25(3):2000045, Table 1,
# read from the Europe PMC full text. The other seven were transcribed from a compilation
# (Dahmani, Azzoune and Boudjema, Mol Biol Res Commun 2024;13(4):223-234, Table 1), because
# their original protocols are PDFs this project could not reach - the NIID manual URL is
# dead, the Chu et al. full text is not open access, and the J-STAGE copy of the NIID paper
# is DRM-encrypted.
#
# A compilation can carry transcription errors, so every candidate was put through three
# independent structural checks before being admitted here, and one failed:
#
#   1. **It must bind NC_045512.2 with zero mismatches.** Every assay here was designed
#      against that genome, so anything else is a transcription error.
#   2. **It must bind in the region its own name claims** - IP2's "12669" label against a
#      binding site at 12,689, N-gene assays inside the N gene, and so on.
#   3. **Its oligos must form an amplicon**: forward before reverse, 57 to 132 bases apart,
#      with the probe strictly between them.
#
# The Japan NIID N2 set failed the first check and is **deliberately absent**. Its reverse
# primer as compiled, TGGCAGCTGTGTAGGTCAAC, has one mismatch to the reference at base 6;
# the reference there reads TGGCACCTGTGTAGGTCAAC. One of those is what NIID published and
# this project could not establish which, so neither is here. Editing the base to match the
# reference would have produced a sequence that passes every check in this file and might
# not be the assay anyone ran - which is the failure mode all three checks exist to catch,
# arrived at from the other direction.
ADDITIONAL_ASSAYS: list[Assay] = [
    Assay(
        name="Charite N",
        target_gene="N (nucleocapsid)",
        source="Corman et al. 2020, Euro Surveill 25(3):2000045, Table 1 (read in full)",
        oligos=[
            Oligo("CharN-F", "CACATTGGCACCCGCAATC", "forward"),
            Oligo("CharN-R", "GAGGAACGAGAAGAGGCTTG", "reverse"),
            Oligo("CharN-P", "ACTTCCTCAAGGAACAACATTGCCA", "probe"),
        ],
    ),
    Assay(
        name="HKU ORF1b",
        target_gene="ORF1b-nsp14",
        source="Chu et al. 2020, Clin Chem 66(4):549-555, via compilation",
        oligos=[
            Oligo("HKU1b-F", "TGGGGYTTTACRGGTAACCT", "forward"),
            Oligo("HKU1b-R", "AACRCGCTTAACAAAGCACTC", "reverse"),
            Oligo("HKU1b-P", "TAGTTGTGATGCWATCATGACTAG", "probe"),
        ],
    ),
    Assay(
        name="HKU N",
        target_gene="N (nucleocapsid)",
        source="Chu et al. 2020, Clin Chem 66(4):549-555, via compilation",
        oligos=[
            Oligo("HKUN-F", "TAATCAGACAAGGAACTGATTA", "forward"),
            Oligo("HKUN-R", "CGAAGGTGTGACTTCCATG", "reverse"),
            Oligo("HKUN-P", "CCGCAAATTGCACAATTTGC", "probe"),
        ],
    ),
    Assay(
        name="China CDC ORF1ab",
        target_gene="ORF1ab",
        source="China CDC protocol, via compilation",
        oligos=[
            Oligo("CN1ab-F", "CCCTGTGGGTTTTACACTTAA", "forward"),
            Oligo("CN1ab-R", "ACGATTGTGCATCAGCTGA", "reverse"),
            Oligo("CN1ab-P", "CCGTCTGCGGTATGTGGAAAGGTTATGG", "probe"),
        ],
    ),
    Assay(
        name="China CDC N",
        target_gene="N (nucleocapsid)",
        source="China CDC protocol, via compilation",
        note=(
            "The assay [@genes2021] reports falling furthest - 63.89% inclusivity by "
            "January 2021, from 86.03% in April 2020 - which makes it this cohort's "
            "closest analogue to CDC N3 above."
        ),
        oligos=[
            Oligo("CNN-F", "GGGGAACTTCTCCTGCTAGAAT", "forward"),
            Oligo("CNN-R", "CAGACATTTTGCTCTCAAGCTG", "reverse"),
            Oligo("CNN-P", "TTGCTGCTGCTTGACAGATT", "probe"),
        ],
    ),
    Assay(
        name="Pasteur IP2",
        target_gene="ORF1ab (nsp9/nsp10)",
        source="Institut Pasteur Paris protocol, via compilation",
        oligos=[
            Oligo("IP2-F", "ATGAGCTTAGTCCTGTTG", "forward"),
            Oligo("IP2-R", "CTCCCTTTGTTGTGTTGT", "reverse"),
            Oligo("IP2-P", "AGATGTCTTGTGCTGCCGGTA", "probe"),
        ],
    ),
    Assay(
        name="Pasteur IP4",
        target_gene="ORF1ab (RdRp)",
        source="Institut Pasteur Paris protocol, via compilation",
        oligos=[
            Oligo("IP4-F", "GGTAACTGGTATGATTTCG", "forward"),
            Oligo("IP4-R", "CTGGTCAAGGTTAATATAGG", "reverse"),
            Oligo("IP4-P", "TCATACAAACCACGCCAGG", "probe"),
        ],
    ),
    Assay(
        name="Thailand N",
        target_gene="N (nucleocapsid)",
        source="Thailand NIH protocol, via compilation",
        note="A 57-base amplicon, the shortest here, with a 16-base probe inside it.",
        oligos=[
            Oligo("THN-F", "CGTTTGGTGGACCCTCAGAT", "forward"),
            Oligo("THN-R", "CCCCACTGCGTTCTCCATT", "reverse"),
            Oligo("THN-P", "CAACTGGCAGTAACCA", "probe"),
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


def all_oligos(assays: list[Assay] | None = None) -> list[tuple[Assay, Oligo]]:
    """Every (assay, oligo) pair. Defaults to the primary cohort, as it always did."""
    return [
        (assay, oligo) for assay in (ASSAYS if assays is None else assays) for oligo in assay.oligos
    ]
