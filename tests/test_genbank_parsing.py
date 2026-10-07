"""Parsing a GenBank record must not change the length of the sequence.

`_one_record` filtered the ORIGIN block with `re.sub(r"[^acgtnACGTN]", "", ...)`, which
DELETED every other IUPAC code instead of keeping it - and deleting a base renumbers the
rest of the genome. 339 of the 2,781 records in the cache behind the published result
carry at least one (1,262 bases in all: Y 541, R 403, K 103, M 89, W 45, H 38, S 31, V 9,
B 2, D 1), so 12% of the corpus was stored a base or more short, with every coordinate
past the first ambiguous position off by one.

Rebuilding the corpus from the same cached responses restored 1,261 bases across 338
records - up to 85 in one - and moved the published result: CDC N1's 2024-Q4 likely-
failing count fell from 13 to 3 (9.2% to 2.3%) as those sequences moved into the excluded
bucket, where an unreadable base belongs. The error ran in the direction of overstating
drift.
"""

from __future__ import annotations

from assaydrift.ncbi import IUPAC_CODES, _one_record

HEADER = """LOCUS       TEST                   60 bp    DNA     linear   VRL 01-JAN-2024
ACCESSION   TEST00001
  ORGANISM  Severe acute respiratory syndrome coronavirus 2
     source          1..60
                     /collection_date="2024-01-15"
                     /geo_loc_name="USA: Testland"
ORIGIN
"""


def _record(sequence: str):
    """A minimal GenBank flatfile whose ORIGIN carries `sequence`, numbered as GenBank
    numbers it - which is the formatting the parser has to strip without touching the
    bases."""
    lines = []
    for start in range(0, len(sequence), 60):
        chunk = sequence[start : start + 60]
        groups = " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10))
        lines.append(f"{start + 1:>9} {groups}")
    return _one_record(HEADER + "\n".join(lines) + "\n//\n")


def test_an_ambiguity_code_is_kept_not_dropped() -> None:
    seq = "ACGTRYSWKMBDHVNACGT" * 3
    rec = _record(seq)
    assert rec is not None
    assert rec.sequence == seq, "the parser altered the sequence"
    assert rec.length == len(seq)


def test_every_iupac_code_survives_a_round_trip() -> None:
    for code in sorted(IUPAC_CODES):
        seq = "ACGT" * 4 + code + "ACGT" * 4
        rec = _record(seq)
        assert rec is not None, code
        assert rec.sequence == seq, f"{code} did not survive parsing"


def test_position_is_preserved_when_an_ambiguous_base_comes_first() -> None:
    """The defect's actual consequence: not a wrong base, a wrong coordinate.

    One deleted base at the front moves every later base one position earlier, so an
    oligo's reported position, the base index inside it, and the mutation calls read off
    those coordinates are all wrong - silently, and only for the 12% of records that
    carry an ambiguity code.
    """
    marker = "GGGGGGGGGG"
    seq = "R" + "ACGT" * 10 + marker
    rec = _record(seq)
    assert rec is not None
    assert rec.sequence.index(marker) == seq.index(marker)
    assert len(rec.sequence) == len(seq)


def test_numbering_and_spacing_are_still_stripped() -> None:
    """Keeping every letter must not mean keeping the flatfile's own formatting."""
    seq = "ACGT" * 30
    rec = _record(seq)
    assert rec is not None
    assert rec.sequence == seq
    assert not any(c.isdigit() or c.isspace() for c in rec.sequence)


def test_a_letter_that_is_not_a_nucleotide_becomes_n_rather_than_vanishing() -> None:
    """Unknown, counted, and still occupying its position.

    Dropping it would renumber the genome again; guessing a base would invent data. `N`
    is the one answer that is both honest and length-preserving, and the rates already
    exclude oligos with an `N` under them.
    """
    seq = "ACGT" * 5 + "Z" + "ACGT" * 5
    rec = _record(seq)
    assert rec is not None
    assert len(rec.sequence) == len(seq)
    assert rec.sequence[20] == "N"
    assert rec.sequence.replace("N", "Z", 1) == seq
