"""What `check` does with a FASTA that is not plain uppercase DNA.

Found by feeding it the reference genome written as RNA, which is what several pipelines
emit and what every textbook writes for a coronavirus. It scored `likely_failing` against
an assay the sequence matches perfectly: `U` is not an IUPAC *DNA* code, so `compare_base`
fell to its "a gap character, or junk in the record" branch and called every base a
mismatch.

The two input paths disagreed about such letters, too. `ncbi.py` stores anything outside
the IUPAC set as `N` - missing data, excluded and reported. `read_fasta`, the path a
user's own file takes, only uppercased, so the same letter became a definite mismatch and
pushed the verdict towards drift.

Three kinds of character, three different right answers, and conflating any two of them
gets a verdict wrong - which the first version of the fix did, by turning alignment gaps
into `N` and replacing a correct `perfect` with `unknown`.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from assaydrift.cli import InputError, check, normalise, read_fasta, read_primers, read_sequences

HERE = Path(__file__).resolve().parent.parent
EXAMPLES = HERE / "examples"


@pytest.fixture(scope="module")
def assay():
    return read_primers(EXAMPLES / "cdc_n2.tsv")


@pytest.fixture(scope="module")
def reference() -> str:
    return dict(read_fasta(EXAMPLES / "two_genomes.fasta"))["NC_045512.2"]


def write(tmp_path: Path, name: str, sequence: str) -> Path:
    path = tmp_path / f"{name}.fasta"
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")
    return path


def verdicts(assay, path: Path) -> list[str]:
    records, _ = read_sequences(path)
    return [row["verdict"] for row in check(assay, records)]


# --- RNA ----------------------------------------------------------------------------


def test_an_rna_sequence_scores_as_the_dna_it_is(tmp_path, assay, reference):
    """The headline defect: a perfect match reported as the most alarming verdict."""
    path = write(tmp_path, "rna", reference.replace("T", "U"))
    assert verdicts(assay, path) == ["perfect"]


def test_the_rna_rewrite_is_reported_and_counted(tmp_path, assay, reference):
    """A tool that silently edits someone's sequence is worse than either behaviour."""
    path = write(tmp_path, "rna", reference.replace("T", "U"))
    _, notes = read_sequences(path)
    assert notes["rna"] == (reference.count("T"), 0, 0)


def test_u_and_t_are_the_same_base_to_the_scorer(reference):
    sequence, rna, gaps, unknown = normalise(reference.replace("T", "U"))
    assert sequence == reference
    assert (rna, gaps, unknown) == (reference.count("T"), 0, 0)


# --- alignment gaps -----------------------------------------------------------------


def test_stripping_gaps_restores_the_sequence_that_was_aligned(reference):
    """Byte for byte, which is the justification for stripping rather than reading."""
    site = reference.find("TTACAAACATTGGCCGCAAA")
    assert site > 0
    gapped = reference[: site + 5] + "---" + reference[site + 5 :]
    sequence, rna, gaps, unknown = normalise(gapped)
    assert sequence == reference
    assert (rna, gaps, unknown) == (0, 3, 0)


def test_gaps_inside_a_primer_site_do_not_become_unknown_bases(tmp_path, assay, reference):
    """The first version of this fix called them `N` and lost a correct verdict.

    `match.py` aligns an oligo across an indel with an affine-gap aligner, so a gapped
    record already scored `perfect`; turning the gaps into unknown bases replaced that
    with "not evidence either way".
    """
    site = reference.find("TTACAAACATTGGCCGCAAA")
    path = write(tmp_path, "aligned", reference[: site + 5] + "---" + reference[site + 5 :])
    records, notes = read_sequences(path)
    assert notes["aligned"] == (0, 3, 0)
    rows = check(assay, records)
    assert rows[0]["verdict"] == "perfect"
    assert all(o["unknown_bases"] == 0 for o in rows[0]["oligos"].values())


def test_a_dot_is_a_gap_character_too(reference):
    """Some formats pad with `.` rather than `-`."""
    sequence, rna, gaps, unknown = normalise(reference[:50] + "..." + reference[50:])
    assert sequence == reference
    assert (rna, gaps, unknown) == (0, 3, 0)


def test_a_record_that_is_nothing_but_padding_is_refused(tmp_path):
    """An empty sequence is not a sequence an assay can be checked against."""
    path = write(tmp_path, "pad", "-" * 80)
    with pytest.raises(InputError, match="empty once gap characters are removed"):
        read_sequences(path)


# --- anything else ------------------------------------------------------------------


def test_a_letter_that_is_not_a_nucleotide_code_is_unknown_not_mismatched(reference):
    """What `ncbi.py` already did with it, so the two readers now agree.

    Length is preserved: a stray letter is a base that was read badly, not padding, and
    deleting it would shift every position after it.
    """
    sequence, rna, gaps, unknown = normalise(reference[:100] + "XZ" + reference[102:])
    assert unknown == 2
    assert len(sequence) == len(reference)
    assert sequence[100:102] == "NN"


def test_an_unknown_base_under_an_oligo_is_unknown_not_drifted(tmp_path, assay, reference):
    site = reference.find("TTACAAACATTGGCCGCAAA")
    broken = reference[: site + 5] + "XXX" + reference[site + 8 :]
    path = write(tmp_path, "junk", broken)
    records, notes = read_sequences(path)
    assert notes["junk"] == (0, 0, 3)
    rows = check(assay, records)
    assert rows[0]["verdict"] == "unknown"
    assert rows[0]["oligos"]["N2-F"]["unknown_bases"] == 3
    assert rows[0]["oligos"]["N2-F"]["mismatches"] == 0


def test_the_three_rewrites_are_counted_separately(reference):
    """Each kind of character in one sequence, each in its own column."""
    mixed = (reference[:50] + "--" + reference[50:100] + "X" + reference[101:]).replace("T", "U")
    sequence, rna, gaps, unknown = normalise(mixed)
    assert gaps == 2
    assert unknown == 1
    assert rna == mixed.count("U")
    assert len(sequence) == len(reference)


# --- the other input shapes a lab will hand it ---------------------------------------


def test_line_wrapped_lowercase_and_crlf_all_read_the_same(tmp_path, assay, reference):
    """Three ways the same genome arrives, one verdict."""
    wrapped = "\n".join(reference[i : i + 60] for i in range(0, len(reference), 60))
    (tmp_path / "wrapped.fasta").write_text(f">w\n{wrapped}\n", encoding="utf-8")
    (tmp_path / "lower.fasta").write_text(f">l\n{reference.lower()}\n", encoding="utf-8")
    (tmp_path / "crlf.fasta").write_bytes(f">c\r\n{reference}\r\n".encode())
    for name in ("wrapped", "lower", "crlf"):
        assert verdicts(assay, tmp_path / f"{name}.fasta") == ["perfect"], name


def test_a_sequence_shorter_than_the_oligos_is_unknown_not_failing(tmp_path, assay):
    """No binding site by construction is not the same fact as a binding site lost."""
    assert verdicts(assay, write(tmp_path, "short", "ACGT")) == ["unknown"]
