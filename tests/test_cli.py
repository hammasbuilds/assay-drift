"""Tests for the `assay-drift check` command on user-supplied files."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.cli import InputError, main, read_fasta, read_primers  # noqa: E402
from assaydrift.match import reverse_complement  # noqa: E402

LEFT = "ACGTTGCAAGCTTGACCTGA"  # forward primer site
PROBE = "GGATCCATTGCAGTCAAGCT"
RIGHT = "TTGACGGCATAGCCTAGCAT"  # reverse primer binds the rc of this
SPACER = "A" * 30
TARGET = SPACER + LEFT + "CC" + PROBE + "CC" + RIGHT + SPACER


def write_primers(tmp_path: Path) -> Path:
    p = tmp_path / "assay.tsv"
    p.write_text(
        "# my assay\n"
        f"F\tforward\t{LEFT}\nP\tprobe\t{PROBE}\nR\treverse\t{reverse_complement(RIGHT)}\n",
        encoding="utf-8",
    )
    return p


def mutate(seq: str, index: int) -> str:
    base = "A" if seq[index] != "A" else "C"
    return seq[:index] + base + seq[index + 1 :]


def test_verdicts_cover_every_case(tmp_path, capsys):
    primers = write_primers(tmp_path)
    left_end = len(SPACER) + len(LEFT) - 1  # forward primer's 3' base
    left_start = len(SPACER)
    gap = TARGET[: len(SPACER) + len(LEFT) + 2] + "N" + TARGET[len(SPACER) + len(LEFT) + 3 :]
    fasta = tmp_path / "s.fasta"
    fasta.write_text(
        f">same\n{TARGET}\n>five\n{mutate(TARGET, left_start)}\n"
        f">three\n{mutate(TARGET, left_end)}\n>gap desc\n{gap}\n",
        encoding="utf-8",
    )
    out = tmp_path / "r.json"
    assert (
        main(["check", "--primers", str(primers), "--sequences", str(fasta), "--json", str(out)])
        == 0
    )
    rows = {r["sequence"]: r["verdict"] for r in json.loads(out.read_text())}
    assert rows == {
        "same": "perfect",
        "five": "drifted",
        "three": "likely_failing",
        "gap": "unknown",
    }
    assert "perfect 1, drifted 1, likely_failing 1, unknown 1" in capsys.readouterr().out


@pytest.mark.parametrize(
    "text, message",
    [
        ("F\tforward\n", "expected name"),
        ("F\tsideways\tACGT\n", "role must be"),
        ("F\tforward\tACGZ\n", "not IUPAC"),
        ("# nothing\n", "no oligos"),
    ],
)
def test_bad_primer_files_are_named_errors(tmp_path, text, message):
    p = tmp_path / "p.tsv"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(InputError, match=message):
        read_primers(p)


@pytest.mark.parametrize(
    "text, message",
    [
        ("ACGT\n>x\nACGT\n", "before the first"),
        ("", "no FASTA"),
        (">x\n>y\nAC\n", "empty sequence"),
    ],
)
def test_bad_fasta_files_are_named_errors(tmp_path, text, message):
    p = tmp_path / "s.fa"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(InputError, match=message):
        read_fasta(p)


def test_missing_file_exits_1(tmp_path, capsys):
    primers = write_primers(tmp_path)
    assert main(["check", "--primers", str(primers), "--sequences", str(tmp_path / "no.fa")]) == 1
    assert "no such file" in capsys.readouterr().err


def test_analyse_without_data_exits_nonzero(tmp_path):
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "analyse.py"), "--data", str(tmp_path / "x.jsonl")],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 1
    assert "run scripts/fetch.py first" in proc.stderr


def test_shipped_example_reproduces_the_readme_output(capsys):
    """The example files in the README's quickstart must keep working."""
    code = main(
        [
            "check",
            "--primers",
            str(ROOT / "examples" / "cdc_n2.tsv"),
            "--sequences",
            str(ROOT / "examples" / "two_genomes.fasta"),
        ]
    )
    out = capsys.readouterr().out
    assert code == 0
    assert "NC_045512.2             perfect" in out
    assert "synthetic_C29215T       likely_failing" in out
    assert "1mm/3'" in out
    assert "perfect 1, drifted 0, likely_failing 1, unknown 0  (of 2)" in out


class TestARecordTooShortToJudge:
    """A sequence shorter than the oligo is not evidence that the assay fails.

    `check` reported a four-base FASTA record as `likely_failing` with all three oligos
    "lost". `find` correctly returns NOT_FOUND when the target is shorter than the oligo,
    but `verdict` reached `likely_failing` through `any_oligo_lost`, and its only route to
    "unknown" was `not usable and not any_oligo_lost` - which a lost oligo closed off.

    It matters because `check` takes an arbitrary FASTA and partial or amplicon-only
    GenBank records are routine. The `drift` pipeline was never affected: the shipped
    corpus has a 29,469-base minimum.
    """

    def test_a_four_base_record_is_unknown_not_failing(self, tmp_path):
        from assaydrift.cli import check
        from assaydrift.primers import by_name

        assay = by_name("CDC N2")
        rows = check(assay, [("tiny", "ACGT")])
        assert rows[0]["verdict"] == "unknown", (
            "a record too short to hold any oligo was judged, not refused"
        )

    def test_a_record_long_enough_for_some_oligos_is_still_refused(self):
        """Partial coverage is still no basis for a verdict about the whole assay."""
        from assaydrift.cli import check
        from assaydrift.primers import by_name

        assay = by_name("CDC N2")
        shortest = min(len(o.sequence) for o in assay.oligos)
        longest = max(len(o.sequence) for o in assay.oligos)
        assert shortest != longest, "this test needs oligos of differing length"
        rows = check(assay, [("partial", "A" * ((shortest + longest) // 2))])
        assert rows[0]["verdict"] == "unknown"

    def test_a_full_length_genome_is_still_judged(self):
        """The guard must not swallow the real case."""
        from assaydrift.cli import check
        from assaydrift.primers import by_name
        from assaydrift.reference import genome

        rows = check(by_name("CDC N2"), [("ref", genome())])
        assert rows[0]["verdict"] == "perfect"
