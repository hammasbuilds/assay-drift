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
