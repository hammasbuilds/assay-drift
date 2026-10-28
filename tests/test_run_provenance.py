"""The `_run` block must describe the run that produced the file, not the default one.

Three of its fields have now been wrong in the same way: each named what a plain
`python scripts/analyse.py` would have used rather than what this run used.

- `manifest` was hardcoded to `data/accessions.tsv`, so a run over 291,827 sequences
  recorded `manifest_records: 2765` beside `sequences_read: 291827` - and
  `corpus_identity_sha256`, the field whose whole purpose is to say which corpus produced
  these rates, carried the small corpus's digest.
- `reproduce` was a constant string, so the same file told a reader to rebuild it with a
  command that rebuilds a different corpus. It is the field a reader trusts most, because
  it claims to be executable.
- `tool_version` is None whenever the package is not installed, which the README makes
  the normal way to run this, so a published file recorded nothing about its own code.

These are metadata on a file that takes hours to produce, so none of them may raise -
`data/` in this checkout is a Windows junction and an earlier run over the full corpus
died at its last step formatting a path.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import analyse  # noqa: E402


def _ns(data: Path | None, granularity: str = "quarter") -> argparse.Namespace:
    return argparse.Namespace(data=data, granularity=granularity)


def test_reproduce_names_the_corpus_this_run_read() -> None:
    manifest = ROOT / "data" / "accessions_big.tsv"
    command = analyse._reproduce(_ns(ROOT / "data" / "sequences_big.jsonl"), manifest)
    assert "sequences_big.jsonl" in command
    assert "accessions_big.tsv" in command
    assert "accessions.tsv" not in command.replace("accessions_big.tsv", ""), (
        "the command still names the default manifest alongside the one that was read"
    )


def test_reproduce_omits_arguments_that_were_not_given() -> None:
    """A default run must still print a command someone can paste."""
    command = analyse._reproduce(_ns(None), ROOT / "data" / "accessions.tsv")
    assert "--data" not in command
    assert command.startswith("python scripts/fetch.py --from-manifest")
    assert "python scripts/analyse.py" in command


def test_reproduce_does_not_raise_on_a_path_outside_the_repository() -> None:
    """The reason this is a test and not an assumption.

    `Path.relative_to(ROOT)` raises for anything outside the repository, and `resolve()`
    makes that worse here because `data/` is a junction pointing elsewhere on the drive -
    so every path under the one directory this script exists to read resolves to outside
    the repository. A run over 291,827 genomes reached its last step and died on a line
    whose only job was to make a filename look tidy.
    """
    outside = Path("D:/_census_scratch/elsewhere.jsonl")
    command = analyse._reproduce(_ns(outside), ROOT / "data" / "accessions.tsv")
    assert "elsewhere.jsonl" in command


def test_reproduce_handles_a_manifest_that_is_not_there() -> None:
    command = analyse._reproduce(_ns(None), ROOT / "data" / "no_such_manifest.tsv")
    assert "--manifest" not in command


def test_the_run_records_which_code_produced_it() -> None:
    """Either an installed version or a checkout revision - not nothing.

    Outside a git checkout and without an install both are None, which is honest; what
    this forbids is the pair being silently absent from the payload.
    """
    assert hasattr(analyse, "_source_revision")
    revision = analyse._source_revision()
    assert revision is None or revision.strip(), "an empty string is not a revision"
    if revision is not None:
        assert " " not in revision
