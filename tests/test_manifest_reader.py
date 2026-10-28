"""The manifest is the repository's claim to reproducibility, so it must be committed.

`scripts/fetch.py --from-manifest` exists because the year-window path takes the newest
deposits in each year and the recent windows are still filling, so it cannot reproduce an
earlier run. The manifest can: it is the identity of every record behind the committed
results.

That claim did not hold for the large corpus. `.gitignore` keeps `data/*` except
`data/accessions.tsv`, and `results/drift_big.json` - published from 291,827 genomes -
was built from `data/accessions_big.tsv`, which is 13.4 MB and was not in git. The
results file carried a `reproduce` command naming a manifest nobody who cloned the
repository would have.

It is now committed gzipped, 1.54 MB, and both scripts read either form.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import analyse  # noqa: E402
import fetch  # noqa: E402

READERS = pytest.mark.parametrize(
    "read", [fetch._manifest_lines, analyse._manifest_lines], ids=["fetch", "analyse"]
)

MANIFEST = (
    "# header\n"
    "# records: 2\n"
    "# identity sha256: abc123\n"
    "accession\tcollected\tcountry\tlength\n"
    "AA000001\t2021-01-02\tUSA: Testland\t29800\n"
    "AA000002\t2021-01-03\tJapan:Shizuoka, west\t29755\n"
)


@READERS
def test_a_plain_manifest_reads(tmp_path: Path, read) -> None:
    path = tmp_path / "m.tsv"
    path.write_text(MANIFEST, encoding="utf-8")
    assert read(path) == MANIFEST.splitlines()


@READERS
def test_a_gzipped_manifest_reads_identically(tmp_path: Path, read) -> None:
    plain = tmp_path / "m.tsv"
    plain.write_text(MANIFEST, encoding="utf-8")
    packed = tmp_path / "m.tsv.gz"
    packed.write_bytes(gzip.compress(MANIFEST.encode("utf-8"), mtime=0))
    assert read(packed) == read(plain)


@READERS
def test_a_wrapped_country_survives_the_round_trip(tmp_path: Path, read) -> None:
    """The defect the manifest format already suffered once.

    A GenBank country qualifier wraps across lines, and kept raw it split 14,247 of
    291,827 rows across two lines in this very file. A reader that is line-oriented must
    see one line per record.
    """
    rows = MANIFEST + "AA000003\t2021-01-04\tPhilippines: NCR, Makati City\t29700\n"
    path = tmp_path / "m.tsv.gz"
    path.write_bytes(gzip.compress(rows.encode("utf-8"), mtime=0))
    lines = [ln for ln in read(path) if ln and not ln.startswith(("#", "accession\t"))]
    assert len(lines) == 3
    assert all(len(ln.split("\t")) == 4 for ln in lines)


def test_the_large_manifest_is_in_the_repository() -> None:
    """The point of the whole change: results/drift_big.json must be rebuildable.

    Checked against git rather than the filesystem, because the uncompressed file is
    present here and ignored - which is exactly the state that made the published result
    irreproducible.
    """
    import subprocess

    big = ROOT / "results" / "drift_big.json"
    if not big.is_file():
        pytest.skip("the large result is not in this checkout")
    out = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "--error-unmatch",
         "data/accessions_big.tsv.gz"],
        capture_output=True, text=True,
    )
    assert out.returncode == 0, (
        "results/drift_big.json is published but data/accessions_big.tsv.gz is not "
        "tracked, so nobody who clones this repository can rebuild its corpus"
    )


def test_the_committed_manifest_counts_what_it_says(tmp_path: Path) -> None:
    """Its own header must match its rows - the check that caught a split-row manifest."""
    path = ROOT / "data" / "accessions_big.tsv.gz"
    if not path.is_file():
        pytest.skip("the large manifest is not in this checkout")
    lines = fetch._manifest_lines(path)
    stated = next(
        int(ln.split(":", 1)[1]) for ln in lines if ln.startswith("# records:")
    )
    rows = [ln for ln in lines if ln and not ln.startswith(("#", "accession\t"))]
    assert len(rows) == stated
    assert all(len(ln.split("\t")) == 4 for ln in rows), "a row is split across lines"
