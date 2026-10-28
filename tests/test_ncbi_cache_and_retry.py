"""The network layer, exercised without a network.

`suite-auditor coverage` reports 61 of 65 functions reached by a passing test here, and
all four it does not reach are in `ncbi.py`: `_cache_dir`, `_slot`, `_get` and `fetch`.
Being the HTTP layer is not a reason to leave them untested - every one can be driven with
a substituted opener and a temporary directory, and two of them carry docstrings recording
defects that already happened:

  * `_cache_dir` computed its path from `__file__` unconditionally, so an installed wheel
    wrote its HTTP cache inside the virtualenv - somebody else's directory, deleted with
    the environment.
  * `_get` let a truncated response escape the retry loop and end the whole run. The
    comment there says that over a few thousand requests it is rare and over the hundreds
    of thousands a large corpus needs it is certain. The run behind the published result
    made 291,827 of them.
"""

from __future__ import annotations

import gzip
import http.client
import io
import os
import sys
import time
import urllib.error
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift import ncbi  # noqa: E402


class _Response(io.BytesIO):
    """Just enough of an HTTP response for `with urlopen(...) as r: r.read()`."""

    def __enter__(self):  # noqa: ANN204
        return self

    def __exit__(self, *exc) -> None:  # noqa: ANN002
        self.close()


@pytest.fixture
def no_sleep(monkeypatch):
    """The retry loop sleeps 2, 4 then 6 seconds; a test must not."""
    monkeypatch.setattr(ncbi.time, "sleep", lambda _seconds: None)


class TestCacheDir:
    def test_the_environment_variable_wins(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path / "chosen"))
        assert ncbi._cache_dir() == tmp_path / "chosen"

    def test_a_tilde_in_the_override_is_expanded(self, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", "~/assay-cache")
        assert ncbi._cache_dir() == Path.home() / "assay-cache"

    def test_a_checkout_caches_beside_itself(self, monkeypatch) -> None:
        """Only when both markers are present, which is what tells a checkout apart."""
        monkeypatch.delenv("ASSAY_DRIFT_CACHE", raising=False)
        assert ncbi._cache_dir() == ROOT / "data" / "cache"

    def test_an_installed_package_does_not_cache_inside_its_environment(
        self, tmp_path, monkeypatch
    ) -> None:
        """The defect the docstring records.

        An installed wheel has no pyproject.toml two directories up, so the path must
        come from the user's cache directory rather than from `__file__`.
        """
        monkeypatch.delenv("ASSAY_DRIFT_CACHE", raising=False)
        fake = tmp_path / "venv" / "Lib" / "site-packages" / "assaydrift" / "ncbi.py"
        fake.parent.mkdir(parents=True)
        monkeypatch.setattr(ncbi, "__file__", str(fake))
        home = tmp_path / "home"
        monkeypatch.setenv("LOCALAPPDATA", str(home))
        monkeypatch.setenv("XDG_CACHE_HOME", str(home))
        where = ncbi._cache_dir()
        assert where == home / "assay-drift"
        assert (tmp_path / "venv") not in where.parents, (
            "an installed package is caching inside its own virtualenv"
        )


class TestSlot:
    def test_a_url_maps_to_one_stable_filename(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        once = ncbi._slot("https://example.org/a")
        assert once == ncbi._slot("https://example.org/a")
        assert once.suffix == ".gz"
        assert once.parent == tmp_path

    def test_different_urls_do_not_share_a_slot(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        assert ncbi._slot("https://example.org/a") != ncbi._slot("https://example.org/b")


class TestGet:
    def test_a_fetched_page_is_returned_and_cached(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        calls = []

        def opener(request, timeout=None):  # noqa: ANN001, ANN202
            calls.append(request.full_url)
            return _Response(b"PAYLOAD")

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", opener)
        assert ncbi._get("https://example.org/x") == "PAYLOAD"
        assert len(calls) == 1
        with gzip.open(ncbi._slot("https://example.org/x"), "rt", encoding="utf-8") as fh:
            assert fh.read() == "PAYLOAD"

    def test_a_cached_page_is_not_fetched_again(self, tmp_path, monkeypatch) -> None:
        """The reason the cache exists: a re-analysis must cost NCBI nothing."""
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        slot = ncbi._slot("https://example.org/y")
        slot.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(slot, "wt", encoding="utf-8") as fh:
            fh.write("FROM CACHE")

        def refuse(request, timeout=None):  # noqa: ANN001, ANN202
            raise AssertionError("a cached page was fetched again")

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", refuse)
        assert ncbi._get("https://example.org/y") == "FROM CACHE"

    def test_a_stale_cache_entry_is_refetched(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        slot = ncbi._slot("https://example.org/z")
        slot.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(slot, "wt", encoding="utf-8") as fh:
            fh.write("OLD")
        stale = time.time() - 10
        os.utime(slot, (stale, stale))
        monkeypatch.setattr(
            ncbi.urllib.request, "urlopen", lambda r, timeout=None: _Response(b"NEW")
        )
        assert ncbi._get("https://example.org/z", max_age=1) == "NEW"

    def test_a_truncated_response_is_retried_rather_than_ending_the_run(
        self, tmp_path, monkeypatch, no_sleep
    ) -> None:
        """The defect the comment above RETRIES records.

        `http.client.IncompleteRead` is not an OSError, so it used to escape the retry
        loop. Over the 291,827 requests behind the published result it is a certainty.
        """
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        attempts = []

        def flaky(request, timeout=None):  # noqa: ANN001, ANN202
            attempts.append(1)
            if len(attempts) == 1:
                raise http.client.IncompleteRead(b"half")
            return _Response(b"WHOLE")

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", flaky)
        assert ncbi._get("https://example.org/t") == "WHOLE"
        assert len(attempts) == 2

    def test_a_persistent_failure_raises_after_the_retries(
        self, tmp_path, monkeypatch, no_sleep
    ) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        attempts = []

        def always_fails(request, timeout=None):  # noqa: ANN001, ANN202
            attempts.append(1)
            raise urllib.error.URLError("down")

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", always_fails)
        with pytest.raises(ncbi.NCBIError, match="example.org"):
            ncbi._get("https://example.org/dead")
        assert len(attempts) == ncbi.RETRIES, (
            f"tried {len(attempts)} times, RETRIES says {ncbi.RETRIES}"
        )

    def test_nothing_is_cached_when_every_attempt_fails(
        self, tmp_path, monkeypatch, no_sleep
    ) -> None:
        """A cached failure would be served for a month."""
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))

        def always_fails(request, timeout=None):  # noqa: ANN001, ANN202
            raise urllib.error.URLError("down")

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", always_fails)
        with pytest.raises(ncbi.NCBIError):
            ncbi._get("https://example.org/dead2")
        assert not ncbi._slot("https://example.org/dead2").exists()


GENBANK = """LOCUS       AA000001
ACCESSION   AA000001
  ORGANISM  Severe acute respiratory syndrome coronavirus 2
                     /collection_date="2021-03-14"
                     /country="USA: Testland"
ORIGIN
        1 acgtacgt
//
"""


class TestFetch:
    def test_records_come_back_parsed(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        monkeypatch.setattr(
            ncbi.urllib.request,
            "urlopen",
            lambda r, timeout=None: _Response(GENBANK.encode()),
        )
        records = ncbi.fetch(["AA000001"])
        assert [r.accession for r in records] == ["AA000001"]
        assert records[0].collected == "2021-03-14"
        assert records[0].country == "USA: Testland"
        # Upper-cased by the parser, which is how every comparison downstream reads it.
        assert records[0].sequence == "ACGTACGT"

    def test_ids_are_requested_in_batches(self, tmp_path, monkeypatch) -> None:
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))
        urls = []

        def opener(request, timeout=None):  # noqa: ANN001, ANN202
            urls.append(request.full_url)
            return _Response(GENBANK.encode())

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", opener)
        ncbi.fetch([f"AA{i:06d}" for i in range(5)], batch=2)
        assert len(urls) == 3, f"5 ids at a batch of 2 is 3 requests, not {len(urls)}"

    def test_one_failed_batch_does_not_end_the_run(
        self, tmp_path, monkeypatch, no_sleep
    ) -> None:
        """The reason fetch() swallows NCBIError: a long run must not die on one batch."""
        monkeypatch.setenv("ASSAY_DRIFT_CACHE", str(tmp_path))

        def opener(request, timeout=None):  # noqa: ANN001, ANN202
            if "AA000000" in request.full_url:
                raise urllib.error.URLError("down")
            return _Response(GENBANK.encode())

        monkeypatch.setattr(ncbi.urllib.request, "urlopen", opener)
        records = ncbi.fetch(["AA000000", "AA000001"], batch=1)
        assert [r.accession for r in records] == ["AA000001"], (
            "the surviving batch must still be returned"
        )
