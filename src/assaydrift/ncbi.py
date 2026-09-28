"""Sequences from NCBI, with the dates that make drift measurable.

Everything here goes through E-utilities, which is the documented way to
query GenBank programmatically and is meant to be polled. The rules NCBI
publishes are followed rather than worked around: an identifying User-Agent,
three requests a second without a key, and `usehistory` so a large result set
is paged on the server instead of being re-searched.

**The collection date is the point.** A drift study that cannot say *when* a
sequence was collected can only report a rate, and a rate without a trend
cannot distinguish "this assay was always imperfect" from "this assay is going
blind". So a record with no usable date is dropped and counted, never
defaulted to today.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _cache_dir() -> Path:
    """Where downloaded GenBank pages are kept.

    `ASSAY_DRIFT_CACHE` wins; otherwise `data/cache` beside the checkout when
    the package is running from one, and the user's cache directory when it is
    not. Computing it from `__file__` unconditionally meant an installed wheel
    wrote its HTTP cache inside the virtualenv, which is somebody else's
    directory and gets deleted with the environment.
    """
    override = os.environ.get("ASSAY_DRIFT_CACHE")
    if override:
        return Path(override).expanduser()
    root = Path(__file__).resolve().parents[2]
    if (root / "pyproject.toml").exists() and (root / "src").is_dir():
        return root / "data" / "cache"
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local"
    else:
        base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "assay-drift"


# Resolved on each use rather than pinned at import, so setting the
# environment variable in a test or a shell actually takes effect.
CACHE = _cache_dir()

# NCBI asks for a real address so they can contact you before blocking you.
# Override it when you run this yourself, so a rate-limit warning reaches you
# and not the person who happened to write the tool.
UA = os.environ.get(
    "ASSAY_DRIFT_USER_AGENT",
    "assay-drift/0.1 (research tool; https://github.com/hammasbuilds/assay-drift)",
)

# Three per second is the documented limit without an API key. The sleep is
# deliberate and not configurable downward: this is somebody else's free
# public service and the whole project depends on continuing to be allowed in.
MIN_INTERVAL = 0.34
_last_call = 0.0

TIMEOUT = 60.0
RETRIES = 3


class NCBIError(RuntimeError):
    """E-utilities could not be read."""


@dataclass
class Record:
    """One GenBank sequence, reduced to what a drift study needs."""

    accession: str
    organism: str
    sequence: str
    collected: str = ""  # ISO-ish; "" when the record does not say
    country: str = ""
    length: int = 0

    @property
    def year(self) -> int | None:
        found = re.search(r"(19|20)\d{2}", self.collected)
        return int(found.group(0)) if found else None


class DataError(ValueError):
    """A local sequence dump could not be read."""


REQUIRED_FIELDS = ("accession", "organism", "sequence", "collected")


def read_jsonl(path: Path, *, on_bad_line=None):
    """Records from a JSONL dump, one per line, reporting where it went wrong.

    `scripts/fetch.py` writes this format and users point the tools at their own
    files, so a malformed line is an ordinary event rather than a bug. It is
    reported as `path:line: what` and, when `on_bad_line` is given, handed over
    and skipped instead of ending the run - a single truncated line should not
    cost 2,000 good ones.

    `length` and `country` are optional; the four fields above are not, because
    a record without a sequence or a collection date cannot enter the study at
    all and silently defaulting either is how a drift study invents a trend.
    """
    if not path.exists():
        raise DataError(f"no data at {path} - run scripts/fetch.py first")
    empty = True
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
                if not isinstance(row, dict):
                    raise DataError(f"expected a JSON object, got {type(row).__name__}")
                missing = [field for field in REQUIRED_FIELDS if field not in row]
                if missing:
                    raise DataError(f"missing field(s) {', '.join(missing)}")
                record = Record(
                    accession=row["accession"],
                    organism=row["organism"],
                    sequence=row["sequence"],
                    collected=row["collected"],
                    country=row.get("country", ""),
                    length=int(row.get("length") or len(row["sequence"])),
                )
            except (json.JSONDecodeError, DataError, TypeError, ValueError) as exc:
                problem = DataError(f"{path}:{number}: {exc}")
                if on_bad_line is None:
                    raise problem from exc
                on_bad_line(problem)
                continue
            empty = False
            yield record
    if empty:
        raise DataError(f"{path} contains no usable records")


def _slot(url: str) -> Path:
    return _cache_dir() / f"{hashlib.blake2b(url.encode(), digest_size=10).hexdigest()}.gz"


def _get(url: str, max_age: float = 30 * 86400) -> str:
    """Fetch, cached for a month.

    GenBank records do not change once deposited, so a long cache is correct
    rather than merely convenient — and it means a re-analysis costs NCBI
    nothing at all.
    """
    global _last_call

    slot = _slot(url)
    if slot.exists() and time.time() - slot.stat().st_mtime < max_age:
        with gzip.open(slot, "rt", encoding="utf-8") as handle:
            return handle.read()

    last: Exception | None = None
    for attempt in range(RETRIES):
        wait = MIN_INTERVAL - (time.time() - _last_call)
        if wait > 0:
            time.sleep(wait)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                text = response.read().decode("utf-8", "replace")
            _last_call = time.time()
            slot.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(slot, "wt", encoding="utf-8") as handle:
                handle.write(text)
            return text
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            _last_call = time.time()
            time.sleep(2.0 * (attempt + 1))
    raise NCBIError(f"{url}: {last}")


def count(term: str) -> int:
    """How many records match, without downloading any of them."""
    url = f"{BASE}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(term)}&retmax=0"
    found = re.search(r"<Count>(\d+)</Count>", _get(url))
    return int(found.group(1)) if found else 0


def search(term: str, limit: int = 500) -> list[str]:
    """Accession ids for a query, newest first.

    Sorted by date so a capped run samples *recent* sequences. Taking the
    first N in default order would sample whatever GenBank happens to return,
    which for a drift study means mostly old records — exactly the ones that
    cannot show recent drift.
    """
    ids: list[str] = []
    step = min(limit, 500)
    while len(ids) < limit:
        url = (
            f"{BASE}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(term)}"
            f"&retmax={step}&retstart={len(ids)}&sort=date"
        )
        page = re.findall(r"<Id>(\d+)</Id>", _get(url))
        if not page:
            break
        ids += page
        if len(page) < step:
            break
    return ids[:limit]


_ORGANISM = re.compile(r"^\s{2}ORGANISM\s+(.+)$", re.M)
_COLLECTED = re.compile(r'/collection_date="([^"]+)"')
_COUNTRY = re.compile(r'/(?:country|geo_loc_name)="([^"]+)"')
_ACCESSION = re.compile(r"^ACCESSION\s+(\S+)", re.M)
_ORIGIN = re.compile(r"\nORIGIN\s*\n(.*?)\n//", re.S)


def _one_record(block: str) -> Record | None:
    accession = _ACCESSION.search(block)
    origin = _ORIGIN.search(block)
    if not accession or not origin:
        return None
    sequence = re.sub(r"[^acgtnACGTN]", "", origin.group(1)).upper()
    if not sequence:
        return None
    organism = _ORGANISM.search(block)
    collected = _COLLECTED.search(block)
    country = _COUNTRY.search(block)
    return Record(
        accession=accession.group(1),
        organism=organism.group(1).strip() if organism else "",
        sequence=sequence,
        collected=collected.group(1) if collected else "",
        country=country.group(1) if country else "",
        length=len(sequence),
    )


def fetch(ids: list[str], batch: int = 100) -> list[Record]:
    """Full GenBank records for these ids, in batches."""
    out: list[Record] = []
    for start in range(0, len(ids), batch):
        chunk = ids[start : start + batch]
        url = f"{BASE}/efetch.fcgi?db=nuccore&id={','.join(chunk)}&rettype=gb&retmode=text"
        try:
            text = _get(url)
        except NCBIError:
            continue  # one bad batch must not end a 500-record run
        for block in text.split("\nLOCUS       "):
            record = _one_record(block if block.startswith("LOCUS") else "LOCUS       " + block)
            if record:
                out.append(record)
    return out
