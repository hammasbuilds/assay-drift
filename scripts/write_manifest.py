"""Write the accession manifest for the corpus that produced `results/drift.json`.

    python scripts/write_manifest.py

`data/sequences.jsonl` is 79 MB of somebody else's database and is rightly gitignored, so
the .gitignore says `scripts/fetch.py` rebuilds it. It rebuilds *a* corpus, not *the*
corpus: `fetch.py` asks NCBI for one deposit-year window at a time and takes the newest N,
and for deposit years 2025 and 2026 that window is still filling - 70,841 and 11,361
records were available when it last ran. Re-running today returns a different 767
sequences for those two years, which is 28% of the corpus and the whole basis of the CDC
N2 "2026" figures.

The manifest is the missing 30 KB: the accession and collection date of every record that
went into the committed results, so anyone can fetch exactly those and get exactly those
numbers. `fetch.py --from-manifest` reads it back.

It also records the submitter-cluster structure per period, because the published rates
rest on far fewer independent observations than their sequence counts suggest - see
`effective_n` in `results/drift.json`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift.analyze import period_of  # noqa: E402

CORPUS = ROOT / "data" / "sequences.jsonl"
MANIFEST = ROOT / "data" / "accessions.tsv"


def submitter_block(accession: str, country: str) -> str:
    """A proxy for "same submitting laboratory, same batch".

    GenBank hands out consecutive accessions to one submission, so the accession prefix
    plus the thousand-block, together with the reported country, groups records that
    arrived together. It is a proxy and nothing more - the records carry no submitter
    field - but it is enough to show that a quarter with 159 sequences can hold four
    groups, and that one of them can be 97% of it.
    """
    match = re.match(r"([A-Za-z]+)(\d+)", accession or "")
    stem = f"{match.group(1)}{int(match.group(2)) // 1000}" if match else (accession or "?")
    return f"{country or '?'}|{stem}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--corpus", type=Path, default=CORPUS)
    ap.add_argument("--out", type=Path, default=MANIFEST)
    args = ap.parse_args()

    if not args.corpus.is_file():
        print(f"missing {args.corpus}: run scripts/fetch.py first", file=sys.stderr)
        return 1

    rows: list[tuple[str, str, str, str]] = []
    digest = hashlib.sha256()
    with args.corpus.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            accession = str(record.get("accession") or "")
            collected = str(record.get("collected") or "")
            country = str(record.get("country") or "")
            rows.append((accession, collected, country, str(record.get("length") or "")))
            # Over the identity of the records, not the sequence bytes, so the hash is
            # stable against a GenBank record being revised in ways that do not change
            # which records these are.
            digest.update(f"{accession}\t{collected}\n".encode())

    rows.sort()
    header = (
        "# Every record in the corpus behind results/drift.json, so it can be rebuilt\n"
        "# exactly. scripts/fetch.py --from-manifest reads this file.\n"
        f"# records: {len(rows)}\n"
        f"# identity sha256: {digest.hexdigest()}\n"
        "accession\tcollected\tcountry\tlength\n"
    )
    args.out.write_text(header + "\n".join("\t".join(r) for r in rows) + "\n", encoding="utf-8")

    periods: dict[str, set[str]] = {}
    counts: dict[str, int] = {}
    for accession, collected, country, _ in rows:
        period = period_of(collected, "quarter")
        if not period:
            continue
        periods.setdefault(period, set()).add(submitter_block(accession, country))
        counts[period] = counts.get(period, 0) + 1

    print(f"wrote {args.out.relative_to(ROOT)}: {len(rows):,} records")
    print(f"identity sha256: {digest.hexdigest()}")
    print(f"\n{'period':<10}{'records':>9}{'clusters':>10}")
    for period in sorted(counts):
        print(f"{period:<10}{counts[period]:>9}{len(periods[period]):>10}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
