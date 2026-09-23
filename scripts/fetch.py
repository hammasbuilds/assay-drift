"""Download a time-stratified sample of genomes from GenBank.

    python scripts/fetch.py --per-year 400

**Stratified by year on purpose.** `ncbi.search` sorts newest first, so taking
the first N of a single query returns almost entirely recent sequences - and a
drift study built on one year of data cannot show drift at all. This runs one
query per deposit year and takes an equal slice from each, so every period in
the final table has sequences behind it.

Deposit year is a proxy for collection year, and not a perfect one: a sequence
collected in December is often deposited the following spring. That is fine
here, because nothing downstream trusts it - every record is bucketed by the
`collection_date` written in the record itself, and one with no usable date is
dropped and counted. The deposit year only decides what gets *downloaded*.

Everything is cached under data/cache, so a second run costs NCBI nothing.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from assaydrift import ncbi  # noqa: E402
from assaydrift.analyze import period_of  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "sequences.jsonl"

QUERY = "SARS-CoV-2[orgn] AND complete genome[title]"
YEARS = range(2020, 2027)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--per-year", type=int, default=400, help="genomes to download per deposit year"
    )
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    written = 0
    undated = 0
    by_period: Counter = Counter()

    with args.out.open("w", encoding="utf-8") as handle:
        for year in YEARS:
            term = f'{QUERY} AND ("{year}/01/01"[PDAT] : "{year}/12/31"[PDAT])'
            available = ncbi.count(term)
            if not available:
                print(f"  {year}: nothing deposited", flush=True)
                continue
            ids = ncbi.search(term, limit=args.per_year)
            records = ncbi.fetch(ids, batch=100)

            kept = 0
            for record in records:
                period = period_of(record.collected)
                if period is None:
                    undated += 1
                    continue
                handle.write(
                    json.dumps(
                        {
                            "accession": record.accession,
                            "organism": record.organism,
                            "collected": record.collected,
                            "country": record.country,
                            "length": record.length,
                            "sequence": record.sequence,
                        }
                    )
                    + "\n"
                )
                by_period[period] += 1
                kept += 1
                written += 1
            print(
                f"  {year}: {available:>9,} available, {len(ids):>4} requested, "
                f"{len(records):>4} fetched, {kept:>4} dated",
                flush=True,
            )

    print()
    print(f"{written:,} sequences written to {args.out}")
    print(f"{undated:,} dropped for having no usable collection date")
    print(f"{time.time() - started:.0f}s")
    print()
    print("collection periods actually covered:")
    for period, count in sorted(by_period.items()):
        print(f"  {period}  {count:>5}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
