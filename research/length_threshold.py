"""Where do the corpus's record lengths actually fall, either side of 25,000 bases?

`region_absent` turns on `MIN_WHOLE_GENOME = 25_000`: an assay whose every oligo is
missing is a lost target region if the record is a whole genome, and a fragment that
never contained the region if it is not. That is a threshold with a number in it, chosen
by judgement, in a paper whose section 4 is about exactly that kind of threshold being
uncalibrated. The threat section should say how much the number matters, and the honest
way to say it is to look.

This reads only the `length` field, so it is minutes rather than hours, and reports the
distribution near the cut together with how the classification would move at 20,000 and
30,000. If the nearest record on either side is far from the boundary, the threshold is
inert and the paper can say so with a number; if it is not, that is a finding.

    .venv/Scripts/python research/length_threshold.py --data data/sequences_big.jsonl
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "research" / "results"

# The alternatives, chosen before looking: a fifth either side of the published value.
CANDIDATES = (20_000, 25_000, 30_000)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "sequences_big.jsonl")
    parser.add_argument("--out", type=Path, default=RESULTS / "length_threshold.json")
    args = parser.parse_args()

    lengths: list[int] = []
    with args.data.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            value = record.get("length")
            if not value:
                sequence = record.get("sequence") or ""
                value = len(sequence)
            lengths.append(int(value))

    lengths.sort()
    total = len(lengths)
    below = {c: sum(1 for v in lengths if v < c) for c in CANDIDATES}

    published = 25_000
    under = [v for v in lengths if v < published]
    over = [v for v in lengths if v >= published]
    payload = {
        "what": (
            f"Record lengths in the {total:,}-record corpus, to say how much the "
            f"{published:,}-base MIN_WHOLE_GENOME threshold decides. Only the length "
            "field is read."
        ),
        "records": total,
        "published_threshold": published,
        "below_published": len(under),
        "largest_below": max(under) if under else None,
        "smallest_at_or_above": min(over) if over else None,
        "gap_around_the_threshold": (min(over) - max(under)) if under and over else None,
        "count_below_each_candidate": below,
        "classification_moves_between_candidates": {
            f"{a}->{b}": below[b] - below[a] for a, b in itertools.pairwise(CANDIDATES)
        },
        "percentiles": {
            "min": lengths[0],
            "p01": lengths[max(0, total // 100 - 1)],
            "p05": lengths[max(0, total // 20 - 1)],
            "median": lengths[total // 2],
            "max": lengths[-1],
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(json.dumps(payload, indent=1))
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
