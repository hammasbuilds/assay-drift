"""How often is a real mismatch unscorable, because no published parameters cover it?

    .venv-research/Scripts/python research/tandem_gap.py [--data data/sequences_big.jsonl]

Section 4.6 scores a mismatch thermodynamically with the nearest-neighbour model, and names
a limitation it does not size: **there are no published parameters for a tandem mismatch** -
two mismatches at adjacent positions of one duplex. Allawi & SantaLucia's internal-mismatch
tables and SantaLucia & Peyret's terminal ones each describe a mismatch whose neighbours are
Watson-Crick pairs. Put two mismatches side by side and the stack between them appears in no
table, so `thermo.duplex` returns `None` and the alignment has no cost.

The limitation was discovered rather than anticipated: the branch that handles it was
written with a comment calling it unreachable, and the corpus contains the case
(C28311T/C28312T under the CDC N1 probe). A limitation that has been named but not measured
is still a hole in the section - "some alignments cannot be scored" is compatible with one
alignment and with half of them, and those are different papers. This counts them.

Three denominators, because the honest share depends on which question is being asked:

  * over **every** genome-oligo alignment, which is small by construction, since most
    alignments are perfect;
  * over the alignments that carry **any** mismatch, which is the population section 4.6's
    model is applied to and the figure the section needs;
  * over the alignments the section would actually **score**, which excludes the ones the
    published analysis already drops for unknown bases - an alignment with an N under the
    oligo is not evidence either way and was never going to be scored.

Alignments that needed an indel are counted and excluded: those are not a window against a
window, so "adjacent" is not defined on them without choosing a gap placement, and the
nearest-neighbour model has nothing to say about them either way.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

RESULTS = ROOT / "research" / "results"


def mismatch_positions(window: str, oligo: str) -> list[int] | None:
    """Indices where the target *definitely* disagrees with the oligo, or None if unscorable.

    The comparison is `match.compare_base`, imported rather than written again. A plain
    `target != wanted` is wrong here and wrong in the direction that inflates the result:
    the published oligos carry deliberate degenerate bases, and a target `A` under a primer
    `R` - which means A or G - is a match, not a mismatch. Counting it as one would have
    manufactured adjacent mismatch pairs at exactly the positions the primer designers
    chose to make tolerant, and this file exists to count a *refusal*, so inventing
    mismatches inflates the very number it reports. The repository has had to correct an
    overstated-drift bug before; `compare_base` is where that correction lives.

    `None` for a window that is short or where any position is AMBIGUOUS - an unknown base
    under the oligo is not evidence either way, which is the predicate the published rates
    already use. Keeping those apart from the tandem count matters because they are two
    different reasons an alignment cannot be scored.
    """
    from assaydrift.match import AMBIGUOUS, MISMATCH, compare_base  # noqa: PLC0415

    if len(window) != len(oligo):
        return None
    out = []
    for index, (target, wanted) in enumerate(zip(window, oligo, strict=True)):
        verdict = compare_base(wanted, target)
        if verdict == AMBIGUOUS:
            return None
        if verdict == MISMATCH:
            out.append(index)
    return out


def has_tandem(positions: list[int]) -> bool:
    """Two mismatches at adjacent indices: the case no published table covers."""
    return any(b - a == 1 for a, b in zip(positions, positions[1:], strict=False))


def scan(data: Path, limit: int = 0, checkpoint: Path | None = None, every: int = 20_000) -> dict:
    from analyse import load  # noqa: PLC0415

    from assaydrift.match import find_oligo, reverse_complement  # noqa: PLC0415
    from assaydrift.primers import ASSAYS  # noqa: PLC0415

    per_oligo: dict[str, dict] = {}
    for assay in ASSAYS:
        for oligo in assay.oligos:
            per_oligo[f"{assay.name}|{oligo.name}"] = {
                "assay": assay.name,
                "oligo": oligo.name,
                "role": oligo.role,
                "length": len(oligo.sequence),
                "alignments": 0,
                "with_a_mismatch": 0,
                "scorable": 0,
                "tandem": 0,
                "indel": 0,
                "unknown_or_short": 0,
                "worst_run": 0,
            }
    seen = 0
    for record in load(data):
        sequence = record.sequence
        if not sequence:
            continue
        seen += 1
        for assay in ASSAYS:
            for oligo in assay.oligos:
                row = per_oligo[f"{assay.name}|{oligo.name}"]
                wanted = (
                    reverse_complement(oligo.sequence.upper())
                    if oligo.role == "reverse"
                    else oligo.sequence.upper()
                )
                hit = find_oligo(sequence, oligo.sequence, oligo.role)
                if not hit.found:
                    continue
                row["alignments"] += 1
                if hit.indels:
                    row["indel"] += 1
                    continue
                window = sequence[hit.position : hit.position + len(wanted)].upper()
                positions = mismatch_positions(window, wanted)
                if positions is None:
                    row["unknown_or_short"] += 1
                    continue
                row["scorable"] += 1
                if not positions:
                    continue
                row["with_a_mismatch"] += 1
                # The longest run of adjacent mismatches, because the gap is not only
                # about pairs: three in a row has no parameters either, and a model that
                # refused pairs while silently scoring triples would be worse than one
                # that refuses both.
                run = best = 1
                for a, b in zip(positions, positions[1:], strict=False):
                    run = run + 1 if b - a == 1 else 1
                    best = max(best, run)
                row["worst_run"] = max(row["worst_run"], best)
                if best >= 2:
                    row["tandem"] += 1
        if limit and seen >= limit:
            break
        if seen % 5000 == 0:
            print(f"  {seen:,} genomes", flush=True)
        # A partial result on disk every so often. One pass over 291,827 genomes is hours
        # of CPU, and this machine is shared: a run that dies at the fifth hour with
        # nothing written has measured nothing. The checkpoint is explicitly marked
        # partial, because a share computed over 60,000 genomes is a real number about a
        # different population and must not be mistaken for the corpus.
        if checkpoint is not None and seen % every == 0:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text(
                json.dumps(_payload(data.name, seen, per_oligo, partial=True), indent=1),
                encoding="utf-8",
            )
    print(f"  {seen:,} genomes scanned", flush=True)

    return _payload(data.name, seen, per_oligo, partial=False)


def _payload(corpus: str, seen: int, per_oligo: dict, *, partial: bool) -> dict:
    def share(numerator: int, denominator: int) -> float | None:
        return round(numerator / denominator, 6) if denominator else None

    totals = {
        key: sum(row[key] for row in per_oligo.values())
        for key in (
            "alignments",
            "with_a_mismatch",
            "scorable",
            "tandem",
            "indel",
            "unknown_or_short",
        )
    }
    return {
        "corpus": corpus,
        # True while the scan is still running. Every share below is then over the genomes
        # read so far, which is a population nobody chose.
        "partial": partial,
        "genomes": seen,
        "oligos": len(per_oligo),
        "totals": totals,
        "worst_run_anywhere": max((row["worst_run"] for row in per_oligo.values()), default=0),
        "shares": {
            "of_all_alignments": share(totals["tandem"], totals["alignments"]),
            # The figure section 4.6 needs: among the alignments the model is asked to
            # score, how many does it refuse?
            "of_mismatched_alignments": share(totals["tandem"], totals["with_a_mismatch"]),
            "of_scorable_alignments": share(totals["tandem"], totals["scorable"]),
        },
        "denominators": {
            "of_all_alignments": "every alignment found, perfect ones included",
            "of_mismatched_alignments": (
                "alignments carrying at least one definite mismatch - the population "
                "section 4.6's model is applied to"
            ),
            "of_scorable_alignments": (
                "alignments with no unknown base and no indel, mismatched or not"
            ),
        },
        "per_oligo": {
            key: {
                **row,
                "tandem_share_of_mismatched": share(row["tandem"], row["with_a_mismatch"]),
            }
            for key, row in sorted(
                per_oligo.items(),
                key=lambda kv: -(kv[1]["tandem"] / max(1, kv[1]["with_a_mismatch"])),
            )
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    parser.add_argument("--limit", type=int, default=0, help="stop after this many genomes")
    parser.add_argument("--out", type=Path, default=RESULTS / "tandem_gap.json")
    args = parser.parse_args()
    if not args.data.exists():
        raise SystemExit(f"no sequence store at {args.data}")
    print(f"tandem-mismatch scan over {args.data.name}")
    payload = scan(args.data, args.limit, args.out.with_suffix(".partial.json"))
    if payload["partial"]:
        raise SystemExit("the scan returned a partial payload, which should be impossible")

    totals = payload["totals"]
    print()
    print(f"  alignments found           {totals['alignments']:>12,}")
    print(f"    needed an indel          {totals['indel']:>12,}")
    print(f"    unknown base or short    {totals['unknown_or_short']:>12,}")
    print(f"    scorable                 {totals['scorable']:>12,}")
    print(f"      carrying a mismatch    {totals['with_a_mismatch']:>12,}")
    print(f"      unscorable (tandem)    {totals['tandem']:>12,}")
    print()
    for name, value in payload["shares"].items():
        print(f"  {name:<28}{'-' if value is None else f'{value:.6f}'}")
    print(f"  longest run of adjacent mismatches anywhere: {payload['worst_run_anywhere']}")
    print()
    print("  per oligo, worst first:")
    print(f"    {'oligo':24}{'mismatched':>12}{'tandem':>9}{'share':>10}{'run':>6}")
    for key, row in payload["per_oligo"].items():
        share = row["tandem_share_of_mismatched"]
        # Padded in the same width whether or not there is a share. Without this an
        # oligo with no mismatched alignment printed "0-" where the column should be,
        # which reads as a value rather than as an absent one.
        cell = "-" if share is None else f"{share:.4f}"
        print(
            f"    {key:24}{row['with_a_mismatch']:>12,}{row['tandem']:>9,}"
            f"{cell:>10}{row['worst_run']:>6}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwritten to {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
