"""Check your own primers against your own sequences.

    assay-drift check --primers my_assay.tsv --sequences genomes.fasta [--json out.json]
    python -m assaydrift check ...

Primers file: tab-separated, one oligo per line, `name<TAB>role<TAB>sequence`,
role being forward, reverse or probe. Blank lines and `#` comments are skipped.
All oligos in the file are treated as one assay.

Prints, per sequence, each oligo's mismatches, unknown bases and 3' primer
mismatches, and a verdict using the same rules as the GenBank study:
`perfect`, `drifted` (mismatches but still expected to amplify),
`likely_failing`, or `unknown` (an `N` under an oligo - not evidence either way).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .analyze import AssayResult
from .match import find_oligo
from .primers import Assay, Oligo

ROLES = ("forward", "reverse", "probe")


class InputError(ValueError):
    """A primers or sequence file that cannot be read as described."""


def read_primers(path: Path) -> Assay:
    oligos: list[Oligo] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) != 3:
            raise InputError(f"{path}:{number}: expected name<TAB>role<TAB>sequence")
        name, role, sequence = parts
        if role.lower() not in ROLES:
            raise InputError(f"{path}:{number}: role must be one of {ROLES}, got {role!r}")
        if not sequence:
            raise InputError(f"{path}:{number}: empty sequence")
        try:
            oligos.append(Oligo(name, sequence.upper(), role.lower()))
        except ValueError as exc:
            raise InputError(f"{path}:{number}: {exc}") from exc
    if not oligos:
        raise InputError(f"{path}: no oligos found")
    return Assay(name=path.stem, target_gene="", source=str(path), oligos=oligos)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].split()[0] if line[1:].strip() else f"seq{len(records) + 1}"
            chunks = []
        elif name is None:
            raise InputError(f"{path}:{number}: sequence data before the first '>' header")
        else:
            chunks.append(line.upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise InputError(f"{path}: no FASTA records found")
    empty = [n for n, s in records if not s]
    if empty:
        raise InputError(f"{path}: empty sequence for {empty[0]}")
    return records


def too_short(sequence: str, assay: Assay) -> int:
    """How many oligos cannot physically fit in this sequence.

    A record shorter than an oligo has no binding site for it by construction, which is
    not the same fact as a binding site having been lost. Partial and amplicon-only
    GenBank records are routine input to `check`.
    """
    return sum(1 for oligo in assay.oligos if len(sequence) < len(oligo.sequence))


def verdict(result: AssayResult, sequence: str = "", assay: Assay | None = None) -> str:
    # FIRST, because "there was nowhere for the oligo to bind" is not evidence about the
    # assay. `find` returns NOT_FOUND when the target is shorter than the oligo, and
    # `any_oligo_lost` then carried that into `likely_failing` - so a four-base record was
    # reported as an assay that would no longer amplify, with three oligos "lost".
    if assay is not None and sequence and too_short(sequence, assay):
        return "unknown"
    if not result.usable and not result.any_oligo_lost:
        return "unknown"
    if result.likely_failing():
        return "likely_failing"
    return "perfect" if result.perfect() else "drifted"


def check(assay: Assay, sequences: list[tuple[str, str]]) -> list[dict]:
    rows = []
    for name, sequence in sequences:
        result = AssayResult(accession=name, assay=assay.name, period="")
        for oligo in assay.oligos:
            result.hits[oligo.name] = find_oligo(sequence, oligo.sequence, oligo.role)
        rows.append(
            {
                "sequence": name,
                "verdict": verdict(result, sequence, assay),
                "oligos": {
                    o: {
                        "found": h.found,
                        "role": h.role,
                        "position": h.position,
                        "mismatches": h.mismatches if h.found else None,
                        "unknown_bases": h.ambiguous,
                        # `usable` is a property, so it was invisible to anyone reading
                        # this file - and its own docstring is "whether this measurement
                        # can go into a rate", because an alignment with unknown bases
                        # under it is not evidence either way. Leaving it out meant a
                        # consumer computing a rate from this JSON would count exactly the
                        # rows the property exists to exclude.
                        "usable": h.usable if h.found else False,
                        # The count, not only the boolean below. On a probe
                        # `blocks_extension` is 0 even with 3' mismatches, because a probe
                        # is never extended - so the boolean alone cannot distinguish "no
                        # 3' mismatches" from "3' mismatches that do not block".
                        "three_prime_mismatches": h.three_prime_mismatches if h.found else None,
                        "three_prime_blocking": h.blocks_extension if h.found else None,
                        # Non-zero means the site was explained by a deletion or insertion
                        # rather than by a run of mismatches, which changes what the
                        # mismatch count above is a count of.
                        "indels": h.indels if h.found else None,
                    }
                    for o, h in result.hits.items()
                },
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="assay-drift",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = ap.add_subparsers(dest="command", required=True)
    c = sub.add_parser("check", help="score your primers against your FASTA sequences")
    c.add_argument("--primers", type=Path, required=True)
    c.add_argument("--sequences", type=Path, required=True)
    c.add_argument("--json", type=Path, help="also write the per-sequence results here")
    args = ap.parse_args(argv)

    for path in (args.primers, args.sequences):
        if not path.is_file():
            print(f"assay-drift: no such file: {path}", file=sys.stderr)
            return 1
    try:
        assay = read_primers(args.primers)
        sequences = read_fasta(args.sequences)
    except (InputError, UnicodeDecodeError) as exc:
        print(f"assay-drift: {exc}", file=sys.stderr)
        return 1

    rows = check(assay, sequences)
    names = [o.name for o in assay.oligos]
    print(f"{'sequence':<24}{'verdict':<16}" + "".join(f"{n:>12}" for n in names))
    for row in rows:
        cells = []
        for n in names:
            o = row["oligos"][n]
            if not o["found"]:
                cells.append("lost")
            elif o["unknown_bases"]:
                cells.append(f"{o['unknown_bases']}N")
            else:
                cells.append(f"{o['mismatches']}mm" + ("/3'" if o["three_prime_blocking"] else ""))
        print(
            f"{row['sequence'][:23]:<24}{row['verdict']:<16}" + "".join(f"{c:>12}" for c in cells)
        )
    counts = {
        v: sum(r["verdict"] == v for r in rows)
        for v in ("perfect", "drifted", "likely_failing", "unknown")
    }
    print("\n" + ", ".join(f"{k} {v}" for k, v in counts.items()) + f"  (of {len(rows)})")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=1), encoding="utf-8")
        print(f"written to {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
