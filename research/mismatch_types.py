"""Which *kinds* of mismatch does the corpus actually produce near a 3' end?

    .venv-research/Scripts/python research/mismatch_types.py [--data data/sequences_big.jsonl]

Section 4 asks where a mismatch lands and how many there are. It does not ask **what the
mismatch is**, and the wet-lab literature says that is the larger effect of the three.
Stadhouders et al. (2010) measured primer-template mismatches in the 3'-end region with the
5'-nuclease assay and report that single mismatches range

    "from minor (<1.5 cycle threshold, eg, A-C, C-A, T-G, G-T) to severe impact
     (>7.0 cycle threshold, eg, A-A, G-A, A-G, C-C) on PCR amplification"

- a fivefold spread across types, for one mismatch, against the roughly one cycle that
Otwell et al. measure for a single 3'-proximal mismatch irrespective of type. A rule that
counts mismatches and measures their distance from the 3' end is blind to the dimension
with the widest range.

That matters for this paper's own withdrawal, in a direction that could go either way. §4.3
withdraws Charite RdRp's rate on the grounds that its flagged genomes carry a *single*
3'-proximal mismatch, which the bench says costs about a cycle. If the mismatches the
corpus actually produces there are Stadhouders' severe types, the withdrawal is weaker than
stated and should say so. If they are the benign ones, it is better supported than the
argument that was made for it. Either way the question is answerable from data already in
the repository, and leaving it unasked is the kind of omission §5 exists to list.

Two independent classifications are reported for each observed mismatch, which is the point
of the exercise:

  * **Stadhouders' own grouping**, from the sentence quoted above. It names four minor and
    four severe ordered pairs, which collapse to the unordered types A-C, G-T (minor) and
    A-A, A-G, C-C (severe). Three unordered types are named in neither group - C-T, G-G,
    T-T - and are reported as unclassified rather than assigned to the nearer one.
  * **This repository's own nearest-neighbour model**, as a destabilisation in kcal/mol
    from `assaydrift.thermo`. It is an independent route to the same question: published
    melting thermodynamics against published Ct measurements, neither derived from the
    other. Whether they agree about which types are cheap is worth knowing before either is
    leaned on.

**On how long this takes, and why a progress estimate from the start of the corpus is
wrong.** The first 100,000 genomes scan at about 0.019 s each; genomes around index
100,000 take about 0.14 s, seven times more. The corpus is ordered, later deposits are
further from the 2019 reference, and `match.find` returns immediately only on a *perfect*
site - so a genome with a mutation under an oligo falls through to the seed scan over the
whole sequence. Cost therefore rises with exactly the quantity this script measures, and
an ETA taken from the first tenth of the corpus understates the rest by most of a day.

The same fact rules out a shortcut. A prefix of the corpus is not a sample of it: over the
2,765-genome corpus G-T is 78% of 3'-proximal primer mismatches and over the first 100,000
of the large one it is 96%, so stopping early would report a composition that belongs to
the earliest deposits rather than to the corpus. Either the whole thing is scanned or a
sample is drawn across all of it; a checkpoint is for resuming, not for publishing.

The mismatch is reported as the unordered pair {primer base, template base}, because
Stadhouders' own lists are symmetric - A-C and C-A are both minor, A-G and G-A both severe -
so the pair carries the classification without a convention for which side comes first.
Ambiguous target bases are excluded through `match.compare_base`, as everywhere else here.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

RESULTS = ROOT / "research" / "results"

# Stadhouders et al. 2010, J Mol Diagn 12(1):109-117, from the abstract. Ordered pairs as
# printed, collapsed here to unordered types; both orders appear in their own lists, which
# is why the collapse loses nothing.
# Outer flank for the thermodynamic contexts below: long enough that the mismatch is an
# internal one rather than a near-terminal one, and mixed so it carries no strong bias of
# its own. Only the immediately neighbouring bases enter a nearest-neighbour calculation,
# so its composition does not affect the mismatch term - its length does.
FLANK = "GCATGC"

MINOR = frozenset({frozenset("AC"), frozenset("GT")})
SEVERE = frozenset({frozenset("AA"), frozenset("AG"), frozenset("CC")})


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def duplex_pair(oligo_base: str, genome_base: str) -> frozenset[str] | None:
    """The bases that actually face each other, which is not the pair you first write down.

    The aligner compares an oligo against the strand the oligo's own sequence appears in,
    so a disagreement is recorded as (oligo base, genome base). The oligo does not anneal
    to that strand: it anneals to its complement, and `thermo.destabilisation` takes the
    complement itself for exactly this reason. The mismatch in the duplex is therefore the
    oligo base against the **complement** of the genome base.

    Getting this wrong is silent and it inverts the answer. A primer wanting `A` over a
    genome `G` reads as "A-G", which Stadhouders lists as severe; the duplex mismatch is
    A against C, which they list as minor. The first version of this script classified the
    genome-sense pair and reported 78% of 3'-proximal mismatches as severe. The giveaway
    was in its own output: "A-G" and "C-T" came back with identical mean destabilisation
    to four decimal places, because both are the same A-C duplex seen from opposite
    strands.

    `None` where either base is a degenerate code, since the mismatch is then not one pair.
    """
    if oligo_base not in COMPLEMENT or genome_base not in COMPLEMENT:
        return None
    return frozenset(oligo_base + COMPLEMENT[genome_base])


def classify(pair: frozenset[str] | None) -> str:
    """Stadhouders' grouping, with everything they did not name left unnamed."""
    if pair is None:
        return "degenerate"
    if pair in MINOR:
        return "minor"
    if pair in SEVERE:
        return "severe"
    return "unclassified"


def _label(pair: frozenset[str] | None) -> str:
    """The duplex mismatch, written as an unordered pair."""
    if pair is None:
        return "degenerate"
    if len(pair) == 1:
        base = next(iter(pair))
        return f"{base}·{base}"
    return "·".join(sorted(pair))


def scan(data: Path, limit: int = 0, checkpoint: Path | None = None, every: int = 20_000) -> dict:
    """Every definite mismatch under every oligo, by type and by 3' distance."""
    from analyse import load  # noqa: PLC0415

    from assaydrift.match import (  # noqa: PLC0415
        AMBIGUOUS,
        MATCH,
        THREE_PRIME_WINDOW,
        compare_base,
        find_oligo,
        reverse_complement,
    )
    from assaydrift.primers import ASSAYS  # noqa: PLC0415
    from assaydrift.thermo import destabilisation  # noqa: PLC0415

    # type -> counts, split by whether the mismatch is inside the 3' window of a primer.
    by_type: dict[str, dict] = {}
    # How many genome-oligo alignments carry exactly one 3'-proximal mismatch and nothing
    # else: the configuration §4.3's withdrawal is about.
    lone_three_prime: collections.Counter[str] = collections.Counter()
    lone_alignments = 0
    seen = skipped_ambiguous = 0

    for record in load(data):
        sequence = record.sequence
        if not sequence:
            continue
        seen += 1
        if limit and seen > limit:
            seen -= 1
            break
        for assay in ASSAYS:
            for oligo in assay.oligos:
                reverse = oligo.role == "reverse"
                wanted = (
                    reverse_complement(oligo.sequence.upper())
                    if reverse
                    else oligo.sequence.upper()
                )
                hit = find_oligo(sequence, oligo.sequence, oligo.role)
                if not hit.found or hit.indels:
                    continue
                window = sequence[hit.position : hit.position + len(wanted)].upper()
                if len(window) != len(wanted):
                    continue
                last = len(wanted) - 1
                found: list[tuple[frozenset[str], int]] = []
                unknown = False
                for index, (base, actual) in enumerate(zip(wanted, window, strict=True)):
                    verdict = compare_base(base, actual)
                    if verdict == AMBIGUOUS:
                        unknown = True
                        break
                    if verdict == MATCH:
                        continue
                    # The library's own distance rule, not a second copy of it.
                    distance = index if reverse else last - index
                    found.append((duplex_pair(base, actual), distance))
                if unknown:
                    skipped_ambiguous += 1
                    continue
                for pair, distance in found:
                    label = _label(pair)
                    row = by_type.setdefault(
                        label,
                        {
                            "type": label,
                            "stadhouders": classify(pair),
                            "total": 0,
                            "three_prime": 0,
                            "on_a_primer": 0,
                            "three_prime_on_a_primer": 0,
                        },
                    )
                    row["total"] += 1
                    inside = distance < THREE_PRIME_WINDOW
                    if inside:
                        row["three_prime"] += 1
                    if oligo.role != "probe":
                        row["on_a_primer"] += 1
                        if inside:
                            row["three_prime_on_a_primer"] += 1
                # The §4.3 configuration: exactly one mismatch, and it is 3'-proximal on a
                # primer. Counted per alignment rather than per mismatch, because that is
                # the unit the published rate counts.
                if oligo.role != "probe" and len(found) == 1:
                    pair, distance = found[0]
                    if distance < THREE_PRIME_WINDOW:
                        lone_alignments += 1
                        lone_three_prime[_label(pair)] += 1

        if checkpoint is not None and seen % every == 0:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text(
                json.dumps(
                    _payload(
                        data.name,
                        seen,
                        skipped_ambiguous,
                        by_type,
                        lone_three_prime,
                        lone_alignments,
                        destabilisation,
                        partial=True,
                    ),
                    indent=1,
                ),
                encoding="utf-8",
            )
        if seen % 5000 == 0:
            print(f"  {seen:,} genomes", flush=True)
    print(f"  {seen:,} genomes scanned", flush=True)

    return _payload(
        data.name,
        seen,
        skipped_ambiguous,
        by_type,
        lone_three_prime,
        lone_alignments,
        destabilisation,
        partial=False,
    )


def _payload(
    corpus: str,
    seen: int,
    skipped_ambiguous: int,
    by_type: dict,
    lone_three_prime,
    lone_alignments: int,
    destabilisation,
    *,
    partial: bool,
) -> dict:
    """Build the result from whatever has been counted so far.

    Split out so a checkpoint is the same shape as the finished artefact, flagged
    `partial`. One pass over 291,827 genomes is hours on a shared machine, and the first
    attempt at this scan died with its session having written nothing at all.
    """
    # The thermodynamic second opinion, computed once per type rather than per observation.
    #
    # Both arguments go in the *oligo's own sense*; `destabilisation` takes the complement
    # of the observed slice itself, which is the whole reason it is a function rather than
    # a caller's two lines. Complementing the flanks here would make every flank a mismatch
    # too, and the number would be the cost of three mismatches reported as the cost of one.
    #
    # The context is the whole difficulty. A nearest-neighbour mismatch cost depends on
    # *both* immediate neighbours and on nothing further away, so the honest figure is the
    # mean over all sixteen neighbour pairs rather than one arbitrary flank. The outer
    # `FLANK` exists only to push the mismatch away from the duplex ends: a three-base
    # duplex puts the mismatch one base from a terminus, where fraying dominates the
    # arithmetic, and the first version of this measurement did exactly that and reported
    # every type as costing about the same.
    for label, row in by_type.items():
        if label == "degenerate":
            row["mean_delta_delta_g"] = None
            row["spread_delta_delta_g"] = None
            row["contexts_scored"] = 0
            continue
        a, b = label.split("·")
        # `label` is the pair that faces each other in the duplex, so it has to be turned
        # back into the genome-sense base the function expects: `destabilisation` takes the
        # complement of its second argument. Passing the duplex partner directly makes a
        # self-pair - A against A - into two identical strings with no mismatch at all, and
        # the result is a flat 0.00 kcal/mol for exactly the types the literature calls the
        # most expensive. That is what the first version of this loop reported.
        costs = []
        for left in "ACGT":
            for right in "ACGT":
                oligo_strand = f"{FLANK}{left}{a}{right}{FLANK}"
                genome_strand = f"{FLANK}{left}{COMPLEMENT[b]}{right}{FLANK}"
                result = destabilisation(oligo_strand, genome_strand)
                if result is not None:
                    costs.append(result.delta_delta_g)
                # The mirror orientation too, because the pair is unordered here and the
                # model is not: an oligo `a` facing `b` is a different duplex from an
                # oligo `b` facing `a`.
                mirror = destabilisation(
                    f"{FLANK}{left}{b}{right}{FLANK}",
                    f"{FLANK}{left}{COMPLEMENT[a]}{right}{FLANK}",
                )
                if mirror is not None:
                    costs.append(mirror.delta_delta_g)
        row["mean_delta_delta_g"] = round(sum(costs) / len(costs), 4) if costs else None
        row["spread_delta_delta_g"] = round(max(costs) - min(costs), 4) if len(costs) > 1 else None
        row["contexts_scored"] = len(costs)

    ordered = sorted(by_type.values(), key=lambda r: -r["three_prime_on_a_primer"])
    totals = {
        "genomes": seen,
        "alignments_with_an_unknown_base": skipped_ambiguous,
        "mismatches": sum(r["total"] for r in ordered),
        "three_prime_on_a_primer": sum(r["three_prime_on_a_primer"] for r in ordered),
        "lone_three_prime_primer_alignments": lone_alignments,
    }
    for group in ("minor", "severe", "unclassified", "degenerate"):
        totals[f"three_prime_on_a_primer_{group}"] = sum(
            r["three_prime_on_a_primer"] for r in ordered if r["stadhouders"] == group
        )
        totals[f"lone_{group}"] = sum(
            n
            for label, n in lone_three_prime.items()
            if classify(None if label == "degenerate" else frozenset(label.split("·"))) == group
        )
    return {
        "corpus": corpus,
        "partial": partial,
        "source": (
            "Stadhouders et al. 2010, J Mol Diagn 12(1):109-117, abstract: single "
            "mismatches range from minor (<1.5 Ct; A-C, C-A, T-G, G-T) to severe "
            "(>7.0 Ct; A-A, G-A, A-G, C-C)."
        ),
        "totals": totals,
        "lone_three_prime_by_type": dict(lone_three_prime.most_common()),
        "by_type": ordered,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--out", type=Path, default=RESULTS / "mismatch_types.json")
    args = parser.parse_args()

    payload = scan(args.data, args.limit, args.out.with_suffix(".partial.json"))
    if payload["partial"]:
        raise SystemExit("the scan returned a partial payload, which should be impossible")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    t = payload["totals"]
    print(f"{t['genomes']:,} genomes, {t['mismatches']:,} definite mismatches")
    print(f"{t['three_prime_on_a_primer']:,} of them 3'-proximal on a primer\n")
    print(f"{'type':12}{'Stadhouders':14}{'3prime/primer':>14}{'share':>8}{'ddG':>9}")
    total = t["three_prime_on_a_primer"] or 1
    for row in payload["by_type"]:
        ddg = "n/a" if row["mean_delta_delta_g"] is None else f"{row['mean_delta_delta_g']:.2f}"
        print(
            f"{row['type']:12}{row['stadhouders']:14}{row['three_prime_on_a_primer']:>14,}"
            f"{row['three_prime_on_a_primer'] / total:>8.1%}{ddg:>9}"
        )
    print()
    for group in ("minor", "severe", "unclassified", "degenerate"):
        n = t[f"three_prime_on_a_primer_{group}"]
        print(f"  {group:14}{n:>9,}{n / total:>8.1%}")
    print(
        f"\nalignments whose only mismatch is 3'-proximal on a primer: "
        f"{t['lone_three_prime_primer_alignments']:,}"
    )
    for group in ("minor", "severe", "unclassified"):
        n = t[f"lone_{group}"]
        lone = t["lone_three_prime_primer_alignments"]
        share = n / lone if lone else 0.0
        print(f"  {group:14}{n:>9,}{share:>8.1%}")
    print(f"\nwritten to {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
