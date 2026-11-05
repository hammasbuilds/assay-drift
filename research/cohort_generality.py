"""Is section 4's finding about these five assays, or about the rule?

    .venv-research/Scripts/python research/cohort_generality.py [--data data/sequences.jsonl]

Section 4 takes apart `likely_failing` on five assays: the CDC panel, the two Charite
assays, and N3 as a positive control. Its central claim is structural rather than
assay-specific - that the published rate is produced almost entirely by clause 2 firing on
a *single* 3'-proximal mismatch, which is the one configuration the wet-lab literature says
costs about a cycle and no false negative. Five assays is a thin base for a claim about how
a family of rules behaves, and nothing in the paper tested it anywhere else.

This scores the same rule, unchanged, over eight further published protocols: Charite N,
HKU's ORF1b and N assays, China CDC's ORF1ab and N, Institut Pasteur's IP2 and IP4, and
Thailand NIH's N. They are in `primers.ADDITIONAL_ASSAYS`, kept apart from the five so that
no published per-assay number moves, and every one was admitted only after binding the 2019
reference exactly, in the region its name claims, with its probe inside its own amplicon.

The question is not whether the other assays drift - that is a fact about them. It is
whether the *decomposition* holds: of the genomes each assay flags, what share is flagged by
clause 2 alone, and of those, what share carries exactly one 3'-proximal mismatch rather
than a pair. If the answer is "nearly all, nearly everywhere", section 4 is about the rule.
If some assays are flagged mostly by mismatch load, the rule's weakness is narrower than
claimed and the paper should say which assays it applies to.

Both cohorts are scored in the same run so the comparison is like for like.
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


def score(data: Path, limit: int = 0) -> dict:
    from analyse import baseline_for, load  # noqa: PLC0415

    from assaydrift.analyze import (  # noqa: PLC0415
        CLUSTERED_THREE_PRIME,
        SEVERE_MISMATCH_LOAD,
        evaluate,
    )
    from assaydrift.primers import ADDITIONAL_ASSAYS, ASSAYS  # noqa: PLC0415

    cohorts = {"primary": ASSAYS, "additional": ADDITIONAL_ASSAYS}
    rows: dict[str, dict] = {}
    for cohort, assays in cohorts.items():
        for assay in assays:
            rows[assay.name] = {
                "assay": assay.name,
                "cohort": cohort,
                "target": assay.target_gene,
                "usable": 0,
                "failing": 0,
                "three_prime_any": 0,
                "three_prime_pair": 0,
                "lone_three_prime": 0,
                "severe_load": 0,
                "lost_oligo": 0,
                "clause_two_only": 0,
            }
    baselines = {assay.name: baseline_for(assay) for assays in cohorts.values() for assay in assays}

    seen = 0
    for record in load(data):
        if not record.sequence:
            continue
        seen += 1
        if limit and seen > limit:
            seen -= 1
            break
        # Thirteen assays per genome, and several of the cohort's oligos are short enough
        # that the exact-match path misses more often and the gapped aligner runs. The
        # first run of this script printed nothing for ten minutes, which is
        # indistinguishable from a hang; a line every 500 genomes is the difference
        # between waiting and guessing.
        if seen % 500 == 0:
            print(f"  {seen:,} genomes", flush=True)
        for assays in cohorts.values():
            for assay in assays:
                result = evaluate(record, assay)
                if result is None:
                    continue
                row = rows[assay.name]
                baseline = baselines[assay.name]
                if result.any_oligo_lost:
                    row["lost_oligo"] += 1
                if not result.usable:
                    continue
                row["usable"] += 1
                excess = result.total_mismatches - baseline
                heavy = excess >= SEVERE_MISMATCH_LOAD
                three_prime = result.three_prime_mismatches
                if heavy:
                    row["severe_load"] += 1
                if three_prime > 0:
                    row["three_prime_any"] += 1
                if three_prime >= CLUSTERED_THREE_PRIME:
                    row["three_prime_pair"] += 1
                if result.likely_failing(baseline):
                    row["failing"] += 1
                    # The configuration section 4 withdraws: clause 2 fired, clause 3 did
                    # not, and the 3' evidence is a single mismatch rather than a pair.
                    if three_prime > 0 and not heavy:
                        row["clause_two_only"] += 1
                        if three_prime == 1:
                            row["lone_three_prime"] += 1

    for row in rows.values():
        usable, failing = row["usable"], row["failing"]
        row["failing_rate"] = round(failing / usable, 6) if usable else None
        row["clause_two_share_of_failing"] = (
            round(row["clause_two_only"] / failing, 4) if failing else None
        )
        row["lone_share_of_failing"] = (
            round(row["lone_three_prime"] / failing, 4) if failing else None
        )

    ordered = sorted(rows.values(), key=lambda r: (r["cohort"], -(r["failing_rate"] or 0)))
    summary = {}
    for cohort in cohorts:
        members = [r for r in ordered if r["cohort"] == cohort and r["failing"]]
        summary[cohort] = {
            "assays": sum(1 for r in ordered if r["cohort"] == cohort),
            "assays_with_any_failure": len(members),
            "failing": sum(r["failing"] for r in members),
            "clause_two_only": sum(r["clause_two_only"] for r in members),
            "lone_three_prime": sum(r["lone_three_prime"] for r in members),
            "severe_load": sum(r["severe_load"] for r in members),
            "lost_oligo": sum(r["lost_oligo"] for r in members),
        }
        total = summary[cohort]["failing"]
        summary[cohort]["clause_two_share"] = (
            round(summary[cohort]["clause_two_only"] / total, 4) if total else None
        )
        summary[cohort]["lone_share"] = (
            round(summary[cohort]["lone_three_prime"] / total, 4) if total else None
        )
    return {
        "corpus": data.name,
        "genomes": seen,
        "what": (
            "The published rule, unchanged, over both cohorts. `clause_two_only` counts "
            "genomes a rule flags on the 3' clause alone, with no severe mismatch load; "
            "`lone_three_prime` is the subset carrying exactly one 3'-proximal mismatch, "
            "which is the configuration section 4.3 withdraws."
        ),
        "summary": summary,
        "per_assay": ordered,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data" / "sequences.jsonl")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--out", type=Path, default=RESULTS / "cohort_generality.json")
    args = parser.parse_args()

    payload = score(args.data, args.limit)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    print(f"{payload['genomes']:,} genomes\n")
    print(
        f"{'assay':20}{'cohort':12}{'usable':>8}{'failing':>9}{'rate':>9}{'clause2':>9}{'lone':>8}"
    )
    for row in payload["per_assay"]:
        rate = "-" if row["failing_rate"] is None else f"{row['failing_rate']:.4f}"
        c2 = (
            "-"
            if row["clause_two_share_of_failing"] is None
            else f"{row['clause_two_share_of_failing']:.3f}"
        )
        lone = (
            "-" if row["lone_share_of_failing"] is None else f"{row['lone_share_of_failing']:.3f}"
        )
        print(
            f"{row['assay']:20}{row['cohort']:12}{row['usable']:>8,}{row['failing']:>9,}"
            f"{rate:>9}{c2:>9}{lone:>8}"
        )
    print()
    for cohort, s in payload["summary"].items():
        # Formatted outside the f-string: a nested same-quote f-string is a syntax error
        # on the interpreter this repository targets, and ruff caught exactly that here
        # once before in research/grid_corpora.py.
        clause_two = "n/a" if s["clause_two_share"] is None else f"{s['clause_two_share']:.1%}"
        lone = "n/a" if s["lone_share"] is None else f"{s['lone_share']:.1%}"
        print(
            f"  {cohort:12} {s['assays_with_any_failure']}/{s['assays']} assays flag "
            f"anything; {s['failing']:,} flagged, clause 2 alone {clause_two}, "
            f"of which a lone 3' mismatch {lone}"
        )
    print(f"\nwritten to {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
