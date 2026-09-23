"""Show what this does, in about two seconds, with no network.

    python demo.py

Two halves. First the setup check: every published oligo is aligned against the
2019 reference genome, vendored in tests/data. If those do not match, nothing
downstream means anything, so it is the first thing shown rather than a detail
in a test file.

Then the finding, read from results/drift.json - the output of a real run over
genomes downloaded from GenBank. Rebuild it with:

    python scripts/fetch.py --per-year 400
    python scripts/analyse.py
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from assaydrift import primers  # noqa: E402
from assaydrift.match import find_oligo  # noqa: E402

ROOT = Path(__file__).resolve().parent
REFERENCE = ROOT / "tests" / "data" / "NC_045512.2.json.gz"
RESULTS = ROOT / "results" / "drift.json"


def show_reference_check() -> None:
    with gzip.open(REFERENCE, "rt", encoding="utf-8") as handle:
        genome = json.load(handle)["sequence"]

    print("=" * 78)
    print("SETUP CHECK - every published oligo against the 2019 reference genome")
    print("=" * 78)
    print(f"NC_045512.2, {len(genome):,} bases, the first SARS-CoV-2 genome deposited.")
    print("A primer that does not match this is a typo, not a discovery.\n")

    print(f"  {'oligo':<10}{'role':<9}{'len':>4}{'position':>10}{'mismatch':>10}   note")
    for _assay, oligo in primers.all_oligos():
        hit = find_oligo(genome, oligo.sequence, oligo.role)
        expected = primers.expected_reference_mismatches(oligo.name)
        status = "ok" if hit.mismatches == expected else "UNEXPECTED"
        note = ""
        if expected:
            known = primers.KNOWN_REFERENCE_MISMATCHES[oligo.name]
            note = (
                f"known: {known['oligo_base']} vs {known['reference_base']}, "
                f"{known['bases_from_3_prime_end']} bases from the 3' end"
            )
        print(
            f"  {oligo.name:<10}{oligo.role:<9}{len(oligo):>4}{hit.position:>10}"
            f"{hit.mismatches:>10}   {status if status != 'ok' else note}"
        )

    print("\n  14 of 15 match perfectly. The exception is real and documented:")
    print("  Charite RdRp's reverse primer was designed for SARS-related")
    print("  coronaviruses generally, from sequences that predate SARS-CoV-2, so")
    print("  it has never matched its target exactly. It sits mid-primer rather")
    print("  than at the 3' end, which is why the assay worked anyway.\n")


def show_drift() -> None:
    if not RESULTS.exists():
        print("No results yet. Build them with:")
        print("    python scripts/fetch.py --per-year 400")
        print("    python scripts/analyse.py")
        return

    payload = json.loads(RESULTS.read_text(encoding="utf-8"))
    print("=" * 78)
    print("THE FINDING - would each assay still detect what was circulating?")
    print("=" * 78)
    print("Share of sequences collected in each quarter that the assay should")
    print("still amplify. Sequences with unknown bases under an oligo are")
    print("excluded, never counted as failures - see match.py.\n")

    for name, data in payload.items():
        label = name
        if data["retired"]:
            label += "  [retired by CDC]"
        if data["baseline_mismatches"]:
            label += f"  [shipped with {data['baseline_mismatches']} mismatch]"
        print(label)
        movement = data["trend"]
        if movement.get("enough_to_say"):
            print(
                f"   {movement['first_period']}  {movement['first_rate']:>6.1%}"
                f"  (n={movement['first_n']})"
            )
            print(
                f"   {movement['last_period']}  {movement['last_rate']:>6.1%}"
                f"  (n={movement['last_n']})"
            )
            print(
                f"   change {movement['change']:+.1%}"
                f"    worst quarter {movement['worst_period']} "
                f"at {movement['worst_rate']:.1%}"
            )
        else:
            print("   not enough reliable periods to claim a trend")
        print()

    print("Full per-quarter tables: results/drift.json, or run scripts/analyse.py")


def main() -> int:
    show_reference_check()
    show_drift()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
