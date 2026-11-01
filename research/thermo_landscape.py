"""Does the 3' positional rule track binding energy? Measured, not assumed.

Section 4 is about `likely_failing` resting on `THREE_PRIME_WINDOW = 5`: a mismatch in the
last five bases of a primer counts as stopping the reaction. That is a primer-design rule
of thumb. `assaydrift.thermo` computes a quantity that is not a rule of thumb - the
nearest-neighbour free-energy cost of the observed target relative to a perfect match -
so the two can be compared, which is a thing section 4 could not previously do.

Two measurements, and they answer different questions.

**The landscape** enumerates every single substitution in every oligo of the catalogue -
fifteen oligos, three alternatives at each base - and asks whether distance from the 3' end
predicts dGG. No corpus is involved: this is a property of the published oligos and the
published parameter set, so it is exhaustive rather than sampled, and nothing about which
mutations happened to circulate can bias it.

**The observed mutations** take the named mutations already published in
`results/mutations_big.json` and score each one, so the landscape can be weighted by what
is actually out there rather than by what is possible.

The expected finding, stated before looking, was that dGG would rise towards the 3' end and
so corroborate the rule. Write down the prediction and the measurement separately; this
file prints both.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from assaydrift import primers, thermo  # noqa: E402
from assaydrift.match import THREE_PRIME_WINDOW, reverse_complement  # noqa: E402

RESULTS = ROOT / "research" / "results"
MUTATIONS = ROOT / "results" / "mutations_big.json"

# Written down before the measurement, so that the result cannot be read as whatever it
# turned out to be. The 3' rule exists because polymerase extends from that terminus; if
# the rule is a proxy for binding energy then dGG should be highest where the rule fires.
PREDICTION = (
    "dGG rises towards the 3' terminus, so the positional window and the thermodynamic "
    "cost agree about which mismatches matter."
)

_ALTERNATIVES = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}

# Both published protocols are ordinary TaqMan chemistry with a combined anneal/extend
# step in the high fifties. Which end of that range is assumed decides the answer for
# several oligos, so both ends are reported and neither one is called "the" margin. An
# oligo that clears the low end and not the high one is a result about the protocol
# rather than about the genome, and collapsing the range to a single threshold would have
# resolved that by itself, silently, in whichever direction the threshold was picked.
ANNEAL_FLOOR_C = 55.0
ANNEAL_CEILING_C = 60.0


def _spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation, with ties averaged.

    Spearman rather than Pearson because the question is whether the two orderings agree,
    and neither quantity is expected to be linear in the other.
    """
    if len(xs) != len(ys) or len(xs) < 3:
        return None

    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            shared = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k]] = shared
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def _auc(flagged: list[float], unflagged: list[float]) -> float | None:
    """P(a flagged mismatch costs more than an unflagged one), ties counted as half.

    This is the question "does the rule pick out the expensive mismatches?" in one
    number. 0.5 is a rule that carries no information about binding energy at all.
    """
    if not flagged or not unflagged:
        return None
    wins = 0.0
    for a in flagged:
        for b in unflagged:
            wins += 1.0 if a > b else 0.5 if a == b else 0.0
    return wins / (len(flagged) * len(unflagged))


def landscape() -> dict:
    """Every single substitution in every catalogue oligo, scored."""
    rows = []
    unscorable = []
    for assay in primers.ASSAYS:
        for oligo in assay.oligos:
            sequence = oligo.sequence.upper()
            pool = thermo.resolutions(sequence)
            length = len(sequence)
            for index in range(length):
                # The oligo's own 5'->3' numbering, which is what `mutations.py`
                # reports, so the two halves of this file speak one coordinate system.
                from_three_prime = length - 1 - index
                reference = (pool[0] if pool else sequence)[index]
                for alternative in _ALTERNATIVES.get(reference, ""):
                    observed = (
                        (pool[0] if pool else sequence)[:index]
                        + alternative
                        + (pool[0] if pool else sequence)[index + 1 :]
                    )
                    cost = thermo.destabilisation(sequence, observed)
                    if cost is None:
                        unscorable.append({"oligo": oligo.name, "base_in_oligo": index + 1})
                        continue
                    rows.append(
                        {
                            "assay": assay.name,
                            "oligo": oligo.name,
                            "role": oligo.role,
                            "base_in_oligo": index + 1,
                            "bases_from_3_prime_end": from_three_prime,
                            "substitution": f"{reference}->{alternative}",
                            "delta_delta_g": round(cost.delta_delta_g, 4),
                            "delta_tm": round(cost.delta_tm, 4),
                            # The rule as the code applies it: primers only, and a probe
                            # 3' end is chemically blocked so the window is meaningless
                            # there. Getting this wrong is what once made N1 look dead.
                            "in_three_prime_window": (
                                oligo.role != "probe" and from_three_prime < THREE_PRIME_WINDOW
                            ),
                        }
                    )
    return {"rows": rows, "unscorable": unscorable}


def by_distance(rows: list[dict]) -> dict:
    """dGG summarised at each distance from the 3' terminus, primers only."""
    out = {}
    primer_rows = [r for r in rows if r["role"] != "probe"]
    for distance in range(0, 10):
        values = [
            r["delta_delta_g"] for r in primer_rows if r["bases_from_3_prime_end"] == distance
        ]
        if not values:
            continue
        out[str(distance)] = {
            "n": len(values),
            "mean": round(statistics.fmean(values), 3),
            "median": round(statistics.median(values), 3),
            "min": round(min(values), 3),
            "max": round(max(values), 3),
        }
    return out


def observed_mutations(rows_by_oligo: dict[str, dict]) -> dict:
    """The named mutations of `results/mutations_big.json`, each scored."""
    if not MUTATIONS.exists():
        return {"available": False, "rows": []}
    catalogue = {oligo.name: oligo for assay in primers.ASSAYS for oligo in assay.oligos}
    scored = []
    for assay_name, named in json.loads(MUTATIONS.read_text(encoding="utf-8")).items():
        for mutation in named:
            oligo = catalogue.get(mutation["oligo"])
            if oligo is None:
                continue
            name = mutation["name"]
            # e.g. "C28311T": reference base, genome position, observed base. The
            # bases are genome-sense, and a reverse primer's own sense is the other
            # strand, so the substitution has to be flipped before it is applied to
            # the oligo. Scoring the genome base against a reverse primer would
            # silently measure a different mutation.
            genome_observed = name[-1]
            observed_base = (
                reverse_complement(genome_observed) if oligo.role == "reverse" else genome_observed
            )
            sequence = oligo.sequence.upper()
            pool = thermo.resolutions(sequence)
            index = mutation["base_in_oligo"] - 1
            if not pool or not 0 <= index < len(sequence):
                continue
            base = pool[0]
            observed = base[:index] + observed_base + base[index + 1 :]
            cost = thermo.destabilisation(sequence, observed)
            scored.append(
                {
                    "assay": assay_name,
                    "mutation": name,
                    "oligo": mutation["oligo"],
                    "role": mutation["role"],
                    "bases_from_3_prime_end": mutation["bases_from_3_prime_end"],
                    "in_three_prime_window": mutation["in_three_prime_window"],
                    "overall_rate": mutation["overall_rate"],
                    "total": mutation["total"],
                    "delta_delta_g": round(cost.delta_delta_g, 4) if cost else None,
                    "delta_tm": round(cost.delta_tm, 4) if cost else None,
                    "perfect_tm": round(cost.perfect_tm, 2) if cost else None,
                    "observed_tm": round(cost.observed_tm, 2) if cost else None,
                    # A degraded Tm only matters if it falls towards the temperature the
                    # reaction anneals at, so the margin is the number to read, not the
                    # shift. A large shift that still clears the anneal by twenty degrees
                    # is a loss of headroom, not a failure, and saying otherwise would be
                    # the same overclaim in the alarming direction that section 4 is about.
                    "margin_over_anneal_floor_c": (
                        round(cost.observed_tm - ANNEAL_FLOOR_C, 2) if cost else None
                    ),
                    "margin_over_anneal_ceiling_c": (
                        round(cost.observed_tm - ANNEAL_CEILING_C, 2) if cost else None
                    ),
                    # Three outcomes, not two. "Clears the whole range" and "clears none
                    # of it" are findings about the oligo; "depends on where in the range
                    # the laboratory ran" is a finding about the protocol, and it is the
                    # one a single threshold would have silently resolved by itself.
                    "anneal_verdict": (
                        (
                            "clears"
                            if cost.observed_tm > ANNEAL_CEILING_C
                            else "below"
                            if cost.observed_tm < ANNEAL_FLOOR_C
                            else "depends_on_anneal_temperature"
                        )
                        if cost
                        else None
                    ),
                }
            )
    scored.sort(key=lambda r: -(r["delta_delta_g"] or 0.0))
    return {"available": True, "rows": scored}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=pathlib.Path, default=RESULTS / "thermo_landscape.json")
    args = parser.parse_args()

    whole = landscape()
    rows = whole["rows"]
    primer_rows = [r for r in rows if r["role"] != "probe"]
    flagged = [r["delta_delta_g"] for r in primer_rows if r["in_three_prime_window"]]
    unflagged = [r["delta_delta_g"] for r in primer_rows if not r["in_three_prime_window"]]
    terminal = [r["delta_delta_g"] for r in primer_rows if r["bases_from_3_prime_end"] == 0]
    inner_window = [
        r["delta_delta_g"]
        for r in primer_rows
        if 1 <= r["bases_from_3_prime_end"] < THREE_PRIME_WINDOW
    ]

    payload = {
        "what": (
            "Every single substitution in the 15 published oligos, scored for "
            "nearest-neighbour destabilisation, against the 3' positional rule that "
            "section 4 says is uncalibrated. Exhaustive over the oligos; no corpus."
        ),
        "prediction_written_before_measuring": PREDICTION,
        "three_prime_window": THREE_PRIME_WINDOW,
        "substitutions": len(rows),
        "unscorable": whole["unscorable"],
        "primer_substitutions": len(primer_rows),
        "flagged_by_the_positional_rule": len(flagged),
        "delta_delta_g": {
            "flagged_mean": round(statistics.fmean(flagged), 3) if flagged else None,
            "unflagged_mean": round(statistics.fmean(unflagged), 3) if unflagged else None,
            "terminal_base_mean": round(statistics.fmean(terminal), 3) if terminal else None,
            "window_but_not_terminal_mean": (
                round(statistics.fmean(inner_window), 3) if inner_window else None
            ),
            "terminal_base_max": round(max(terminal), 3) if terminal else None,
            "window_but_not_terminal_min": (round(min(inner_window), 3) if inner_window else None),
        },
        "auc_rule_vs_energy": (
            round(a, 4) if (a := _auc(flagged, unflagged)) is not None else None
        ),
        "spearman_distance_vs_delta_delta_g": (
            round(s, 4)
            if (
                s := _spearman(
                    [float(r["bases_from_3_prime_end"]) for r in primer_rows],
                    [r["delta_delta_g"] for r in primer_rows],
                )
            )
            is not None
            else None
        ),
        "by_distance_from_3_prime_end": by_distance(rows),
        "observed_named_mutations": observed_mutations({}),
        "rows": rows,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    summary = {k: v for k, v in payload.items() if k not in {"rows", "observed_named_mutations"}}
    print(json.dumps(summary, indent=1))
    print("\nTop observed mutations by binding cost:")
    for row in payload["observed_named_mutations"]["rows"][:8]:
        print(
            f"  {row['mutation']:9} {row['oligo']:7} {row['role']:8}"
            f" {row['bases_from_3_prime_end']:3} from 3'"
            f"  window={str(row['in_three_prime_window']):5}"
            f"  ddG={row['delta_delta_g']:5.2f}"
            f"  rate={row['overall_rate'] * 100:5.2f}%"
        )
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
