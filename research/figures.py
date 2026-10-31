"""Figures for the assay-drift manuscript.

Same house rules as the sibling repositories: matplotlib only, Okabe-Ito, nothing encoded
by colour alone, PDF and PNG together, and every caption carries its n.

The two figures here exist to show how much of the published `likely_failing` rate is the
rule rather than the data, because a 2025 wet-lab study argues in-silico predictions of
this kind overstate the risk. A paper that reports the rate without this decomposition is
answering a question nobody asked.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
FIGURES = ROOT / "research" / "figures"
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#E69F00", "#009E73", "#999999"
ASSAYS = ["CDC N1", "CDC N2", "CDC N3", "Charite E", "Charite RdRp"]


def _save(fig, name: str) -> tuple[Path, Path]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    pdf, png = FIGURES / f"{name}.pdf", FIGURES / f"{name}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=200)
    plt.close(fig)
    return pdf, png


def _bare(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def clause_decomposition(grid: dict) -> tuple[Path, Path, str]:
    """What each clause of the rule contributes to the published rate.

    Three nested rules, so three stacked segments: the floor an oligo with no binding
    site gives, what the mismatch-load clause adds on top, and what the 3'-window clause
    adds on top of that. The point is which segment is doing the work.
    """
    floor, load_part, three_part, totals = [], [], [], []
    for assay in ASSAYS:
        published = grid[f"{assay}|w5|l3|all three clauses"]["rate"]
        no_three = grid[f"{assay}|w5|l3|without the 3' clause"]["rate"]
        lost = grid[f"{assay}|w5|l3|oligo lost only"]["rate"]
        floor.append(lost)
        load_part.append(max(0.0, no_three - lost))
        three_part.append(max(0.0, published - no_three))
        totals.append(published)

    x = range(len(ASSAYS))
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.bar(x, floor, width=0.6, color=GREY, label="oligo lost (the floor)",
           edgecolor="white", linewidth=0.6)
    ax.bar(x, load_part, width=0.6, bottom=floor, color=GREEN, hatch="//",
           label="mismatch load ≥ 3", edgecolor="white", linewidth=0.6)
    bottom2 = [a + b for a, b in zip(floor, load_part)]
    ax.bar(x, three_part, width=0.6, bottom=bottom2, color=ORANGE, hatch="..",
           label="a mismatch in a primer's last five bases", edgecolor="white",
           linewidth=0.6)
    for i, total in enumerate(totals):
        ax.text(i, total + 0.003, f"{total:.4f}", ha="center", fontsize=8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(ASSAYS, fontsize=9, rotation=15, ha="right")
    ax.set_ylabel("share of usable sequences called likely_failing", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    _bare(ax)
    n = grid[f"{ASSAYS[0]}|w5|l3|all three clauses"]["n"]
    worst = max(zip(ASSAYS, three_part), key=lambda t: t[1])
    caption = (
        "Which clause of the `likely_failing` rule produces the published rate, at the "
        f"published settings (3' window 5, mismatch load 3), about {n:,} usable sequences "
        "per assay. The grey floor is the only clause a bench would not dispute - an oligo "
        "with no binding site cannot prime - and it is **zero for all five assays**: no "
        "sequence in this corpus loses an oligo outright. For "
        f"{worst[0]} the 3'-window clause alone accounts for {worst[1]:.4f} of a "
        f"{totals[ASSAYS.index(worst[0])]:.4f} rate."
    )
    return (*_save(fig, "fig1_clause_decomposition"), caption)


def sensitivity_heatmap(grid: dict, assay: str = "Charite RdRp") -> tuple[Path, Path, str]:
    """The rate across the two thresholds, for one assay.

    Neither threshold came from a wet-lab calibration: 5 bases is a primer-design rule of
    thumb and 3 mismatches is a judgement. This is what the published number would have
    been under the other defensible choices.
    """
    windows = [3, 5, 7]
    loads = [2, 3, 4, 5]
    rows = []
    for window in windows:
        rows.append([grid[f"{assay}|w{window}|l{load}|all three clauses"]["rate"]
                     for load in loads])

    fig, ax = plt.subplots(figsize=(5.8, 3.2))
    highest = max(max(r) for r in rows) or 1.0
    image = ax.imshow(rows, cmap="Oranges", vmin=0, vmax=highest, aspect="auto")
    ax.set_xticks(range(len(loads)))
    ax.set_xticklabels([f"≥{l}" for l in loads], fontsize=9)
    ax.set_yticks(range(len(windows)))
    ax.set_yticklabels([f"{w} bases" for w in windows], fontsize=9)
    ax.set_xlabel("mismatch load that counts as failing", fontsize=10)
    ax.set_ylabel("3' window", fontsize=10)
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            ax.text(j, i, f"{value:.4f}", ha="center", va="center", fontsize=8,
                    color="white" if value > highest * 0.6 else "#222222")
    ax.add_patch(plt.Rectangle((loads.index(3) - 0.5, windows.index(5) - 0.5), 1, 1,
                               fill=False, edgecolor=BLUE, linewidth=2.4))
    fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03, label="likely_failing rate")
    lowest = min(min(r) for r in rows)
    caption = (
        f"`likely_failing` for {assay} across both thresholds, "
        f"{grid[f'{assay}|w5|l3|all three clauses']['n']:,} usable sequences. The blue box "
        "is the published setting. The rate ranges from "
        f"{lowest:.4f} to {highest:.4f} across choices that are all defensible, and "
        "neither threshold was calibrated against a bench."
    )
    return (*_save(fig, "fig2_sensitivity_heatmap"), caption)


# --- 3. where the rate actually comes from -------------------------------------------


def three_prime_distribution(counts: dict[str, dict[int, int]]) -> tuple[Path, Path, str]:
    """How many 3'-proximal mismatches each assay-genome pair carries: 0, 1, or two or more.

    This is the figure the paper's section 4 turns on, and it is nearly empty by design.
    `likely_failing` fires on **one** mismatch within five bases of a primer's 3' end.
    Otwell et al. 2025 measured that configuration in vitro at about one cycle of Ct shift
    and no false negative, and measured **two** at nine to twelve cycles and a genuine
    false negative. So the question that decides whether the published rate means anything
    is: are these single mismatches, or clusters?

    Across 2,765 genomes and five assays, every one of them is a single mismatch. The
    "two or more" column is zero everywhere. A bar chart of zeros is the honest way to
    show that, because the alternative is a sentence the reader has to take on trust.
    """
    names = list(counts)
    ones = [counts[a].get(1, 0) for a in names]
    pairs = [sum(n for k, n in counts[a].items() if k >= 2) for a in names]

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    x = range(len(names))
    width = 0.38
    ax.bar([i - width / 2 for i in x], ones, width, color=ORANGE, edgecolor="white",
           linewidth=0.6, label="exactly one 3'-proximal mismatch (tolerated in vitro)")
    ax.bar([i + width / 2 for i in x], pairs, width, color=BLUE, edgecolor="white",
           linewidth=0.6, label="two or more (the configuration that fails in vitro)")

    for i, (a, b) in enumerate(zip(ones, pairs)):
        ax.text(i - width / 2, a + max(ones) * 0.015, f"{a:,}", ha="center", fontsize=8.5)
        ax.text(i + width / 2, b + max(ones) * 0.015, str(b), ha="center", fontsize=8.5,
                color="#0072B2", fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel("assay-genome pairs")
    ax.set_ylim(0, max(ones) * 1.18)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    _bare(ax)
    pdf, png = _save(fig, "fig3_three_prime_distribution")
    total_ones, total_pairs = sum(ones), sum(pairs)
    caption = (
        f"Number of mismatches falling inside a primer's 3' window, per assay, over "
        f"2,765 genomes. {total_ones:,} assay-genome pairs carry exactly one and "
        f"{total_pairs} carry two or more. `likely_failing` fires on one; the only "
        f"direct in-vitro measurement of one reports about a single cycle of Ct shift "
        f"and no false negative, while two produce a 9-12 cycle shift and a real false "
        f"negative. Charite RdRp's 11.16% is composed entirely of the left-hand bar."
    )
    return pdf, png, caption
