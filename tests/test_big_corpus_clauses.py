"""Section 4.5's numbers, read off `results/drift_big.json` rather than retyped.

Section 4.2's claim - no genome in the 2,765-genome corpus carries two mismatches in a
3' window - is pinned by `test_the_small_corpus_has_no_three_prime_pair_anywhere`. It is
a claim about a corpus, and the obvious way for it to be wrong is for a hundredfold
larger corpus to contain the configuration. It does: 52 assay-genome pairs. This file
pins the larger corpus's side of that, because the paper and the README both now publish
the counts and a count in prose drifts away from the artefact it came from.

The counts are the ones published in paper/main.md section 4.5 and in the README's
withdrawal table. A regenerated artefact that disagrees should fail here, loudly, rather
than leave two tables quietly wrong.
"""

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
BIG = ROOT / "results" / "drift_big.json"

# assay -> (usable, three_prime hits, of which a pair, severe_load, clause 1)
PUBLISHED = {
    "CDC N1": (289_287, 423, 4, 1_445, 0),
    "CDC N2": (288_108, 2_097, 0, 4, 0),
    "CDC N3": (289_945, 672, 0, 4, 0),
    "Charite E": (285_298, 98, 0, 20, 0),
    "Charite RdRp": (281_456, 51_465, 48, 4, 0),
}

PUBLISHED_PAIRS = 52
PUBLISHED_HITS = 54_755
PUBLISHED_EVALUATIONS = 1_434_094


def _totals():
    if not BIG.exists():
        pytest.skip("results/drift_big.json is not present")
    data = json.loads(BIG.read_text(encoding="utf-8"))
    out = {}
    for assay, body in data.items():
        if assay == "_run":
            continue
        agg = dict.fromkeys(
            (
                "usable",
                "three_prime",
                "three_prime_pair",
                "severe_load",
                "lost_oligo",
                "failing_clustered",
            ),
            0,
        )
        for s in body["periods"].values():
            if "three_prime_pair" not in s:
                pytest.skip("results/drift_big.json predates the clause counters")
            for key in agg:
                agg[key] += s.get(key) or 0
        out[assay] = agg
    return out


def test_the_published_per_assay_counts_match_the_artefact():
    totals = _totals()
    assert set(totals) == set(PUBLISHED), (set(totals), set(PUBLISHED))
    for assay, (usable, hits, pairs, severe, clause_one) in PUBLISHED.items():
        got = totals[assay]
        assert got["usable"] == usable, (assay, "usable", got["usable"], usable)
        assert got["three_prime"] == hits, (assay, "three_prime", got["three_prime"], hits)
        assert got["three_prime_pair"] == pairs, (
            assay,
            "three_prime_pair",
            got["three_prime_pair"],
            pairs,
        )
        assert got["severe_load"] == severe, (assay, "severe_load", got["severe_load"], severe)
        assert got["lost_oligo"] == clause_one, (
            assay,
            "lost_oligo",
            got["lost_oligo"],
            clause_one,
        )


def test_clause_one_fires_for_no_sequence_in_either_corpus():
    """The abstract's strongest sentence, over the larger corpus.

    `lost_oligo` is the one clause of `likely_failing` a bench would not dispute. The
    abstract says it fires for no sequence in either corpus. The one record that used to
    fire it - LR877184, a 2,859-base partial deposit - is counted as `region_absent`
    instead, so this assertion is only honest if that counter is non-zero: otherwise the
    zero could mean the predicate had been widened until it caught nothing.
    """
    totals = _totals()
    data = json.loads(BIG.read_text(encoding="utf-8"))
    for assay, agg in totals.items():
        assert agg["lost_oligo"] == 0, (assay, agg["lost_oligo"])
        absent = sum(s.get("region_absent") or 0 for s in data[assay]["periods"].values())
        assert absent == 1, (
            f"{assay}: region_absent={absent}, expected the single partial deposit. A zero "
            f"here means lost_oligo==0 was reached by excluding records, not by the rule."
        )


def test_the_three_prime_pair_exists_in_the_larger_corpus():
    """The point of section 4.5: the configuration section 4.2 found none of is reachable."""
    totals = _totals()
    pairs = sum(a["three_prime_pair"] for a in totals.values())
    hits = sum(a["three_prime"] for a in totals.values())
    evaluations = sum(a["usable"] for a in totals.values())
    assert pairs == PUBLISHED_PAIRS, (pairs, PUBLISHED_PAIRS)
    assert hits == PUBLISHED_HITS, (hits, PUBLISHED_HITS)
    assert evaluations == PUBLISHED_EVALUATIONS, (evaluations, PUBLISHED_EVALUATIONS)
    assert pairs > 0, "section 4.5 is written around pairs existing here"
    share = pairs / hits
    assert f"{share * 100:.3f}" == "0.095", f"{share * 100:.3f}% is published as 0.095%"


def test_two_assays_in_six_quarters_carry_a_pair():
    """The README says 'reached by two assays in six quarters'; count them."""
    if not BIG.exists():
        pytest.skip("results/drift_big.json is not present")
    data = json.loads(BIG.read_text(encoding="utf-8"))
    assays, quarters = set(), set()
    for assay, body in data.items():
        if assay == "_run":
            continue
        for period, s in body["periods"].items():
            if s.get("three_prime_pair"):
                assays.add(assay)
                quarters.add(period)
    assert len(assays) == 2, sorted(assays)
    assert len(quarters) == 6, sorted(quarters)


def test_the_rdrp_pairs_are_concentrated_in_the_delta_quarters():
    """Section 4.5 attributes 46 of RdRp's 48 pairs to 2021-Q3 and 2021-Q4."""
    if not BIG.exists():
        pytest.skip("results/drift_big.json is not present")
    data = json.loads(BIG.read_text(encoding="utf-8"))
    periods = data["Charite RdRp"]["periods"]
    delta = sum(
        periods[q].get("three_prime_pair") or 0 for q in ("2021-Q3", "2021-Q4") if q in periods
    )
    assert delta == 46, delta


def test_requiring_two_three_prime_mismatches_drops_rdrp_by_three_orders():
    """The withdrawal's arithmetic, at the published settings, on the larger corpus."""
    totals = _totals()
    rdrp = totals["Charite RdRp"]
    published = rdrp["three_prime"] / rdrp["usable"]
    # `failing_clustered`, not the two counters added: the rule is a union and the sum
    # double-counts a genome in both clauses - 1 of RdRp's at this corpus size. Both
    # round to 0.00018, which is why the published table is unchanged, and is also
    # exactly why nothing noticed. See test_clustered_is_the_union_not_the_sum.py.
    clustered = rdrp["failing_clustered"] / rdrp["usable"]
    assert f"{published:.4f}" == "0.1829", f"{published:.4f}"
    assert clustered < published / 500, (clustered, published)
    assert f"{clustered:.5f}" == "0.00018", f"{clustered:.5f}"


def test_the_readme_withdrawal_table_matches_the_artefact():
    """The README publishes five clustered rates on its front page; recompute all five.

    The README's other figures are checked by `test_published_numbers.py`, which was
    written before this table existed and does not know about it. A number on a front
    page that nothing recomputes is the kind that survives a regeneration unchanged.
    """
    totals = _totals()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    expected = {
        "Charite RdRp": "0.00018",
        "CDC N1": "0.0050",
        "CDC N2": "0.000014",
        "CDC N3": "0.000014",
        "Charite E": "0.00007",
    }
    for assay, shown in expected.items():
        agg = totals[assay]
        rate = agg["failing_clustered"] / agg["usable"]
        places = len(shown.split(".")[1])
        assert f"{rate:.{places}f}" == shown, (assay, f"{rate:.{places}f}", shown)
        assert f"| {shown} |" in readme or f"| **{shown}** |" in readme, (assay, shown)

    published = {
        "Charite RdRp": "0.1829",
        "CDC N1": "0.0062",
        "CDC N2": "0.0073",
        "CDC N3": "0.0023",
        "Charite E": "0.0004",
    }
    for assay, shown in published.items():
        agg = totals[assay]
        # clause 1 is zero everywhere, so the published rate is clause 2 at one
        # mismatch unioned with clause 3; the artefact's own `failing` count is what
        # the table reports and is read back here rather than reconstructed.
        data = json.loads(BIG.read_text(encoding="utf-8"))
        failing = sum(s.get("failing") or 0 for s in data[assay]["periods"].values())
        assert f"{failing / agg['usable']:.4f}" == shown, (
            assay,
            f"{failing / agg['usable']:.4f}",
            shown,
        )
        assert f"| {shown} |" in readme, (assay, shown)


def test_section_4_5_publishes_the_same_rates_as_the_readme():
    """Two tables, one artefact. The duplication is deliberate; the drift would not be."""
    paper = (ROOT / "paper" / "main.md").read_text(encoding="utf-8")
    start = paper.index("## 4.5 The same re-scoring")
    section = paper[start : paper.index("\n# 5. Threats to validity")]
    for shown in ("0.00018", "0.0050", "0.000014", "0.000070", "0.1829"):
        assert shown in section, shown
    for count in ("51,465", "289,287", "281,456", "54,755", "1,434,094"):
        assert count in section, count
    assert "52 pairs in 54,755" in section
