"""The nearest-neighbour model, against an independent implementation and its own claims.

`assaydrift.thermo` vendors four papers' worth of measured parameters and walks them. Two
things can go wrong and neither shows up as an exception: a parameter table with a wrong
digit in it, and a walk that scores the right duplex the wrong way round. A reversal is the
worst of them, because a palindromic stack would still come out right and the error would
only bite on real oligos.

So the suite checks the walk against Biopython's `MeltingTemp.Tm_NN`, which is a separate
implementation of the same published model, over every single substitution at every
position of randomly drawn oligos. That is the guarantee the module's docstring claims.
Biopython is a test-only dependency: where it is absent the oracle tests skip and the
invariant tests - which need nothing but the module - still run.

The invariants are the claims the paper leans on, written as tests so that a later edit to
the tables or the walk cannot quietly move them:

* a perfect match costs nothing
* a mismatch at the 3' terminus costs far less than the same mismatch one base in
* an unknown target base is refused rather than guessed
* a degenerate oligo is scored as the pool member that binds what is there
"""

from __future__ import annotations

import json
import pathlib
import random
import statistics
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from assaydrift import primers, thermo  # noqa: E402

LANDSCAPE = ROOT / "research" / "results" / "thermo_landscape.json"

# The reaction `thermo`'s defaults describe, spelled out so the oracle is called with the
# same conditions rather than with Biopython's own defaults.
CONDITIONS = {
    "saltcorr": 5,
    "Na": thermo.DEFAULT_NA_MM,
    "K": 0,
    "Tris": 0,
    "Mg": thermo.DEFAULT_MG_MM,
    "dNTPs": thermo.DEFAULT_DNTP_MM,
    "dnac1": thermo.DEFAULT_OLIGO_NM,
    "dnac2": 0,
}

# Fixed so a failure is reproducible. A seed that changes per run turns a real
# disagreement into one that cannot be bisected.
SEED = 20261101


def _oracle():
    return pytest.importorskip(
        "Bio.SeqUtils.MeltingTemp",
        reason="Biopython is the cross-check oracle and is test-only",
    )


def _catalogue():
    return [oligo for assay in primers.ASSAYS for oligo in assay.oligos]


# --------------------------------------------------------------------------- the oracle


def test_every_single_substitution_agrees_with_biopython():
    """The guarantee the module's docstring makes, over tens of thousands of duplexes."""
    oracle = _oracle()
    rng = random.Random(SEED)
    worst = 0.0
    worst_case = None
    checked = 0
    for _ in range(60):
        length = rng.randint(18, 26)
        oligo = "".join(rng.choice("ACGT") for _ in range(length))
        for index in range(length):
            for base in "ACGT":
                observed = oligo[:index] + base + oligo[index + 1 :]
                template = thermo.complement(observed)
                mine = thermo.melting_temperature(oligo, template)
                theirs = oracle.Tm_NN(oligo, c_seq=template, **CONDITIONS)
                checked += 1
                gap = abs(mine - theirs)
                if gap > worst:
                    worst, worst_case = gap, (oligo, observed, mine, theirs)
    assert checked > 4_000, checked
    assert worst < 1e-9, (worst, worst_case)


def test_the_published_oligos_agree_with_biopython():
    """The random draw above cannot contain a degenerate base; the catalogue does."""
    oracle = _oracle()
    checked = 0
    for oligo in _catalogue():
        pool = thermo.resolutions(oligo.sequence)
        assert pool, oligo.name
        for member in pool:
            template = thermo.complement(member)
            mine = thermo.melting_temperature(member, template)
            theirs = oracle.Tm_NN(member, c_seq=template, **CONDITIONS)
            assert abs(mine - theirs) < 1e-9, (oligo.name, member, mine, theirs)
            checked += 1
    assert checked >= len(_catalogue())


def test_multiple_mismatches_agree_with_biopython_or_are_refused_by_both():
    """Single substitutions could agree while the zipping mis-handles adjacent ones.

    Adjacent mismatches make a stack the published tables do not cover, and the oracle
    raises on exactly those. So the agreement checked here is two-sided: where the
    oracle produces a number, the numbers match; where it refuses, this module refuses
    too. A module that returned a value for a duplex the model cannot score would pass
    a one-sided test and be inventing thermodynamics.
    """
    oracle = _oracle()
    rng = random.Random(SEED + 1)
    agreed = 0
    both_refused = 0
    for _ in range(400):
        length = rng.randint(18, 26)
        oligo = "".join(rng.choice("ACGT") for _ in range(length))
        observed = list(oligo)
        for index in rng.sample(range(length), rng.randint(2, 5)):
            observed[index] = rng.choice("ACGT")
        template = thermo.complement("".join(observed))
        mine = thermo.melting_temperature(oligo, template)
        try:
            theirs = oracle.Tm_NN(oligo, c_seq=template, **CONDITIONS)
        except ValueError:
            assert mine is None, (oligo, "".join(observed), mine)
            both_refused += 1
            continue
        assert mine is not None, (oligo, "".join(observed), theirs)
        assert abs(mine - theirs) < 1e-9, (oligo, "".join(observed), mine, theirs)
        agreed += 1
    # Both arms have to have been exercised, or this test proves only one of them.
    assert agreed > 0 and both_refused > 0, (agreed, both_refused)


def test_a_tandem_mismatch_is_refused_rather_than_scored():
    """The published parameters stop at one mismatch between matched pairs.

    This branch was documented as unreachable and is not: two adjacent mispaired
    positions make a stack no table covers. The corpus contains the case - C28311T and
    C28312T are adjacent bases of the same N1 probe.
    """
    oligo = "ACGTACGTACGTACGTACGT"
    # One mismatch in the middle: scorable.
    assert thermo.destabilisation(oligo, "ACGTACGTAGGTACGTACGT") is not None
    # Two adjacent: not.
    assert thermo.destabilisation(oligo, "ACGTACGTAGCTACGTACGT") is None
    # Directly at the duplex level. `TACA` under `ACGT` mismatches at one position only
    # and is scorable; `TAAA` mismatches at two adjacent ones and is not.
    assert thermo.duplex("ACGT", "TACA") is not None
    assert thermo.duplex("ACGT", "TAAA") is None


def test_adjacent_named_mutations_exist_in_the_corpus():
    """The refusal above is a real limitation, not a hypothetical one."""
    data = _landscape()["observed_named_mutations"]
    if not data["available"]:
        pytest.skip("results/mutations_big.json is not present")
    positions = {}
    for row in data["rows"]:
        site = int("".join(c for c in row["mutation"] if c.isdigit()))
        positions.setdefault(row["oligo"], set()).add(site)
    adjacent = {
        oligo: sorted(sites)
        for oligo, sites in positions.items()
        if any(site + 1 in sites for site in sites)
    }
    assert adjacent, positions
    assert 28311 in adjacent.get("N1-P", []) and 28312 in adjacent.get("N1-P", [])


# ------------------------------------------------------------------------- the invariants


def test_a_perfect_match_costs_nothing():
    for oligo in _catalogue():
        pool = thermo.resolutions(oligo.sequence)
        for member in pool:
            cost = thermo.destabilisation(oligo.sequence, member)
            assert cost is not None, (oligo.name, member)
            assert abs(cost.delta_delta_g) < 1e-9, (oligo.name, member, cost)
            assert abs(cost.delta_tm) < 1e-9, (oligo.name, member, cost)
            assert cost.costs_nothing


def test_the_complement_does_not_reverse_the_strand():
    """A reversal here would score every duplex backwards and still look plausible."""
    assert thermo.complement("AAAC") == "TTTG"
    assert thermo.complement("AAAC") != "GTTT"
    assert thermo.complement("ACGT") == "TGCA"
    assert thermo.complement("N") is None
    assert thermo.complement("ACGU") is None


def test_a_mismatch_at_the_three_prime_terminus_costs_far_less_than_one_base_in():
    """The mechanism the paper's new section rests on, over the real catalogue.

    A terminal pair frays and gives up one stack; one base in, the same substitution
    breaks two. If this ever stopped holding, the section 4.6 argument would be wrong.
    """
    terminal = []
    adjacent = []
    for oligo in _catalogue():
        pool = thermo.resolutions(oligo.sequence)
        member = pool[0]
        for offset, bucket in ((0, terminal), (1, adjacent)):
            index = len(member) - 1 - offset
            base = member[index]
            for alternative in "ACGT":
                if alternative == base:
                    continue
                observed = member[:index] + alternative + member[index + 1 :]
                cost = thermo.destabilisation(oligo.sequence, observed)
                assert cost is not None
                bucket.append(cost.delta_delta_g)
    assert statistics.fmean(terminal) < statistics.fmean(adjacent) / 2.0, (
        statistics.fmean(terminal),
        statistics.fmean(adjacent),
    )


def test_an_unknown_target_base_is_refused_rather_than_guessed():
    oligo = _catalogue()[0].sequence
    for bad in ("N", "R", "-", "?"):
        observed = bad + oligo[1:]
        assert thermo.destabilisation(oligo, observed) is None, bad


def test_a_length_disagreement_is_refused():
    oligo = _catalogue()[0].sequence
    assert thermo.destabilisation(oligo, oligo[:-1]) is None
    assert thermo.destabilisation(oligo, oligo + "A") is None
    assert thermo.duplex("AC", "GTG") is None
    assert thermo.duplex("A", "T") is None


def test_an_oversized_degenerate_pool_is_refused_rather_than_truncated():
    """A truncated pool would report the best of an arbitrary subset."""
    assert thermo.resolutions("N" * 4) is None  # not a code this accepts at all
    small = thermo.resolutions("R" * 6)
    assert small is not None and len(small) == 64
    assert thermo.resolutions("R" * 7) is None
    assert thermo.MAX_RESOLUTIONS == 64


def test_a_degenerate_oligo_is_scored_as_the_member_that_binds_what_is_there():
    """The pool is a mixture the manufacturer synthesises, not an unknown."""
    oligo = "GGACTAGGCATTAGCATTR"
    pool = thermo.resolutions(oligo)
    assert pool is not None and len(pool) == 2
    for member in pool:
        cost = thermo.destabilisation(oligo, member)
        assert cost is not None
        assert cost.resolution == member
        assert cost.costs_nothing, (member, cost.delta_delta_g)


def test_a_mismatch_lowers_the_melting_temperature():
    """Direction, on every oligo: the sign of delta_tm is load-bearing for the paper."""
    for oligo in _catalogue():
        member = thermo.resolutions(oligo.sequence)[0]
        middle = len(member) // 2
        base = member[middle]
        alternative = next(b for b in "ACGT" if b != base)
        observed = member[:middle] + alternative + member[middle + 1 :]
        cost = thermo.destabilisation(oligo.sequence, observed)
        assert cost is not None, oligo.name
        assert cost.delta_delta_g > 0, (oligo.name, cost)
        assert cost.delta_tm < 0, (oligo.name, cost)
        assert cost.observed_tm < cost.perfect_tm


def test_the_salt_correction_refuses_a_zero_ion_concentration():
    with pytest.raises(ValueError, match="monovalent"):
        thermo.salt_correction("ACGTACGT", na_mm=0.0, mg_mm=0.0, dntp_mm=0.0)


def test_magnesium_at_or_below_the_dntp_concentration_is_chelated_away():
    """dNTPs bind magnesium, so only the excess counts - and the code must say so."""
    oligo = "ACGTACGTACGT"
    chelated = thermo.salt_correction(oligo, na_mm=50.0, mg_mm=0.8, dntp_mm=0.8)
    none_at_all = thermo.salt_correction(oligo, na_mm=50.0, mg_mm=0.0, dntp_mm=0.0)
    free = thermo.salt_correction(oligo, na_mm=50.0, mg_mm=3.0, dntp_mm=0.8)
    assert chelated == none_at_all
    assert free != chelated


def test_free_energy_is_the_published_combination():
    """dG = dH - T dS, with dS in cal and dH in kcal - the unit mix is the trap."""
    assert thermo.free_energy(-10.0, 0.0, celsius=37.0) == pytest.approx(-10.0)
    assert thermo.free_energy(0.0, -1000.0, celsius=36.85) == pytest.approx(310.0)
    hotter = thermo.free_energy(-80.0, -220.0, celsius=60.0)
    cooler = thermo.free_energy(-80.0, -220.0, celsius=37.0)
    assert hotter > cooler  # a duplex is less stable when it is warmer


# -------------------------------------------------------------- the landscape's numbers


def _landscape():
    if not LANDSCAPE.exists():
        pytest.skip("research/results/thermo_landscape.json is not present")
    return json.loads(LANDSCAPE.read_text(encoding="utf-8"))


def test_the_landscape_is_exhaustive_over_the_catalogue():
    """Three alternatives at every base of every oligo, with nothing unscorable."""
    data = _landscape()
    expected = sum(len(oligo.sequence) for oligo in _catalogue()) * 3
    assert data["substitutions"] == expected
    assert data["unscorable"] == []


def test_the_positional_rule_does_not_rank_mismatches_by_binding_cost():
    """Section 4.6's headline, and the refutation of its own written prediction.

    An AUC near one half says the rule carries no information about dGG. This is pinned
    because it is the kind of number a later edit could move without anyone noticing,
    and the whole section rests on it being near a half rather than near one.
    """
    data = _landscape()
    assert data["auc_rule_vs_energy"] == pytest.approx(0.4605, abs=5e-4)
    assert abs(data["spearman_distance_vs_delta_delta_g"]) < 0.1
    assert "dGG rises towards the 3' terminus" in data["prediction_written_before_measuring"]


def test_the_terminal_base_is_the_cheap_one_in_the_published_landscape():
    data = _landscape()["delta_delta_g"]
    assert data["terminal_base_mean"] == pytest.approx(0.811, abs=5e-4)
    assert data["window_but_not_terminal_mean"] == pytest.approx(3.704, abs=5e-4)
    # The flagged mismatches are on average cheaper than the unflagged ones, which is
    # the opposite of what the rule would predict if it were an energy proxy.
    assert data["flagged_mean"] < data["unflagged_mean"]


def test_the_anneal_verdict_keeps_the_three_way_split():
    """Collapsing the range to one threshold is the error this guards."""
    data = _landscape()["observed_named_mutations"]
    if not data["available"]:
        pytest.skip("results/mutations_big.json is not present")
    verdicts = {}
    for row in data["rows"]:
        verdicts[row["anneal_verdict"]] = verdicts.get(row["anneal_verdict"], 0) + 1
    assert verdicts == {"clears": 18, "depends_on_anneal_temperature": 6, "below": 1}
    # Every row carries both ends of the range, not one "margin".
    for row in data["rows"]:
        assert row["margin_over_anneal_floor_c"] > row["margin_over_anneal_ceiling_c"]


def test_the_mismatches_that_lose_their_margin_are_mostly_invisible_to_the_rule():
    """The gap section 4.6 reports: five of seven sit outside the 3' window."""
    data = _landscape()["observed_named_mutations"]
    if not data["available"]:
        pytest.skip("results/mutations_big.json is not present")
    at_risk = [r for r in data["rows"] if r["anneal_verdict"] != "clears"]
    assert len(at_risk) == 7
    assert sum(1 for r in at_risk if not r["in_three_prime_window"]) == 5
    # And none of them is a probe, which is why the rule is right to ignore probe 3' ends.
    assert [r for r in at_risk if r["role"] == "probe"] == []


def test_no_probe_mutation_loses_its_annealing_margin():
    """The overclaim this section nearly made, kept refuted.

    A 6.6 degree Tm drop on a probe carried by 57% of the corpus reads as alarming until
    the margin is computed: hydrolysis probes are designed hot and every one of them
    still clears the range with degrees to spare.
    """
    data = _landscape()["observed_named_mutations"]
    if not data["available"]:
        pytest.skip("results/mutations_big.json is not present")
    probes = [r for r in data["rows"] if r["role"] == "probe"]
    assert len(probes) >= 10
    for row in probes:
        assert row["anneal_verdict"] == "clears", row
        assert row["margin_over_anneal_ceiling_c"] > 4.0, row


# ------------------------------------------------------- section 4.6's table vs the file


PAPER = ROOT / "paper" / "main.md"


def _section_46_table():
    """The seven rows of section 4.6's table, parsed out of the manuscript.

    A table typed into prose beside an artefact drifts from it. This repository's habit is
    to read the prose back and compare, and a markdown table is parseable enough for that.
    """
    text = PAPER.read_text(encoding="utf-8")
    start = text.index("## 4.6 A second opinion that is not a rule of thumb")
    end = text.index("\n# 5. Threats to validity", start)
    rows = {}
    for line in text[start:end].splitlines():
        if not line.startswith("| "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) != 9 or cells[0] in {"mutation", "---"} or cells[0].startswith("-"):
            continue
        rows[cells[0]] = cells
    return rows


def test_the_manuscript_table_has_the_seven_rows_the_artefact_has():
    data = _landscape()["observed_named_mutations"]
    if not data["available"]:
        pytest.skip("results/mutations_big.json is not present")
    at_risk = {r["mutation"] for r in data["rows"] if r["anneal_verdict"] != "clears"}
    assert set(_section_46_table()) == at_risk


def test_every_number_in_the_manuscript_table_matches_the_artefact():
    data = _landscape()["observed_named_mutations"]
    if not data["available"]:
        pytest.skip("results/mutations_big.json is not present")
    measured = {r["mutation"]: r for r in data["rows"]}
    table = _section_46_table()
    assert table, "section 4.6's table was not found"
    for mutation, cells in table.items():
        row = measured[mutation]
        _, oligo, from_three, in_window, ddg, dtm, tm, _verdict, share = cells
        assert oligo == row["oligo"], mutation
        assert int(from_three) == row["bases_from_3_prime_end"], mutation
        assert (in_window == "yes") == row["in_three_prime_window"], mutation
        assert float(ddg) == pytest.approx(row["delta_delta_g"], abs=5e-3), mutation
        # The manuscript writes a minus sign, not a hyphen.
        assert float(dtm.replace("−", "-")) == pytest.approx(row["delta_tm"], abs=5e-3), mutation
        assert float(tm) == pytest.approx(row["observed_tm"], abs=5e-3), mutation
        assert float(share.rstrip("%")) == pytest.approx(row["overall_rate"] * 100, abs=5e-3), (
            mutation
        )


def test_the_manuscript_does_not_claim_the_model_predicts_failure():
    """The caution the section rests on; a later edit must not quietly drop it."""
    text = PAPER.read_text(encoding="utf-8")
    start = text.index("## 4.6 A second opinion that is not a rule of thumb")
    section = " ".join(text[start : text.index("\n# 5. Threats", start)].split())
    assert "is not a failure prediction" in section
    assert "Neither is a calibration of the other" in section
    # The refuted prediction has to stay in the section that refutes it.
    assert "written before the measurement" in section


README = ROOT / "README.md"


def _readme_table():
    """Section 4.6's table as the README repeats it, keyed by mutation.

    The README's copy has eight columns where the manuscript's has nine - it carries no
    verdict column - so the two are parsed separately and compared on what they share.
    """
    text = README.read_text(encoding="utf-8")
    start = text.index("### Does the 3' rule track binding energy? No")
    end = text.index("\n## Tech stack", start)
    rows = {}
    for line in text[start:end].splitlines():
        if not line.startswith("| "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) != 8 or cells[0] == "mutation" or cells[0].startswith("-"):
            continue
        rows[cells[0]] = cells
    return rows


def test_the_readme_table_and_the_manuscript_table_say_the_same_thing():
    """Two copies of one table is two chances to drift; nothing else compares them."""
    readme = _readme_table()
    paper = _section_46_table()
    assert readme, "the README table was not found"
    assert set(readme) == set(paper)
    for mutation, cells in readme.items():
        # README: mutation, oligo, from 3', in window, ddG, dTm, Tm, share
        # paper:   mutation, oligo, from 3', in window, ddG, dTm, Tm, verdict, share
        other = paper[mutation]
        assert cells[1:7] == [c.replace("−", "-") for c in other[1:7]], mutation
        assert cells[7] == other[8], mutation


def test_the_readme_quotes_the_landscape_headline_numbers():
    data = _landscape()
    text = README.read_text(encoding="utf-8")
    start = text.index("### Does the 3' rule track binding energy? No")
    section = " ".join(text[start : text.index("\n## Tech stack", start)].split())
    assert f"**{data['auc_rule_vs_energy']}**" in section
    # Quoted to three places in prose where the artefact holds four, so the comparison is
    # made at the precision the prose uses. Asserting the artefact's own repr would fail
    # on a correct rounding, and relaxing it to a substring would pass on a wrong digit.
    rho = data["spearman_distance_vs_delta_delta_g"]
    assert f"**{rho:.3f}**" in section, (rho, f"{rho:.3f}")
    assert f"**{data['delta_delta_g']['terminal_base_mean']} kcal/mol**" in section
    assert f"**{data['delta_delta_g']['window_but_not_terminal_mean']}**" in section
