"""The size of section 4.6's unscorable case: adjacent mismatches, which no table covers.

Section 4.6 scores a mismatch with the nearest-neighbour model and names a limitation it
does not size: there are no published parameters for a **tandem** mismatch, two mismatches
at adjacent positions of one duplex. Allawi & SantaLucia's internal tables and SantaLucia &
Peyret's terminal ones each describe a mismatch whose neighbours are Watson-Crick pairs, so
putting two side by side leaves a stack that appears in no table and `thermo.duplex` returns
`None`.

`research/tandem_gap.py` counts how often that happens over a real corpus. The counting has
one way to go badly wrong and it goes wrong in the direction that inflates the answer, so
most of this file is about that: **the comparison must be the library's own IUPAC rule**.
Three of the fifteen published oligos carry a deliberate degenerate base - `CDC N3`'s probe
an `Y`, `Charité RdRp`'s primers an `R` and an `S` - and a plain `target != wanted` counts a
target `A` under a primer `R` as a mismatch when `R` means A or G. Written that way, the
scan reported 400 of 400 genomes as mismatched at `N3-P` where the real number is 1, and
every one of those invented mismatches was a candidate for an invented adjacent pair.

The denominators are the other half. A share of tandem alignments is only meaningful against
a stated population, and three are reported: all alignments, the mismatched ones (the
population the model is asked about), and the scorable ones.
"""

from __future__ import annotations

import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "research"))

tandem_gap = pytest.importorskip("tandem_gap")

MEASURED = ROOT / "research" / "results" / "tandem_gap.json"


def _positions(window: str, oligo: str):
    return tandem_gap.mismatch_positions(window, oligo)


# ------------------------------------------------- the comparison, which is the whole risk


def test_a_degenerate_oligo_base_is_not_a_mismatch():
    """`R` means A or G, so a target A under it matches. The bug this file exists for."""
    assert _positions("ACGT", "RCGT") == []
    assert _positions("GCGT", "RCGT") == []
    # And a target outside what the code allows still is one.
    assert _positions("CCGT", "RCGT") == [0]
    assert _positions("TCGT", "RCGT") == [0]


def test_every_degenerate_code_in_the_published_oligos_is_handled():
    """The three that actually appear: Y in N3-P, R and S in the RdRp primers."""
    from assaydrift.match import IUPAC

    for code, allowed in (("Y", "CT"), ("R", "AG"), ("S", "CG")):
        assert set(IUPAC[code]) == set(allowed), code
        for base in allowed:
            assert _positions(base + "CGT", code + "CGT") == [], (code, base)
        for base in set("ACGT") - set(allowed):
            assert _positions(base + "CGT", code + "CGT") == [0], (code, base)


def test_the_published_oligos_really_do_carry_degenerate_bases():
    """If they stopped, the test above would be guarding nothing."""
    from assaydrift.primers import ASSAYS

    degenerate = {
        f"{assay.name}|{oligo.name}": sorted(set(oligo.sequence.upper()) - set("ACGT"))
        for assay in ASSAYS
        for oligo in assay.oligos
        if set(oligo.sequence.upper()) - set("ACGT")
    }
    assert len(degenerate) == 3, degenerate
    assert set(code for codes in degenerate.values() for code in codes) == {"Y", "R", "S"}


def test_an_ambiguous_target_base_makes_the_alignment_unscorable():
    """An unknown base under the oligo is not evidence either way."""
    assert _positions("NCGT", "ACGT") is None
    assert _positions("RCGT", "ACGT") is None


def test_a_short_window_is_unscorable_rather_than_partly_scored():
    assert _positions("ACG", "ACGT") is None
    assert _positions("", "ACGT") is None


def test_an_exact_match_has_no_mismatches():
    assert _positions("ACGTACGT", "ACGTACGT") == []


# ------------------------------------------------------------------------ adjacency


def test_adjacent_mismatches_are_the_case_with_no_parameters():
    assert tandem_gap.has_tandem([0, 1]) is True
    assert tandem_gap.has_tandem([5, 6]) is True
    assert tandem_gap.has_tandem([0, 2]) is False
    assert tandem_gap.has_tandem([0]) is False
    assert tandem_gap.has_tandem([]) is False


def test_a_run_of_three_is_also_unscorable():
    """Three in a row has no parameters either, so a pairs-only rule would miss it."""
    assert tandem_gap.has_tandem([4, 5, 6]) is True


def test_the_positions_and_the_adjacency_test_compose():
    assert tandem_gap.has_tandem(_positions("AACGT", "TTCGT")) is True
    assert tandem_gap.has_tandem(_positions("ATCGA", "TTCGT")) is False


def test_the_model_really_cannot_score_a_tandem_mismatch():
    """The premise. If `thermo` ever gains the parameters, this file is measuring nothing."""
    from assaydrift import thermo

    assert thermo.duplex("ACGT", "TAAA") is None, "a tandem mismatch is now scorable"
    # While a single internal mismatch is scorable, which is what makes the gap specific.
    assert thermo.duplex("ACGTACGT", "TGCATGCT") is not None


# ------------------------------------------------------------ the artefact, when present


def test_the_shares_have_their_denominators_named():
    if not MEASURED.exists():
        pytest.skip("research/results/tandem_gap.json is not present")
    payload = json.loads(MEASURED.read_text(encoding="utf-8"))
    assert set(payload["shares"]) == set(payload["denominators"]), (
        "a share without a stated denominator is not a measurement"
    )
    for name, text in payload["denominators"].items():
        assert len(text) > 20, name


def test_a_finished_scan_is_not_marked_partial():
    if not MEASURED.exists():
        pytest.skip("research/results/tandem_gap.json is not present")
    payload = json.loads(MEASURED.read_text(encoding="utf-8"))
    assert payload["partial"] is False, (
        "this file is a checkpoint from a run that did not finish; its shares are over the "
        "genomes read so far, which is a population nobody chose"
    )


def test_the_counts_add_up():
    if not MEASURED.exists():
        pytest.skip("research/results/tandem_gap.json is not present")
    totals = json.loads(MEASURED.read_text(encoding="utf-8"))["totals"]
    assert (
        totals["indel"] + totals["unknown_or_short"] + totals["scorable"] == (totals["alignments"])
    )
    # A tandem alignment is a mismatched one, so it cannot exceed that count.
    assert totals["tandem"] <= totals["with_a_mismatch"] <= totals["scorable"]


def test_every_oligo_row_sums_into_the_totals():
    if not MEASURED.exists():
        pytest.skip("research/results/tandem_gap.json is not present")
    payload = json.loads(MEASURED.read_text(encoding="utf-8"))
    rows = payload["per_oligo"].values()
    for key in ("alignments", "with_a_mismatch", "scorable", "tandem", "indel"):
        assert sum(row[key] for row in rows) == payload["totals"][key], key
    assert len(payload["per_oligo"]) == payload["oligos"] == 15
