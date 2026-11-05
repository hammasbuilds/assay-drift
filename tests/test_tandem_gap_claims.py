"""Section 4.6's tandem-mismatch figures, against `research/results/tandem_gap.json`.

`assaydrift.thermo` refuses to score adjacent mismatches, because no published
nearest-neighbour table covers that stack, and the manuscript argued the refusal mattered
from a single example - C28311T and C28312T on CDC N1's probe. That establishes the case is
not hypothetical and says nothing about its size, so the size is now measured over all
291,827 genomes.

Every figure below is recomputed from the per-oligo rows rather than read out of the
artefact's own totals, so a bug in the summing cannot make the summary agree with itself.
The one exception is `worst_run_anywhere`, which is not derivable from the per-oligo rows
and is checked against them for consistency instead.
"""

from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "main.md"
RESULT = ROOT / "research" / "results" / "tandem_gap.json"


def _result() -> dict:
    if not RESULT.is_file():
        pytest.skip("run research/tandem_gap.py --data data/sequences_big.jsonl")
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    if payload.get("partial"):
        pytest.skip("the scan is still running; its artefact is partial")
    return payload


def _rows() -> list[dict]:
    return list(_result()["per_oligo"].values())


def _paper() -> str:
    """The manuscript with runs of whitespace flattened to single spaces.

    Every check below is for a phrase, and the manuscript hard-wraps its paragraphs, so a
    phrase that spans a line break would not be found in the raw text. Flattening is the
    difference between testing the prose and testing where the author pressed return.
    """
    if not PAPER.is_file():
        pytest.skip("the manuscript is not shipped with the tests")
    return re.sub(r"\s+", " ", PAPER.read_text(encoding="utf-8"))


def test_the_scan_covered_the_whole_large_corpus():
    payload = _result()
    assert payload["genomes"] == 291827, payload["genomes"]
    assert payload["oligos"] == 15, payload["oligos"]
    assert "sequences_big" in payload["corpus"], payload["corpus"]


def test_the_per_oligo_rows_sum_to_the_published_totals():
    """The summing this file then relies on."""
    rows, totals = _rows(), _result()["totals"]
    for key in ("alignments", "with_a_mismatch", "scorable", "tandem", "indel"):
        assert sum(r[key] for r in rows) == totals[key], key


def test_the_headline_rate_is_the_measured_one():
    """3,481 of 736,212 mismatched alignments, 0.47%."""
    rows = _rows()
    tandem = sum(r["tandem"] for r in rows)
    mismatched = sum(r["with_a_mismatch"] for r in rows)
    assert tandem == 3481, tandem
    assert mismatched == 736212, mismatched
    assert f"{tandem / mismatched * 100:.2f}" == "0.47"
    paper = _paper()
    assert f"{tandem:,} of {mismatched:,}" in paper
    assert "0.47%" in paper


def test_indels_are_counted_and_kept_out_of_the_tandem_count():
    """A gap is a different kind of unscorable from a stack with no parameters."""
    rows = _rows()
    indels = sum(r["indel"] for r in rows)
    assert indels == 794, indels
    assert f"({indels:,})" in _paper() or f"excluded ({indels})" in _paper()
    # And it is not folded into the tandem figure.
    assert indels != sum(r["tandem"] for r in rows)


def test_the_case_is_concentrated_on_one_probe():
    """97.6% of it, which is why the paper's single example came from that oligo."""
    rows = _rows()
    total = sum(r["tandem"] for r in rows)
    worst = max(rows, key=lambda r: r["tandem"])
    assert (worst["assay"], worst["oligo"]) == ("CDC N1", "N1-P"), worst
    assert worst["tandem"] == 3398, worst["tandem"]
    assert f"{worst['tandem'] / total * 100:.1f}" == "97.6"
    assert f"{worst['tandem_share_of_mismatched'] * 100:.2f}" == "2.03"
    paper = _paper()
    assert f"{worst['tandem']:,} of the {total:,}" in paper
    assert "97.6%" in paper and "2.03%" in paper


def test_no_other_oligo_comes_close():
    """The claim that the concentration is real and not an artefact of one big denominator."""
    rows = _rows()
    others = [r for r in rows if r["oligo"] != "N1-P"]
    assert max(r["tandem_share_of_mismatched"] for r in others) < 0.005
    assert "no other oligo exceeds 0.5%" in _paper()


def test_only_eight_of_the_fifteen_oligos_hit_the_case_at_all():
    rows = _rows()
    assert sum(1 for r in rows if r["tandem"]) == 8
    assert "Eight of the fifteen" in _paper()


def test_the_primer_figures_the_section_leans_on():
    """The 3' rule is about primers, so this is the number that bounds §4.6's argument."""
    rows = _rows()
    primers = [r for r in rows if r["role"] != "probe"]
    probes = [r for r in rows if r["role"] == "probe"]
    assert len(primers) == 10 and len(probes) == 5, (len(primers), len(probes))
    tandem = sum(r["tandem"] for r in primers)
    mismatched = sum(r["with_a_mismatch"] for r in primers)
    assert tandem == 77, tandem
    assert mismatched == 558723, mismatched
    assert f"{tandem / mismatched * 100:.4f}" == "0.0138"
    paper = _paper()
    assert f"**{tandem} times in {mismatched:,} mismatched" in paper
    assert "0.0138%" in paper
    assert "all ten primers" in paper, "the paper miscounts the primers"
    assert (
        f"{sum(r['tandem'] for r in probes):,} in {sum(r['with_a_mismatch'] for r in probes):,}"
        in paper
    )
    # "about one alignment in seven thousand"
    assert 7000 <= round(mismatched / tandem) < 8000, round(mismatched / tandem)


def test_the_longest_adjacent_run_is_consistent_with_the_per_oligo_rows():
    """Not derivable from them, so checked for consistency rather than recomputed."""
    payload = _result()
    rows = _rows()
    assert payload["worst_run_anywhere"] == max(r["worst_run"] for r in rows)
    assert payload["worst_run_anywhere"] == 3
    assert "longest run of adjacent mismatches anywhere in the corpus is **three**" in _paper()


def test_an_oligo_with_no_tandem_case_still_reports_its_denominator():
    """Otherwise a zero could mean 'none found' or 'never looked'."""
    rows = _rows()
    quiet = [r for r in rows if r["tandem"] == 0]
    assert quiet, "every oligo hit the case, which the concentration claim denies"
    for row in quiet:
        assert row["alignments"] > 0, row
        assert row["scorable"] > 0, row
