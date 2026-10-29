"""docs/COMPARISON.md agrees with docs/aligner-comparison.json.

The comparison against Biopython is the only document here that quotes another tool's
output, and the one with the most room to flatter this one. So every figure in it is
checked against the artefact the comparison script wrote, and the two claims that carry
the page - that the full-length disagreements are explained by IUPAC degeneracy, and that
the trimmed ones are a local aligner dropping an end - are checked against the recorded
rows rather than against the prose.

What is not checked here: whether Biopython's alignments are right. Its column is its own
output, and the document says so.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "COMPARISON.md"
ARTEFACT = ROOT / "docs" / "aligner-comparison.json"

pytestmark = pytest.mark.skipif(
    not (DOC.exists() and ARTEFACT.exists()),
    reason="docs/ is not shipped with the tests",
)


@pytest.fixture(scope="module")
def data() -> dict:
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def doc() -> str:
    return DOC.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def full_length(data) -> list[dict]:
    """Disagreements where Biopython aligned every base of the oligo."""
    return [r for r in data["disagreements"] if not (r.get("biopython_unaligned") or 0)]


@pytest.fixture(scope="module")
def trimmed(data) -> list[dict]:
    return [r for r in data["disagreements"] if (r.get("biopython_unaligned") or 0) > 0]


def test_the_run_the_document_describes_is_the_run_in_the_artefact(data, doc):
    assert f"Biopython {data['biopython_version']}" in doc
    assert f"seed {data['seed']}" in doc
    assert f"**{data['oligos']} oligos**" in doc
    assert f"**{data['genomes']} genomes**" in doc
    assert f"**{data['pairs']} pairs**" in doc


def test_the_result_table_is_the_artefacts(data, doc):
    for label, key in (
        ("pairs", "pairs"),
        ("both put the oligo at the same base", "both_found_the_same_site"),
        ("…and counted the same mismatches", "and_agreed_on_the_mismatch_count"),
        ("disagreed", "disagreed"),
        ("sites found only by Biopython", "found_only_by_biopython"),
        ("Biopython left part of the oligo unaligned", "biopython_left_oligo_bases_unaligned"),
    ):
        value = data[key]
        assert re.search(rf"\| {re.escape(label)} \| \*{{0,2}}{value}\*{{0,2}} \|", doc), (
            f"the document's {label!r} row does not read {value}"
        )


def test_the_seconds_are_rounded_from_the_artefact(data, doc):
    """Both tools' time, as whole seconds, because they move run to run."""
    ours, theirs = round(data["assay_drift_seconds"]), round(data["biopython_seconds"])
    assert f"**{ours}** / {theirs}" in doc
    assert data["assay_drift_seconds"] < data["biopython_seconds"], (
        "the document presents this tool as the faster one"
    )


def test_biopython_found_no_site_this_tool_missed(data, doc):
    assert data["found_only_by_biopython"] == 0
    assert not data["found_only_by_biopython_rows"]
    assert "Biopython never found a binding site this tool missed" in doc


def test_the_disagreements_fall_on_the_number_of_oligos_claimed(data, doc):
    oligos = {r["oligo"] for r in data["disagreements"]}
    assert f"**{len(oligos)} of\nthe {data['oligos']} oligos**" in doc or \
           f"**{len(oligos)} of the {data['oligos']} oligos**" in doc
    assert len(data["disagreements"]) == data["disagreed"]


# --- the degeneracy half -------------------------------------------------------------


def test_every_full_length_disagreement_is_biopython_counting_higher(full_length, doc):
    """The direction is the claim: higher, never lower."""
    assert full_length
    for row in full_length:
        assert (row.get("biopython_mismatches") or 0) > (row.get("mismatches") or 0), (
            f"{row['oligo']} on {row['accession']} has Biopython no higher"
        )
    assert f"In **all {len(full_length)}** of these Biopython is *higher*, never lower" in doc


def test_the_gap_never_exceeds_the_oligos_degenerate_base_count(full_length, doc):
    for row in full_length:
        degenerate = row.get("degenerate_bases_in_oligo") or []
        gap = (row["biopython_mismatches"] or 0) - (row["mismatches"] or 0)
        assert 0 < gap <= len(degenerate), (
            f"{row['oligo']}: Biopython is {gap} higher with {len(degenerate)} "
            "degenerate base(s)"
        )
    assert "never exceeds\nthat oligo's degenerate-base count" in doc or \
           "never exceeds that oligo's degenerate-base count" in doc


def test_plain_letter_equality_accounts_for_every_one_of_them(data, full_length, doc):
    """The cause, proved rather than asserted.

    At the site this tool reports, the oligo compared to the window with plain equality -
    no IUPAC - must give Biopython's number. That is the one variable between them.
    """
    explained = [
        r for r in full_length
        if r.get("mismatches_plain_identity") == r.get("biopython_mismatches")
        and r.get("mismatches_plain_identity") != r.get("mismatches")
    ]
    assert len(explained) == len(full_length), (
        f"only {len(explained)} of {len(full_length)} full-length disagreements are "
        "explained by plain letter equality"
    )
    assert data["iupac_explains"] == len(full_length)
    assert f"**All {len(full_length)}** of Biopython's counts\nequal that plain count" in doc \
        or f"**All {len(full_length)}** of Biopython's counts equal that plain count" in doc


def test_the_degenerate_oligos_in_the_table_are_the_ones_with_degenerate_bases(doc):
    """Read from the catalogue, not from the artefact, so the sequences are the real ones."""
    sys.path.insert(0, str(ROOT / "src"))
    from assaydrift.primers import ASSAYS  # noqa: PLC0415

    degenerate = {
        o.name: sorted(set(o.sequence) - set("ACGT"))
        for assay in ASSAYS
        for o in assay.oligos
    }
    named = {name for name, codes in degenerate.items() if codes}
    assert named == {"RdRp-F", "RdRp-R", "N3-P"}, named
    for name in named:
        sequence = next(
            o.sequence for assay in ASSAYS for o in assay.oligos if o.name == name
        )
        assert f"`{sequence}`" in doc, f"{name}'s sequence is not the one in the document"
    assert "Three of the fifteen oligos" in doc


# --- the trimming half ---------------------------------------------------------------


def test_every_trimmed_disagreement_lost_exactly_two_bases(trimmed, doc):
    assert trimmed
    assert {r["biopython_unaligned"] for r in trimmed} == {2}
    assert "always by exactly two bases" in doc


def test_the_trimming_table_names_the_right_direction(trimmed, doc):
    """`E-F` understates and `N3-P` does not, which is the point of that table."""
    by_oligo: dict[str, list[dict]] = {}
    for row in trimmed:
        by_oligo.setdefault(row["oligo"], []).append(row)
    assert set(by_oligo) == {"E-F", "N3-P"}, set(by_oligo)

    ef = by_oligo["E-F"]
    assert {r["mismatches"] for r in ef} == {1}
    assert {r["biopython_mismatches"] for r in ef} == {0}
    assert f"| **understates** | {len(ef)} |" in doc

    n3 = by_oligo["N3-P"]
    assert {r["mismatches"] for r in n3} == {0}
    assert {r["biopython_mismatches"] for r in n3} == {0}
    assert f"| equal | {len(n3)} |" in doc


def test_the_understating_oligo_has_no_degenerate_base(trimmed, doc):
    """Which is why it is the one worth reading: a real mismatch, trimmed away."""
    ef = [r for r in trimmed if r["oligo"] == "E-F"]
    assert all(not (r.get("degenerate_bases_in_oligo") or []) for r in ef)
    assert "carries no degenerate base" in doc


# --- the honesty sections -------------------------------------------------------------


def test_the_corrections_section_keeps_all_three(doc):
    section = doc.split("## Three corrections to this comparison", 1)
    assert len(section) == 2, "the corrections section is gone"
    for phrase in ("reverse-complements a reverse primer", "19 mismatches of 24",
                   "Unterminated string"):
        assert phrase in section[1], f"a correction has been dropped: {phrase}"


def test_the_limits_section_still_names_what_it_does_not_show(doc):
    limits = doc.split("## What this does not show", 1)
    assert len(limits) == 2, "the limits section is gone"
    for phrase in ("not a drift tool", "chosen by me", "not a corpus",
                   "is a recall comparison"):
        assert phrase in limits[1], f"the limits no longer say: {phrase}"


def test_the_document_does_not_claim_a_site_biopython_alone_found(data, doc):
    """A guard against the easiest way for this page to become flattering."""
    assert "sites found only by Biopython | **0**" in doc
    assert data["biopython_errors"] == 0, (
        "Biopython errored on some pairs; the document reports none"
    )
