"""Every number the README quotes must match the committed results.

Written because a claim that these matched was made on the strength of a comparison that
read `overall_rate` as `overall`, got `None` on both sides for all 28 rows, and reported
"0 moved". The rates had in fact shifted in 17 of 28 rows and one README figure was stale
by 0.2 points. A check that cannot fail is worse than no check, because it is quoted.

These tests read `results/*.json` - the committed output of a real run over the 2,765
genomes in the manifest - and never the network.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
README = ROOT / "README.md"
DRIFT = ROOT / "results" / "drift.json"
MUTATIONS = ROOT / "results" / "mutations.json"


def _mutation_rates() -> dict[str, float]:
    data = json.loads(MUTATIONS.read_text(encoding="utf-8"))
    out: dict[str, float] = {}
    for rows in data.values():
        if not isinstance(rows, list):
            continue
        for row in rows:
            name, rate = row.get("name"), row.get("overall_rate")
            if name and rate is not None:
                out[name] = rate * 100
    return out


def test_the_mutation_table_has_the_rows_we_think_it_has() -> None:
    """Guards the key name. The vacuous comparison went unnoticed because it read a key
    that does not exist, so nothing here may assume a key without asserting it."""
    data = json.loads(MUTATIONS.read_text(encoding="utf-8"))
    rows = [r for v in data.values() if isinstance(v, list) for r in v]
    assert rows, "no mutation rows at all"
    for row in rows:
        assert "overall_rate" in row, f"no overall_rate in {sorted(row)}"
        assert "name" in row and "oligo" in row
    assert _mutation_rates(), "rates could not be read"


def test_every_mutation_percentage_in_the_readme_matches_the_results() -> None:
    rates = _mutation_rates()
    text = README.read_text(encoding="utf-8")
    checked = 0
    for match in re.finditer(r"([ACGT]\d{4,5}[ACGT])[^\n]*?(\d{1,2}\.\d)%", text):
        name, claimed = match.group(1), float(match.group(2))
        if name not in rates:
            continue
        checked += 1
        actual = rates[name]
        assert abs(actual - claimed) < 0.05, (
            f"{name}: README says {claimed}%, results/mutations.json says {actual:.1f}%"
        )
    assert checked >= 3, f"only {checked} mutation figures found - has the README changed?"


def test_every_assay_rate_in_the_readme_table_matches_the_results() -> None:
    """The quarter table's first, latest and change columns."""
    drift = json.loads(DRIFT.read_text(encoding="utf-8"))
    # Only the QUARTER table. The by-year section below it regroups the same sequences and
    # its figures legitimately differ, so checking both against `trend` would compare a
    # year rate to a quarter rate and fail for the right reason in the wrong place.
    whole = README.read_text(encoding="utf-8")
    assert "### The same data by year" in whole, "the by-year heading moved"
    text = whole.split("### The same data by year", 1)[0]
    names = {
        "CDC N1": "CDC N1",
        "Charité E": "Charite E",
        "Charité RdRp": "Charite RdRp",
        "CDC N3": "CDC N3",
        "CDC N2": "CDC N2",
    }
    checked = 0
    for line in text.splitlines():
        if not line.startswith("| **") and not line.startswith("| CDC N3"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        label = cells[0].replace("*", "").split("(")[0].strip()
        key = names.get(label)
        if key is None or key not in drift:
            continue
        trend = drift[key].get("trend") or {}
        if not trend.get("enough_to_say"):
            continue
        for index, field in ((1, "first_rate"), (2, "last_rate")):
            claimed = float(cells[index].replace("*", "").replace("%", ""))
            actual = trend[field] * 100
            assert abs(actual - claimed) < 0.05, (
                f"{label} {field}: README {claimed}%, results {actual:.1f}%"
            )
            checked += 1
    assert checked >= 8, f"only {checked} table cells checked"


def test_the_results_record_which_corpus_produced_them() -> None:
    """A rate with no provenance cannot be re-derived, and this corpus is gitignored."""
    drift = json.loads(DRIFT.read_text(encoding="utf-8"))
    run = drift.get("_run")
    assert run, "results/drift.json has no _run provenance block"
    assert run.get("corpus_identity_sha256"), "no corpus identity recorded"
    assert run.get("manifest") == "data/accessions.tsv"
    assert run.get("sequences_read") == run.get("manifest_records"), (
        "the run read a different number of sequences than the manifest lists"
    )
    manifest = ROOT / "data" / "accessions.tsv"
    assert manifest.is_file(), "the manifest the provenance points at is not committed"
    listed = [
        line
        for line in manifest.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#") and not line.startswith("accession\t")
    ]
    assert len(listed) == run["manifest_records"]


def test_every_period_rate_reports_how_independent_it_is() -> None:
    """`usable` counts sequences; sequences from one submission are not independent.

    2024-Q1 holds 159 sequences from 4 submitter groups - an effective n near 1 - and
    passes the 25-sequence floor. A rate published without that figure invites being read
    as 159 observations.
    """
    drift = json.loads(DRIFT.read_text(encoding="utf-8"))
    worst = None
    for assay, entry in drift.items():
        if assay.startswith("_"):
            continue
        for period, row in (entry.get("periods") or {}).items():
            assert "clusters" in row, f"{assay} {period} has no cluster count"
            assert "effective_n" in row, f"{assay} {period} has no effective_n"
            fine = row.get("reliable") and row.get("effective_n") is not None
            if fine and (worst is None or row["effective_n"] < worst[2]):
                worst = (assay, period, row["effective_n"])
    assert worst is not None
    # Not an assertion about the data being good - an assertion that the weakest case is
    # still visible, because this is the number the headline rests on.
    assert worst[2] < 2.0, (
        "the corpus no longer contains a reliable period with a tiny effective n - if "
        "that is real, the README's independence caveat needs rewriting"
    )


def test_the_readme_independence_table_matches_the_results() -> None:
    """The caveat's own numbers are numbers, so they are checked like the others."""
    drift = json.loads(DRIFT.read_text(encoding="utf-8"))
    periods = drift["CDC N2"]["periods"]
    text = README.read_text(encoding="utf-8")
    checked = 0
    for line in text.splitlines():
        match = re.match(
            r"\|\s*(\d{4}-Q\d)[^|]*\|\s*(\d+)\s*\|\s*\*{0,2}(\d+)\*{0,2}\s*\|"
            r"\s*\*{0,2}([\d.]+)\*{0,2}\s*\|",
            line,
        )
        if not match:
            continue
        period, sequences, clusters, effective = match.groups()
        if period not in periods:
            continue
        row = periods[period]
        assert int(sequences) == row["total"], f"{period} sequences"
        assert int(clusters) == row["clusters"], f"{period} clusters"
        # Compared at the precision the README prints, not with a tolerance: 1.05 displays
        # as 1.1 and |1.1 - 1.05| is exactly 0.05, which no sane tolerance straddles.
        assert float(effective) == round(row["effective_n"], 1), (
            f"{period} effective_n: README {effective}, results {row['effective_n']}"
        )
        checked += 1
    assert checked >= 4, f"only {checked} independence rows checked"


def test_demo_runs_against_the_committed_results() -> None:
    """Catches the class, not just the instance.

    Adding the `_run` provenance block to drift.json broke demo.py, which iterated the
    top level and found a key with no `trend`. The README's whole convention is that
    demo.py reprints the published numbers, so a consumer that dies on a new metadata key
    is a broken promise - and nothing was checking it.
    """
    import os
    import subprocess
    import sys

    # An empty PATH is the point: demo.py must not reach for anything on it. But a
    # wholly stripped environment does not start CPython on Windows at all - it needs
    # SystemRoot to seed hash randomisation, and without it 3.10 dies with
    # "_Py_HashRandomization_Init: failed to get random numbers" before running a line.
    # So the interpreter's own variables are passed through and nothing else.
    env = {"PYTHONIOENCODING": "utf-8", "PATH": ""}
    for needed in ("SystemRoot", "SYSTEMROOT", "windir", "TEMP", "TMP"):
        if needed in os.environ:
            env[needed] = os.environ[needed]

    done = subprocess.run(
        [sys.executable, str(ROOT / "demo.py")],
        capture_output=True,
        cwd=ROOT,
        timeout=600,
        check=False,
        env=env,
    )
    output = done.stdout.decode("utf-8", "replace")
    assert done.returncode == 0, done.stderr.decode("utf-8", "replace")[-800:]
    # It must print the finding AND say which corpus produced it.
    assert "exact match" in output
    assert "Corpus identity:" in output, "demo.py does not state the corpus provenance"
    for assay in ("CDC N1", "CDC N2", "Charite RdRp"):
        assert assay in output, f"{assay} missing from demo output"


def test_the_c29215t_quarter_table_matches_the_results():
    """The per-quarter carriage table, which no test covered.

    Its 2026-Q1 row read 30.6% over n=108 while a paragraph twenty lines above said
    28.7% over n=115 for the same quarter - the README disagreeing with itself, because
    one figure was updated after a re-run and the other was not. Both now come from
    `results/`, so neither can move without this failing.
    """
    drift = json.loads(DRIFT.read_text(encoding="utf-8"))
    mutations = json.loads(MUTATIONS.read_text(encoding="utf-8"))
    row = next(r for r in mutations["CDC N2"] if (r.get("mutation") or r.get("name")) == "C29215T")
    periods = drift["CDC N2"]["periods"]
    readme = README.read_text(encoding="utf-8")
    for quarter, carrying in row["by_period"].items():
        usable = periods[quarter]["usable"]
        rate = 100 * carrying / usable
        line = f"| {quarter} | {rate:.1f}% | {usable} |"
        assert line in readme, (
            f"the README's {quarter} row is not {line!r} - results say {carrying} of "
            f"{usable} readable sequences carry C29215T"
        )


def test_the_pooled_2026_carriage_figure_matches_the_results():
    """ "51 of 333 readable sequences (15.3%)" has to be those two numbers.

    An earlier version of this sentence said 53/333 and 15.9%, and also quoted a
    Wisconsin-versus-UK split for 2026 alone that nothing in `results/` emits - so
    there was no way to check it and it was simply carried forward. The figures quoted
    now are ones the pipeline writes down.
    """
    drift = json.loads(DRIFT.read_text(encoding="utf-8"))
    mutations = json.loads(MUTATIONS.read_text(encoding="utf-8"))
    row = next(r for r in mutations["CDC N2"] if (r.get("mutation") or r.get("name")) == "C29215T")
    periods = drift["CDC N2"]["periods"]
    carrying = sum(n for q, n in row["by_period"].items() if q.startswith("2026"))
    usable = sum(p["usable"] for q, p in periods.items() if q.startswith("2026"))
    readme = README.read_text(encoding="utf-8")
    assert f"{carrying} of {usable} readable" in readme, (
        f"the README does not say {carrying} of {usable} readable sequences"
    )
    assert f"({100 * carrying / usable:.1f}%)" in readme


def test_the_top_places_quoted_are_the_ones_the_pipeline_stored():
    """The clustering argument rests on this list, so it has to be the stored one.

    It quoted `USA: Missouri 4/26`, which the pipeline no longer reports at all - the
    place list changed when the corpus was corrected and the prose did not.
    """
    mutations = json.loads(MUTATIONS.read_text(encoding="utf-8"))
    row = next(r for r in mutations["CDC N2"] if (r.get("mutation") or r.get("name")) == "C29215T")
    readme = README.read_text(encoding="utf-8")
    for place in row["top_places"]:
        quoted = f"{place['place']} {place['carrying']}/{place['sequenced']}"
        assert quoted in readme, f"the README does not quote {quoted!r}"
