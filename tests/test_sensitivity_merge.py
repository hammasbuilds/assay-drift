"""The per-window split of the sensitivity grid, and what the merge has to refuse.

One pass over the 291,827-genome corpus is about three hours. Only the 3' window costs a
pass - the count is produced during alignment and `Hit` keeps the count rather than the
positions, so a window cannot be varied after the fact, while the load and the shape are
re-scored from the aligned rows in memory. So the grid is three passes, and `--window` runs
one of them while `--merge` joins the parts.

The headline the grid carries is a **range** across its cells, and that is what makes the
merge's refusals load-bearing rather than tidy:

  * a set of parts covering two of three windows still produces a range, narrower than the
    truth and indistinguishable from it in the output;
  * two parts measured over different corpora produce rates that are comparable in shape
    and about different populations;
  * two parts disagreeing about a shared cell mean one of them is stale, and silently
    keeping either would publish a number from a run nobody can identify.

Each is asserted here. The merge is also checked to reproduce the committed single-process
grid cell for cell, which is the property that says the split changed the schedule and not
the measurement.
"""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research"))

sensitivity = pytest.importorskip("sensitivity")

COMMITTED = ROOT / "research" / "results" / "sensitivity.json"


def _part(tmp_path: pathlib.Path, name: str, windows, grid: dict, corpus="sequences.jsonl"):
    path = tmp_path / name
    path.write_text(
        json.dumps(sensitivity._payload(corpus, windows, grid), indent=1), encoding="utf-8"
    )
    return path


def _cells(window: int) -> dict:
    """One cell per load and shape, as a real part would carry for that window."""
    return {
        f"Assay|w{window}|l{load}|{shape}": {"rate": 0.1, "n": 100, "failing": 10}
        for load in sensitivity.LOADS
        for shape in sensitivity.SHAPES
    }


def test_the_window_choices_are_the_grid_s_windows():
    """`--window` is restricted to them, so a part can never cover a window the grid lacks."""
    assert sensitivity.WINDOWS == (3, 5, 7)
    assert 5 in sensitivity.WINDOWS, "the published cell is at window 5"


def test_a_merge_missing_a_window_is_refused(tmp_path):
    parts = [
        _part(tmp_path, "w3.json", (3,), _cells(3)),
        _part(tmp_path, "w5.json", (5,), _cells(5)),
    ]
    with pytest.raises(SystemExit) as caught:
        sensitivity.merge(parts, tmp_path / "out.json")
    assert "cover windows" in str(caught.value)
    assert not (tmp_path / "out.json").exists(), "a refused merge must write nothing"


def test_parts_from_different_corpora_are_refused(tmp_path):
    parts = [
        _part(tmp_path, "w3.json", (3,), _cells(3), corpus="sequences.jsonl"),
        _part(tmp_path, "w5.json", (5,), _cells(5), corpus="sequences_big.jsonl"),
        _part(tmp_path, "w7.json", (7,), _cells(7), corpus="sequences.jsonl"),
    ]
    with pytest.raises(SystemExit) as caught:
        sensitivity.merge(parts, tmp_path / "out.json")
    assert "corpus" in str(caught.value)


def test_parts_disagreeing_about_a_shared_cell_are_refused(tmp_path):
    shared = _cells(5)
    stale = dict(shared)
    key = next(iter(stale))
    stale[key] = {"rate": 0.9, "n": 100, "failing": 90}
    parts = [
        _part(tmp_path, "w3.json", (3,), {**_cells(3), **shared}),
        _part(tmp_path, "w5.json", (5,), stale),
        _part(tmp_path, "w7.json", (7,), _cells(7)),
    ]
    with pytest.raises(SystemExit) as caught:
        sensitivity.merge(parts, tmp_path / "out.json")
    assert "disagree about cell" in str(caught.value)


def test_a_missing_part_is_refused_rather_than_skipped(tmp_path):
    parts = [
        _part(tmp_path, "w3.json", (3,), _cells(3)),
        tmp_path / "absent.json",
        _part(tmp_path, "w7.json", (7,), _cells(7)),
    ]
    with pytest.raises(SystemExit) as caught:
        sensitivity.merge(parts, tmp_path / "out.json")
    assert "missing part" in str(caught.value)


def test_a_complete_merge_keeps_every_cell_and_records_its_parts(tmp_path):
    parts = [
        _part(tmp_path, "w3.json", (3,), _cells(3)),
        _part(tmp_path, "w5.json", (5,), _cells(5)),
        _part(tmp_path, "w7.json", (7,), _cells(7)),
    ]
    out = tmp_path / "out.json"
    assert sensitivity.merge(parts, out) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["windows"] == list(sensitivity.WINDOWS)
    assert payload["merged_from"] == ["w3.json", "w5.json", "w7.json"]
    expected = len(sensitivity.WINDOWS) * len(sensitivity.LOADS) * len(sensitivity.SHAPES)
    assert len(payload["grid"]) == expected


def test_the_split_reproduces_the_committed_grid(tmp_path):
    """The property that matters: a different schedule, not a different measurement.

    Run as three subprocesses over the small corpus, because `features` sets a module
    global before the first alignment and a second value of it inside one process would be
    measuring the wrong window. That is also why the windows are separate processes in the
    real run rather than a loop.
    """
    if not COMMITTED.exists():
        pytest.skip("research/results/sensitivity.json is not present")
    data = ROOT / "data" / "sequences.jsonl"
    if not data.exists():
        pytest.skip("the small corpus is not on disk")
    python = ROOT / ".venv-research" / "Scripts" / "python.exe"
    if not python.exists():
        pytest.skip("the research environment is not present")
    parts = []
    for window in sensitivity.WINDOWS:
        part = tmp_path / f"w{window}.json"
        done = subprocess.run(  # noqa: S603
            [
                str(python),
                str(ROOT / "research" / "sensitivity.py"),
                "--data",
                str(data),
                "--window",
                str(window),
                "--out",
                str(part),
            ],
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
        )
        assert done.returncode == 0, done.stderr[-2000:]
        parts.append(part)
    out = tmp_path / "merged.json"
    assert sensitivity.merge(parts, out) == 0
    merged = json.loads(out.read_text(encoding="utf-8"))["grid"]
    committed = json.loads(COMMITTED.read_text(encoding="utf-8"))["grid"]
    assert set(merged) == set(committed)
    differing = {k: (committed[k], merged[k]) for k in committed if committed[k] != merged[k]}
    assert not differing, differing
