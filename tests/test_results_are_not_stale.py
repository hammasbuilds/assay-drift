"""The committed results reproduce from the code that is in the repository now.

`tests/test_published_numbers.py` checks that the README agrees with
`results/drift.json`. That is necessary and it is not sufficient, and the gap between the
two is a defect this repository actually shipped.

On 2026-10-28, commit 277cc50 changed how trend endpoints are chosen - from a sequence
count alone to a Kish effective sample size. `results/drift.json` was last generated on
2026-10-07 and `results/drift_big.json` about three hours before that commit landed.
Neither was regenerated. The README matched the stale files, so the existing test passed,
while the current code produced different endpoints from the same corpus: CDC N1's trend
ran 2020-Q1 to 2026-Q3 in the committed file and 2022-Q4 to 2025-Q4 when recomputed, on a
corpus whose identity hash was byte-identical.

A README-versus-artefact test cannot catch that by construction, because both sides of
its comparison are stale together. This test compares the artefact against a **fresh
run**, which is the only comparison that notices.

It runs on the 2,765-genome corpus, which takes seconds. The 291,827-genome corpus takes
hours, so `drift_big.json` is checked for the weaker property that its run metadata is
not older than the last change to the code that produces it - enough to fail loudly when
it goes stale, without a three-hour test.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SMALL = ROOT / "data" / "sequences.jsonl"
DRIFT = ROOT / "results" / "drift.json"
DRIFT_BIG = ROOT / "results" / "drift_big.json"
ANALYSE = ROOT / "scripts" / "analyse.py"

pytestmark = pytest.mark.skipif(
    not (SMALL.exists() and DRIFT.exists() and ANALYSE.exists()),
    reason="the corpus and results are not shipped with the tests",
)


@pytest.fixture(scope="module")
def fresh(tmp_path_factory) -> dict:
    """Run the published command on the published corpus, into a scratch directory."""
    out = tmp_path_factory.mktemp("fresh")
    result = subprocess.run(
        [
            sys.executable,
            str(ANALYSE),
            "--data",
            str(SMALL),
            "--json",
            str(out / "drift.json"),
            "--mutations-json",
            str(out / "mutations.json"),
        ],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
        timeout=900,
    )
    assert result.returncode == 0, (
        f"analyse.py failed:\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}"
    )
    path = out / "drift.json"
    assert path.is_file(), (
        f"analyse.py exited 0 but wrote no {path}. stdout:\n{result.stdout[-2000:]}"
    )
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def committed() -> dict:
    return json.loads(DRIFT.read_text(encoding="utf-8"))


def test_the_same_corpus_was_measured(fresh, committed):
    """If the corpus differs, every other assertion here is meaningless."""
    assert fresh["_run"]["corpus_identity_sha256"] == committed["_run"]["corpus_identity_sha256"]
    assert fresh["_run"]["sequences_read"] == committed["_run"]["sequences_read"]


def assays(payload: dict) -> list[str]:
    return sorted(k for k in payload if not k.startswith("_"))


def test_the_same_assays_are_reported(fresh, committed):
    assert assays(fresh) == assays(committed)


@pytest.mark.parametrize("field", ["first_period", "last_period", "enough_to_say"])
def test_trend_endpoints_reproduce(fresh, committed, field):
    """The field that actually went stale, checked by name.

    `first_period` is the one that moved from 2020-Q1 to 2022-Q4 and went unnoticed for
    three days across two artefacts.
    """
    for assay in assays(committed):
        got = (fresh[assay].get("trend") or {}).get(field)
        want = (committed[assay].get("trend") or {}).get(field)
        assert got == want, (
            f"{assay}.trend.{field}: results/drift.json says {want!r}, a fresh run of "
            f"the current code says {got!r}. The committed artefact is stale - "
            f"regenerate it with `python scripts/analyse.py`."
        )


@pytest.mark.parametrize("field", ["first_rate", "last_rate", "worst_rate"])
def test_trend_rates_reproduce(fresh, committed, field):
    for assay in assays(committed):
        got = (fresh[assay].get("trend") or {}).get(field)
        want = (committed[assay].get("trend") or {}).get(field)
        if want is None and got is None:
            continue
        assert got is not None and want is not None, (
            f"{assay}.trend.{field}: committed {want!r}, fresh {got!r}"
        )
        assert abs(got - want) < 1e-9, (
            f"{assay}.trend.{field}: committed {want}, fresh {got}. Regenerate results/drift.json."
        )


def test_every_period_reproduces(fresh, committed):
    """Not just the endpoints - a stale middle is as wrong as a stale end."""
    for assay in assays(committed):
        a = committed[assay]["periods"]
        b = fresh[assay]["periods"]
        assert sorted(a) == sorted(b), f"{assay}: period set differs"
        for period, want in a.items():
            got = b[period]
            for key in ("total", "usable", "perfect", "failing", "lost_oligo", "three_prime"):
                assert got.get(key) == want.get(key), (
                    f"{assay} {period} {key}: committed {want.get(key)}, fresh {got.get(key)}"
                )


def test_the_clustered_rate_is_present_in_the_committed_results(committed):
    """A field the code computes but the artefact predates is the same staleness bug."""
    for assay in assays(committed):
        for name, period in committed[assay]["periods"].items():
            assert "failing_clustered_rate" in period, (
                f"{assay} {name} has no failing_clustered_rate; "
                f"results/drift.json predates the clustered rule - regenerate it."
            )


@pytest.mark.skipif(not DRIFT_BIG.exists(), reason="the big results are not shipped")
def test_the_big_results_come_from_the_code_in_this_repository():
    """The 291,827-genome run takes hours, so compare a digest rather than rerun it.

    Three versions of this guard, and the first two were wrong in opposite directions.

    It began as "the artefact is not older than the last commit touching the analysis
    sources". That passed in any working tree where the code had been edited and not yet
    committed - the state in which the artefact is *most* likely to be stale. Dating the
    dirty files by their mtime closed that, and then overshot: running this repository's
    own formatter over the tree changed three files' mtimes without changing a single
    expression, and the guard then demanded a three-hour run to republish byte-identical
    numbers. A guard that cries wolf over whitespace gets switched off.

    So the comparison is a digest of the sources' *syntax trees*, recorded in the artefact
    when it was generated. Reformatting leaves it unchanged; editing a threshold, a
    comparison or a docstring does not.
    """
    run = json.loads(DRIFT_BIG.read_text(encoding="utf-8"))["_run"]
    recorded = run.get("analysis_source_digest")
    if recorded is None:
        pytest.fail(
            "results/drift_big.json carries no analysis_source_digest, so nothing says "
            "whether its numbers come from this code. Regenerate it with "
            "  python scripts/analyse.py --data data/sequences_big.jsonl "
            "--json results/drift_big.json --mutations-json results/mutations_big.json"
        )

    sys.path.insert(0, str(ROOT / "scripts"))
    from analyse import ANALYSIS_SOURCES, source_digest  # noqa: PLC0415

    current = source_digest(ROOT)
    assert current is not None, (
        f"one of {', '.join(ANALYSIS_SOURCES)} could not be read or parsed, so the "
        f"digest cannot be computed - which is not a pass"
    )
    assert current == recorded, (
        f"results/drift_big.json was generated from analysis code whose syntax digest "
        f"was {recorded[:12]}; the code here digests to {current[:12]}, so the published "
        f"big-corpus numbers do not come from it. Regenerate with "
        f"  python scripts/analyse.py --data data/sequences_big.jsonl "
        f"--json results/drift_big.json --mutations-json results/mutations_big.json"
    )


def test_the_digest_ignores_formatting_and_notices_a_changed_threshold(tmp_path):
    """Both directions of the claim the guard above rests on.

    A digest that moved on whitespace would reintroduce the false alarm it was rewritten
    to remove; one that did not move on a different constant would be worthless. Checked
    on copies in a temporary directory, so neither the repository nor the published
    artefact is touched.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    from analyse import ANALYSIS_SOURCES, source_digest  # noqa: PLC0415

    def plant(tag, edit):
        root = tmp_path / tag
        for name in ANALYSIS_SOURCES:
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(
                edit(name, (ROOT / name).read_text(encoding="utf-8")), encoding="utf-8"
            )
        return source_digest(root)

    same = lambda name, text: text  # noqa: E731
    original = plant("plain", same)
    assert original == source_digest(ROOT), "copying the files changed the digest"

    # Blank lines and trailing spaces: the two things the formatter moved, written as
    # chr(10) so the escape survives every editor and shell this file passes through.
    nl = chr(10)
    reformatted = plant(
        "reformatted",
        lambda name, text: text.replace(nl * 3, nl * 4).replace(
            "    return " + nl, "    return" + nl
        ),
    )
    assert reformatted == original, (
        "the digest moved on blank lines, which is the false alarm this guard exists to avoid"
    )

    def bump(name, text):
        if name.endswith("analyze.py"):
            return text.replace("SEVERE_MISMATCH_LOAD = 3", "SEVERE_MISMATCH_LOAD = 4", 1)
        return text

    assert plant("bumped", bump) != original, (
        "the digest did not move when SEVERE_MISMATCH_LOAD went from 3 to 4, which moves "
        "every failing rate in the output"
    )

    def gut(name, text):
        if name.endswith("analyze.py"):
            return text.replace(
                "return bool(self.hits) and all(h.usable for h in self.hits.values())",
                "return True",
                1,
            )
        return text

    gutted = ROOT / "src" / "assaydrift" / "analyze.py"
    assert (
        "return bool(self.hits) and all(h.usable for h in self.hits.values())"
        in gutted.read_text(encoding="utf-8")
    ), "the line this test plants over has moved; the plant would be a no-op"
    assert plant("gutted", gut) != original, (
        "the digest did not move when AssayResult.usable was replaced with True"
    )
