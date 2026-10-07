"""The functions `suite-auditor coverage` found no test reaching at all.

Twelve of them, and the first three matter most: `PeriodSummary.effective_n`, `.clusters`
and `.median_mismatches` are quoted per quarter in the README and in `results/drift.json`,
and nothing exercised the arithmetic behind them. A published number computed by untested
code is the shape of error this repository has had to correct twice.

The rest is `ncbi.py`, which the audit reported at 35% of statements: a network layer is
awkward to test and `_get` is the one seam it needs, so the three functions above it are
driven through a stub. `read_jsonl`, `Record.year`, `_slot` and `Assay.__iter__` are pure
and need nothing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from assaydrift import ncbi
from assaydrift.analyze import PeriodSummary
from assaydrift.ncbi import DataError, Record, read_jsonl
from assaydrift.primers import ASSAYS


def _summary(**kwargs) -> PeriodSummary:
    """A PeriodSummary with only the fields a given property reads."""
    base = {"period": "2021-Q1", "total": 0, "usable": 0}
    base.update(kwargs)
    return PeriodSummary(**base)


class TestEffectiveN:
    """Kish effective sample size: n^2 / sum of squared group sizes.

    The README quotes it beside every rate, and its whole purpose is to say that 159
    sequences from four submissions are not 159 observations. Nothing tested it.
    """

    def test_one_group_is_worth_one_observation(self) -> None:
        assert _summary(submitters={"a": 159}).effective_n == pytest.approx(1.0)

    def test_equal_groups_are_worth_their_count(self) -> None:
        assert _summary(submitters={"a": 10, "b": 10, "c": 10}).effective_n == pytest.approx(3.0)
        assert _summary(submitters={f"g{i}": 5 for i in range(7)}).effective_n == pytest.approx(7.0)

    def test_one_dominant_group_pulls_it_towards_one(self) -> None:
        """The case the measure exists for: a quarter that is mostly one submission."""
        lopsided = _summary(submitters={"big": 155, "a": 1, "b": 1, "c": 1, "d": 1})
        assert lopsided.clusters == 5
        assert lopsided.effective_n == pytest.approx(1.05, abs=0.01)
        assert lopsided.effective_n < 1.1, (
            "159 sequences that are 97% one submission must not read as 5 observations"
        )

    def test_it_never_exceeds_the_number_of_groups_nor_falls_below_one(self) -> None:
        for sizes in ([1], [1, 1], [100, 1], [3, 4, 5], [1] * 20, [50, 30, 20]):
            s = _summary(submitters={f"g{i}": n for i, n in enumerate(sizes)})
            assert 1.0 - 1e-9 <= s.effective_n <= len(sizes) + 1e-9, sizes

    def test_no_submitters_is_none_not_zero(self) -> None:
        """None means "not computed"; 0.0 would read as "no independent observations"."""
        assert _summary(submitters={}).effective_n is None
        assert _summary().effective_n is None

    def test_clusters_counts_groups_not_records(self) -> None:
        assert _summary(submitters={"a": 100, "b": 1}).clusters == 2
        assert _summary(submitters={}).clusters == 0


class TestMedianMismatches:
    def test_the_median_of_the_counts(self) -> None:
        assert _summary(mismatch_counts=[0, 1, 2]).median_mismatches == 1
        assert _summary(mismatch_counts=[0, 0, 1, 3]).median_mismatches == pytest.approx(0.5)

    def test_no_counts_is_none(self) -> None:
        """Not 0: "no mismatches measured" and "a median of zero mismatches" are
        different statements, and one of them would read as a clean assay."""
        assert _summary(mismatch_counts=[]).median_mismatches is None
        assert _summary().median_mismatches is None


class TestRecordYear:
    @pytest.mark.parametrize(
        ("collected", "expected"),
        [
            ("2021-03-14", 2021),
            ("14-Mar-2021", 2021),
            ("2021", 2021),
            ("1999-01", 1999),
            ("", None),
            ("not a date", None),
            ("18-Mar-99", None),  # a two-digit year is not one this study can place
        ],
    )
    def test_the_year_it_can_read(self, collected: str, expected: int | None) -> None:
        record = Record(
            accession="X1", organism="o", sequence="ACGT", collected=collected, length=4
        )
        assert record.year == expected


class TestReadJsonl:
    def _write(self, tmp_path: Path, *lines: str) -> Path:
        path = tmp_path / "seq.jsonl"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def _good(self, **over) -> str:
        row = {
            "accession": "MW1",
            "organism": "SARS-CoV-2",
            "sequence": "ACGT",
            "collected": "2021-03",
        }
        row.update(over)
        return json.dumps(row)

    def test_a_missing_file_says_what_to_run(self, tmp_path: Path) -> None:
        with pytest.raises(DataError, match="run scripts/fetch.py"):
            list(read_jsonl(tmp_path / "nothing.jsonl"))

    def test_blank_lines_are_skipped(self, tmp_path: Path) -> None:
        path = self._write(tmp_path, self._good(), "", "   ", self._good(accession="MW2"))
        assert [r.accession for r in read_jsonl(path)] == ["MW1", "MW2"]

    def test_a_bad_line_names_the_file_and_the_line(self, tmp_path: Path) -> None:
        path = self._write(tmp_path, self._good(), "{not json")
        with pytest.raises(DataError) as caught:
            list(read_jsonl(path))
        assert "seq.jsonl:2" in str(caught.value)

    def test_on_bad_line_keeps_the_good_records(self, tmp_path: Path) -> None:
        """A single truncated line must not cost the other two thousand."""
        path = self._write(
            tmp_path, self._good(), "{not json", self._good(accession="MW3")
        )
        seen: list[str] = []
        got = list(read_jsonl(path, on_bad_line=seen.append))
        assert [r.accession for r in got] == ["MW1", "MW3"]
        # `on_bad_line` is handed the DataError itself, not its text - which is
        # more useful to a caller and is what this test assumed wrongly first.
        assert len(seen) == 1
        assert isinstance(seen[0], DataError)
        assert "seq.jsonl:2" in str(seen[0])

    @pytest.mark.parametrize("field", ["accession", "organism", "sequence", "collected"])
    def test_a_required_field_is_required(self, tmp_path: Path, field: str) -> None:
        """Defaulting any of these is how a drift study invents a trend."""
        row = json.loads(self._good())
        del row[field]
        path = self._write(tmp_path, json.dumps(row))
        with pytest.raises(DataError, match=field):
            list(read_jsonl(path))

    def test_a_line_that_is_not_an_object_is_rejected(self, tmp_path: Path) -> None:
        path = self._write(tmp_path, "[1, 2, 3]")
        with pytest.raises(DataError):
            list(read_jsonl(path))

    def test_an_empty_file_is_an_error_not_an_empty_study(self, tmp_path: Path) -> None:
        path = tmp_path / "seq.jsonl"
        path.write_text("", encoding="utf-8")
        with pytest.raises(DataError):
            list(read_jsonl(path))


class TestNetworkLayerThroughItsOneSeam:
    """`count`, `search` and `fetch` all go through `_get`, so that is the stub point."""

    def test_count_reads_the_count_element(self, monkeypatch) -> None:
        reply = "<eSearchResult><Count>4821</Count></eSearchResult>"
        monkeypatch.setattr(ncbi, "_get", lambda url, **kw: reply)
        assert ncbi.count("anything") == 4821

    def test_count_with_no_count_element_is_zero_not_a_crash(self, monkeypatch) -> None:
        monkeypatch.setattr(ncbi, "_get", lambda url, **kw: "<eSearchResult/>")
        assert ncbi.count("anything") == 0

    def test_search_pages_until_it_has_enough(self, monkeypatch) -> None:
        pages = [
            "".join(f"<Id>{i}</Id>" for i in range(500)),
            "".join(f"<Id>{i}</Id>" for i in range(500, 700)),
        ]
        seen: list[str] = []

        def fake(url, **kw):
            seen.append(url)
            return pages[len(seen) - 1] if len(seen) <= len(pages) else ""

        monkeypatch.setattr(ncbi, "_get", fake)
        ids = ncbi.search("term", limit=700)
        assert len(ids) == 700
        assert ids[0] == "0" and ids[-1] == "699"
        assert "sort=date" in seen[0], "a capped run must sample the NEWEST records"

    def test_search_stops_on_an_empty_page(self, monkeypatch) -> None:
        monkeypatch.setattr(ncbi, "_get", lambda url, **kw: "")
        assert ncbi.search("term", limit=500) == []

    def test_search_never_returns_more_than_asked_for(self, monkeypatch) -> None:
        monkeypatch.setattr(
            ncbi, "_get", lambda url, **kw: "".join(f"<Id>{i}</Id>" for i in range(500))
        )
        assert len(ncbi.search("term", limit=10)) == 10


def test_an_assay_iterates_its_own_oligos() -> None:
    for assay in ASSAYS.values() if isinstance(ASSAYS, dict) else ASSAYS:
        oligos = list(assay)
        assert oligos, f"{assay} iterated to nothing"
        assert all(hasattr(o, "sequence") and hasattr(o, "role") for o in oligos)
        assert len(oligos) == len(assay.oligos)
