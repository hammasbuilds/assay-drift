"""Which base changed, where in the oligo, and how often.

A rate says an assay drifted. It does not say what an assay designer needs to
know next, which is whether to shift a primer by six bases or redesign the
amplicon - and that depends entirely on *where* in the oligo the change landed
and how close it is to the 3' terminus.

So this reports, per oligo, every position that disagrees with the reference
genome: the base in the oligo, the genome coordinate, the substitution, how far
it is from the 3' end, and the share of sequences in each period carrying it.
The names that come out - C28311T, C26270T, G15451A - are produced here, not
looked up; the tool knows nothing about variants.

Two details decide whether a coordinate means anything:

**Coordinates are the reference's, not the sample's.** Each oligo is located
once in NC_045512.2, and a mismatch at offset *j* of that alignment is named at
reference position `start + j`. Naming it from the sample's own alignment would
shift every coordinate in a genome carrying an upstream deletion, so the same
mutation would be reported under several names.

**`N` is never a substitution.** An unknown base is missing data; counting it
as a change would name a mutation that is really a sequencing gap, and gaps are
not evenly spread over time.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from .analyze import period_of
from .match import IUPAC, THREE_PRIME_WINDOW, find_oligo, reverse_complement
from .ncbi import Record
from .primers import ASSAYS, Assay, Oligo
from .reference import genome


@dataclass(frozen=True)
class Site:
    """One position of one oligo, as it sits on the reference genome."""

    assay: str
    oligo: str
    role: str
    oligo_length: int
    # 1-based position within the oligo as written, 5' to 3'.
    base_in_oligo: int
    bases_from_3_prime_end: int
    genome_position: int  # 1-based, on the reference
    reference_base: str

    @property
    def in_three_prime_window(self) -> bool:
        """Close enough to the 3' terminus to stop extension - primers only."""
        return self.role != "probe" and self.bases_from_3_prime_end < THREE_PRIME_WINDOW


@dataclass
class Observation:
    """One substitution at one site, counted over periods."""

    site: Site
    alt_base: str
    by_period: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    @property
    def name(self) -> str:
        """The mutation as a geneticist would write it: C28311T."""
        return f"{self.site.reference_base}{self.site.genome_position}{self.alt_base}"

    @property
    def total(self) -> int:
        return sum(self.by_period.values())


def _reference_sites(oligo: Oligo, assay: Assay) -> list[Site]:
    """Where each base of this oligo sits on the reference genome.

    Index *j* below runs along the *search string* - the oligo as written for a
    forward primer or probe, its reverse complement for a reverse primer - so it
    is also the offset along the genome. The oligo's own numbering runs the
    other way for a reverse primer, which is the same end-for-end flip that
    decides where the 3' window goes.
    """
    hit = find_oligo(genome(), oligo.sequence, oligo.role)
    if not hit.found:
        return []
    length = len(oligo.sequence)
    reverse = oligo.role == "reverse"
    reference = genome()
    sites = []
    for j in range(length):
        sites.append(
            Site(
                assay=assay.name,
                oligo=oligo.name,
                role=oligo.role,
                oligo_length=length,
                base_in_oligo=(length - j) if reverse else (j + 1),
                bases_from_3_prime_end=j if reverse else (length - 1 - j),
                genome_position=hit.position + j + 1,
                reference_base=reference[hit.position + j],
            )
        )
    return sites


def reference_sites(assays: list[Assay] | None = None) -> dict[str, list[Site]]:
    """Per-oligo reference coordinates, keyed by oligo name."""
    return {
        oligo.name: _reference_sites(oligo, assay)
        for assay in (assays if assays is not None else ASSAYS)
        for oligo in assay.oligos
    }


def _search_string(oligo: Oligo) -> str:
    return reverse_complement(oligo.sequence) if oligo.role == "reverse" else oligo.sequence


def observe(
    records: list[Record] | None = None,
    assays: list[Assay] | None = None,
    granularity: str = "year",
    *,
    stream=None,
) -> tuple[dict[tuple[str, int, str], Observation], dict[str, int]]:
    """Count every substitution under every oligo, by period.

    Returns the observations keyed by (oligo name, position in the search
    string, alternative base), and the number of *usable* sequences per period
    per assay - the denominator a rate has to be taken over, since a sequence
    with an unknown base under an oligo is no evidence either way.

    `records` may be a list or, via `stream`, any iterable of Records; the
    caller decides whether 80 MB of sequence is held in memory.
    """
    assays = assays if assays is not None else ASSAYS
    sites = reference_sites(assays)
    searched = {oligo.name: _search_string(oligo) for assay in assays for oligo in assay.oligos}

    observations: dict[tuple[str, int, str], Observation] = {}
    denominators: dict[str, int] = defaultdict(int)

    source = stream if stream is not None else (records or [])
    for record in source:
        period = period_of(record.collected, granularity)
        if period is None:
            continue
        for assay in assays:
            hits = {
                oligo.name: find_oligo(record.sequence, oligo.sequence, oligo.role)
                for oligo in assay.oligos
            }
            if not all(hit.usable for hit in hits.values()):
                continue  # unknown bases under an oligo: no evidence either way
            denominators[f"{assay.name}|{period}"] += 1
            for oligo in assay.oligos:
                hit = hits[oligo.name]
                search = searched[oligo.name]
                window = record.sequence[hit.position : hit.position + len(search)].upper()
                for j, (base, actual) in enumerate(zip(search, window, strict=True)):
                    site = sites[oligo.name][j]
                    # A sample that still reads as the 2019 reference has not
                    # drifted, even where the oligo itself was designed with a
                    # known mismatch against that reference (Charite RdRp's
                    # reverse primer, "S vs T"). Comparing only to the oligo's
                    # own base would report that design choice as a "mutation"
                    # carried by 100% of every sample, forever.
                    if actual == "N" or actual == site.reference_base or actual in IUPAC[base]:
                        continue
                    key = (oligo.name, j, actual)
                    if key not in observations:
                        observations[key] = Observation(site=site, alt_base=actual)
                    observations[key].by_period[period] += 1
    return observations, dict(denominators)


def rank(
    observations: dict[tuple[str, int, str], Observation],
    denominators: dict[str, int],
    period: str | None = None,
) -> list[tuple[Observation, float]]:
    """Observations with their rate, most frequent first.

    The rate is over usable sequences of the *same period*, so a mutation that
    appeared in 2021 is not diluted by six years of data that predate it.
    """
    out = []
    for observation in observations.values():
        if period is None:
            count = observation.total
            prefix = f"{observation.site.assay}|"
            total = sum(n for k, n in denominators.items() if k.startswith(prefix))
        else:
            count = observation.by_period.get(period, 0)
            total = denominators.get(f"{observation.site.assay}|{period}", 0)
        out.append((observation, count / total if total else 0.0))
    out.sort(key=lambda pair: (-pair[1], pair[0].site.genome_position))
    return out


def peak(observation: Observation, denominators: dict[str, int]) -> tuple[str, float]:
    """The period in which this substitution was most common, and its rate.

    Reported instead of an overall share because these are epidemic waves: a
    mutation carried by two thirds of one year's sequences and none of the next
    averages to something that describes neither.
    """
    best_period, best_rate = "", 0.0
    for period, count in observation.by_period.items():
        total = denominators.get(f"{observation.site.assay}|{period}", 0)
        rate = count / total if total else 0.0
        if rate > best_rate:
            best_period, best_rate = period, rate
    return best_period, best_rate
