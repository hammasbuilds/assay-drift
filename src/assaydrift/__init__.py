"""Do published PCR diagnostic assays still match what is circulating?

A diagnostic PCR test works by binding short synthetic DNA sequences - primers
and a probe - to a specific stretch of the target organism's genome. When the
organism mutates under those binding sites, the test starts returning false
negatives, and it does so silently: a failed amplification looks exactly like a
negative sample.

This measures that, from public data. It downloads dated genomes from GenBank,
aligns published assay oligos against each one, and reports per-quarter how many
sequences each assay would still have detected.
"""

from .analyze import AssayResult, PeriodSummary, evaluate, period_of, summarise, trend
from .match import Hit, find, find_oligo, reverse_complement
from .ncbi import NCBIError, Record, count, fetch, search
from .primers import ASSAYS, REFERENCE, Assay, Oligo, all_oligos, by_name

__version__ = "0.1.0"

__all__ = [
    "ASSAYS",
    "REFERENCE",
    "Assay",
    "AssayResult",
    "Hit",
    "NCBIError",
    "Oligo",
    "PeriodSummary",
    "Record",
    "all_oligos",
    "by_name",
    "count",
    "evaluate",
    "fetch",
    "find",
    "find_oligo",
    "period_of",
    "reverse_complement",
    "search",
    "summarise",
    "trend",
]
