# assay-drift next to a general-purpose aligner

[Biopython](https://biopython.org/)'s `Bio.Align.PairwiseAligner` is the aligner most
people reach for when they need to know where a short sequence sits in a long one. It is
the right thing to compare this tool against, because the question underneath every figure
in this repository is exactly that: **where does this oligo sit in this genome, and how
many bases disagree?**

They answer it differently. `assaydrift.match` searches unambiguous seeds, scores a window
with IUPAC codes read as sets on *both* sides, and refuses a site whose mismatches exceed a
fifth of the oligo's length. `PairwiseAligner` runs a local alignment with a substitution
matrix and affine gaps and returns the best-scoring one, wherever it is.

Measured 2026-10-29, Biopython 1.88, seed 20261029: **15 oligos** (the five published
assays) against **20 genomes** from the fetch cache, so no network - **300 pairs**. Every
number here is in [`aligner-comparison.json`](aligner-comparison.json), written by
`scripts/compare_aligner.py`, and `tests/test_aligner_comparison.py` checks this page
against it.

## The result

| | |
|---|---:|
| pairs | 300 |
| both put the oligo at the same base | **270** |
| …and counted the same mismatches | **233** |
| disagreed | 67 |
| sites found only by Biopython | **0** |
| Biopython left part of the oligo unaligned | 33 |
| assay-drift / Biopython, seconds | **38** / 231 |

Biopython never found a binding site this tool missed. The 67 disagreements fall on **4 of
the 15 oligos**, and they split cleanly in two.

## 37 of them are one thing: a published primer's degenerate bases

Three of the fifteen oligos in the catalogue carry IUPAC degenerate bases, because their
designers put them there:

| oligo | sequence | degenerate |
|---|---|---|
| `RdRp-F` | `GTGARATGGTCATGTGTGGCGG` | `R` |
| `RdRp-R` | `CARATGTTAAASACACTATTAGCATA` | `R`, `S` |
| `N3-P` | `AYCACATTGGCACCCGCAATCCTG` | `Y` |

`R` means "A or G". A subject base of `A` under that `R` is what the primer was written to
accept, so this tool scores it a **match** - not even an ambiguity. A substitution matrix
of +1/-1 has no notion of a code standing for a set, so `R` against `A` is two different
letters and costs a mismatch.

| | assay-drift | Biopython | rows |
|---|---:|---:|---:|
| `RdRp-F` at its real site | 0 mm | **1 mm** | 17 |
| `RdRp-R` at its real site | 1 mm | **2 mm** | 20 |

In **all 37** of these Biopython is *higher*, never lower, and the difference never exceeds
that oligo's degenerate-base count. The artefact carries a third column to prove the cause
rather than assert it: at the site this tool reports, the oligo compared to the window with
**plain letter equality** - no IUPAC, no gaps, no scoring. **All 37** of Biopython's counts
equal that plain count. One variable between the two tools, and it is this one.

That matters here more than it would elsewhere. Charité RdRp is the assay carrying this
repository's sharpest result - 86.8% exact match at the end of the series, and 94.9% likely
failing in 2023-Q3 - and it is the assay with the most degenerate bases. An aligner that
counts a designed degeneracy as drift inflates the signal on precisely the assay whose
drift matters most.

## 30 of them are the other thing: a trimmed end

Local alignment returns the best-scoring *subsequence*, so it can buy a better score by
not aligning a base that disagrees. It did, 30 times, always by exactly two bases:

| oligo | assay-drift | Biopython | direction | rows |
|---|---:|---:|---|---:|
| `E-F` | 1 mm at 26205 | **0 mm** at 26207 | **understates** | 10 |
| `N3-P` | 0 mm at 28631 | 0 mm at 28633 | equal | 20 |

`E-F` is the one to look at. It carries no degenerate base; there is a real mismatch near
its end, and the aligner reported a clean match by leaving that end out of the alignment
and starting two bases later. **For a primer the 3′ end is the business end** - it is where
the polymerase extends, which is why this tool applies a separate rule to the last five
bases of a primer and not to a probe. Trimming is leniency in the one place leniency is
least defensible.

Unconstrained, it goes further than two bases. Before the comparison was corrected (below),
free local alignment put `N1-R` **24,800 bases from its real site**, with 9 of its 20 bases
unaligned, because a short near-perfect partial match outscores a full-length match
carrying mismatches. That is a correct local alignment and a useless binding site.

## Three corrections to this comparison, before any of it counted

Each of these was a bug in the harness, not in either tool, and each only showed up by
running the thing end to end.

- **The reverse primers were not the same question on both sides.** `find_oligo`
  reverse-complements a reverse primer, because that is the strand a forward-strand genome
  record carries; Biopython was handed the oligo as written. Five of the fifteen oligos are
  reverse primers and every one of their "disagreements" was the harness asking two
  different questions. Agreement went from 17 of 30 to 27 of 30 on the pilot once both
  sides searched for the same string.
- **An "anchored" second pass was published-adjacent and wrong.** It was meant to force
  every base of the oligo to align, and it returned **19 mismatches of 24** for a probe
  this tool scores at 0 - it anchored on a position trimming had already shifted, then
  counted columns of a global alignment that is degenerate inside a window of the oligo's
  own length. It was deleted rather than reported, and replaced by the plain-equality
  column above, which needs no aligner and is checkable by eye.
- **Nine megabytes of genome went down a pipe.** Each of the 15 oligos carried its own copy
  of the same 30,000-base record, stdin truncated at 300 pairs, and the other interpreter
  died on `Unterminated string ... column 82` - a message about JSON that was really a
  message about a pipe. Genomes are sent once now, keyed, in a file.

## What this does not show

- **Biopython is not a drift tool and is not trying to be.** It returns an alignment. It
  makes no verdict, applies no 3′-end rule, draws no line between a mismatch and an `N`
  read as missing data, and knows nothing about assays. Every column it is missing here is
  a column it never claimed.
- **Two of its parameters were chosen by me.** `open_gap_score` and `extend_gap_score` are
  set to -10 so that the aligner prefers a substitution over a gap, because this tool
  counts mismatches. A different choice moves its numbers, and a gap-tolerant setting would
  change the trimming behaviour above. The values are in the script and in the artefact.
- **20 genomes, not a corpus.** Drawn from the fetch cache with a fixed seed. The
  291,827-genome run in the README is this tool measuring; this page is two tools being
  compared on a sample small enough to inspect by hand.
- **Nothing here is a recall comparison.** Biopython found no site this tool missed, which
  is a fact about these 300 pairs and not a bound on anything. The tools do not share a
  definition of "the binding site" - that is the finding, not a nuisance to be averaged out.
