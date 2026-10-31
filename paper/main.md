---
title: "Where a Mutation Lands: Positional Drift in SARS-CoV-2 PCR Assays, and How Much of It Is the Rule"
author: "Muhammad Hammas, Independent Researcher"
date: 2026-10-31
---

# Abstract

A PCR assay fails silently when its target mutates under an oligo, because a failed
amplification looks exactly like a negative sample. Public tools exist to watch for this —
NEB's Primer Monitor, assayM, an automated validation website, PCR_strainer — and all of
them report the **current** mismatch picture against a sequence database.

We align five published SARS-CoV-2 assays against **291,827 dated GenBank genomes** from
2020 to 2026 and report by collection quarter, with a Kish effective sample size per
quarter so that a rate resting on one submission is visible as such. Two results follow
from reporting over time rather than at a point. **Two assays lost their exact match
entirely** (99.9% and 98.4% down to 0.0%), and neither mutation sits where it would stop
the reaction. And **Charité RdRp ends the series at 86.8% exact match while having been
5.1% exact, with 94.9% of sequences carrying a 3′-terminal mutation, in 2023-Q3** — an
assay that looks healthy at both ends of a series and was unusable in the middle. An
endpoint-only report cannot see that; ours would not have either, before the worst
quarter was reported beside the endpoints.

A 2025 wet-lab study argues that in-silico predictions of PCR failure overstate the risk.
We take that criticism as a measurement problem and apply it to our own rule. The result
is unflattering and is published as a headline rather than a caveat: the only clause a
bench would not dispute — an oligo with no binding site, which cannot prime — **fires for
no sequence in the corpus**, and for four of the five assays the entire `likely_failing`
rate comes from a single uncalibrated clause about a primer's last five bases.

# 1. Introduction

Existing primer-monitoring tools answer "which variants currently carry a mismatch under
this oligo". That is the right question for a laboratory choosing an assay this week, and
the wrong one for deciding whether an assay was dependable over a period. The two come
apart, and the gap is this paper's subject.

**Related work.** An automated validation website screens deployed assays against genome
databases using ThermonucleotideBLAST, reporting per-assay recall and mismatch
heatmaps [@bioinformatics2020]. assayM maps primer and probe regions onto a genome browser
and counts strains with mutations in those regions [@assaym]. NEB's Primer Monitor tracks
GISAID sequences against user-supplied primer sets, released after FDA guidance urged
developers to monitor primer alignment [@primermonitor]. PCR_strainer examined 15
primer/probe sets from eight laboratory-developed tests and found most with no mismatches
in over 98% of genomes, with three sets problematic [@pcrstrainer]. A 2025 wet-lab study
tested assays in vitro and found most performed well despite mismatches, measured by PCR
efficiency and Ct shift — concluding that in-silico false-negative predictions may
overstate risk [@wetlab2025].

None of those reports a time series with a stated independence structure, and the last of
them is a direct challenge to the category this tool's most alarming verdict belongs to.

**Contributions.**

1. A quarterly series over 291,827 dated genomes with **Kish effective sample size** per
   quarter, and an independence floor an endpoint must clear (§2, §3).
2. The finding that **endpoint reporting hides the middle**, with Charité RdRp as the case
   (§3).
3. A **decomposition of our own failure rule** showing that its bench-defensible clause
   never fires, taken directly from the wet-lab criticism (§4).

# 2. Method

Each oligo is aligned with IUPAC codes read as sets **on both sides**: a primer `R` matches
`A` and `G`, and a subject `R` where the primer wants `A` is ambiguous, not a mismatch.
Overstating drift is the error this work has had to correct more than once.

An unknown base under an oligo is **missing data**, excluded from rates and reported, so a
change in sequencing quality does not appear as viral drift. A record with no usable
collection date is dropped and counted, never defaulted to today.

Each sequence is scored `perfect`, `drifted` (mismatched, expected to amplify),
`likely_failing`, or `unknown`. The 3′ window applies only to primers: a hydrolysis probe
is blocked at its 3′ end and never extended.

**Independence.** Sequences from one submission are not independent observations, so each
quarter is grouped by (country, accession block) and given a Kish effective sample size
$n_\text{eff} = n^2 / \sum n_i^2$. A quarter needs 25 usable sequences **and** 10
independent submitter groups to carry an endpoint. Median effective $n$ across the 28
reliable quarters is 82.3, against about 4 on a 2,765-genome subset — which is why the
large corpus exists. Without the floor, the series would have ended at 2026-Q4: 67 usable
sequences from 7 groups, $n_\text{eff} = 1.5$, and CDC N2 reported as down 0.9 points
instead of 19.3.

# 3. Results

Exact match, first to last quarter able to carry an endpoint:

| Assay | 2020-Q1 | 2026-Q3 | Change | Likely failing at the end | Worst quarter |
|---|---:|---:|---:|---:|---|
| Charité E | 99.9% | **0.0%** | −99.9 | 0.0% | 0.2% (2023-Q3) |
| CDC N1 | 98.4% | **0.0%** | −98.4 | 5.8% | 5.8% (2026-Q3) |
| CDC N2 | 99.4% | 80.1% | −19.3 | 7.6% | **30.4% (2026-Q2)** |
| Charité RdRp | 99.9% | 86.8% | −13.0 | 13.1% | **94.9% (2023-Q3)** |
| CDC N3 *(retired)* | 98.8% | 86.1% | −12.7 | 0.7% | 0.8% (2026-Q2) |

**Charité RdRp is the result.** Its endpoints say 86.8% exact and 13.1% likely failing —
slow erosion. In 2023-Q3 it was 5.1% exact and 94.9% likely failing over 506 usable
sequences from 57 submitter groups, so not an artefact of one submission. Delta's NSP12
G671S sits at the 3′ terminus of its forward primer; the rate recovered when Omicron
displaced Delta. **An assay can be fine at both ends of a series and have been unusable in
the middle**, which is what an endpoint-only tool cannot report and why the worst
qualifying quarter is published beside the endpoints.

The same four oligos, at the same base of each, appear independently in a corpus 105 times
smaller, with shares differing because the year mixes differ — not because the finding
moved.

# 4. How much of `likely_failing` is the rule

@wetlab2025 found assays mostly working despite mismatches and concluded in-silico
predictions may overstate risk. `likely_failing` is such a prediction. Its rule has three
clauses: (1) an oligo with no binding site, (2) a mismatch within 5 bases of a primer's 3′
end, (3) ≥3 mismatches across the assay. Clause 1 is the only one a bench would not
dispute. Clauses 2 and 3 carry a threshold each, and **neither came from a wet-lab
calibration**: 5 bases is a primer-design rule of thumb, 3 is a judgement.

Re-scoring the 2,765-genome corpus across 3 window values × 4 load values × 3 rule shapes:

| Assay | published | without clause 2 | clause 1 alone |
|---|---:|---:|---:|
| CDC N1 | 0.0239 | 0.0235 | **0.0000** |
| CDC N2 | 0.0245 | 0.0000 | **0.0000** |
| CDC N3 | 0.0058 | 0.0000 | **0.0000** |
| Charité E | 0.0011 | 0.0000 | **0.0000** |
| Charité RdRp | 0.1116 | 0.0004 | **0.0000** |

Two things follow. **The conservative floor is zero for all five assays** — no sequence in
this corpus loses an oligo outright, so every `likely_failing` verdict is heuristic. And
**for four of the five, clause 2 is the entire rate**; Charité RdRp loses 99.6% of its
rate without it.

The category is kept rather than withdrawn, because reporting only exact-match loss
discards the positional information that is the point: where a mutation lands matters more
than how much of it there is, and G15451A one base from a 3′ terminus is a real, named
observation. What changes is the labelling. `likely_failing` is a **prediction under a
stated rule**, its decomposition is published beside it, and @wetlab2025 is cited at the
point of use.

# 5. Threats to validity

**The failure rule is uncalibrated** (§4). This is the dominant threat and it is the
paper's own section rather than a line here.

**The corpus is GenBank, not a clinical sample set.** Submission is uneven by country and
institution; the Kish correction addresses the arithmetic of that, not the selection.

**Short and partial records.** A record shorter than an oligo has no binding site by
construction, which is not the same fact as a site having been lost; those are `unknown`.

**One alignment implementation.** The matcher is compared against Biopython's
`PairwiseAligner` on 300 oligo–genome pairs: 270 agree on the binding site, and every
full-length disagreement is a published primer's deliberate degenerate base that a
substitution matrix counts as a mismatch. That validates the aligner, not the rule.

**The sensitivity grid is the small corpus.** One pass over 291,827 genomes is about three
hours and the grid is twelve passes. The direction is reported; the magnitudes are the
2,765-genome subset's.

# 6. Reproducibility

Both corpora rebuild from committed manifests — `data/accessions.tsv` and
`data/accessions_big.tsv.gz` — because a published rate whose corpus cannot be rebuilt is
not a result. `tests/test_published_numbers.py` checks every figure in the README against
the result file it came from, per corpus. The figures and the grid rebuild with
`research/sensitivity.py`.

# References

- [@bioinformatics2020] *A Public Website for the Automated Assessment and Validation of SARS-CoV-2 Diagnostic PCR Assays.* Bioinformatics 37(7):1024.
- [@assaym] *assayM: a web application to monitor mutations in COVID-19 diagnostic assays.* bioRxiv 2020.12.18.423467.
- [@primermonitor] *Primer Monitor: an online tool to track SARS-CoV-2 variants that may impact primers used in diagnostic assays.* New England Biolabs.
- [@pcrstrainer] *PCR_strainer* — assessment of 15 primer/probe sets from eight laboratory-developed tests.
- [@wetlab2025] Front. Cell. Infect. Microbiol. 2025, 10.3389/fcimb.2025.1524025 — in vitro performance despite primer and probe mismatches.

> **Citation status.** All five were located by literature search and their claims are
> quoted from search summaries rather than from the full texts. @wetlab2025 is
> load-bearing for §4 and must be read in full before submission — its exact conditions,
> assays and effect sizes determine how strongly §4's conclusion can be stated.
