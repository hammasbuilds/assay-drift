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

A 2025 wet-lab study [@wetlab2025] tested 16 assays against over 200 mismatched templates
in vitro. We take it as a measurement problem and apply it to our own rule. The result is
unflattering and is published as a headline rather than a caveat. The only clause a bench
would not dispute — an oligo with no binding site, which cannot prime — **fires for no
sequence in either corpus**: not one of the 2,765 genomes and not one of the 291,827,
once the single partial deposit that has no such region to lose is counted separately
(§5). Worse, the clause that produces almost all of the rate is
contradicted by the source it derives from: Lefever et al. [@lefever2013] license
tolerating single mismatches *more than* 5 bases from a 3′ end, which our rule inverts
into a failure prediction for single mismatches *within* 5 bases. Across 2,765 genomes
and five assays, **every one of the 395 3′-proximal hits is a single mismatch and none is
a pair** — and a single 3′-proximal mismatch is the one configuration Otwell et al.
measured directly, at about one cycle of Ct shift and no false negative. At 291,827
genomes the pair does appear, **52 times in 54,755 3′-proximal hits, 0.095% of them**
(§4.5), which puts a measured base rate under the surveillance claim instead of an
absence. We therefore
withdraw our most alarming number as a failure prediction — Charité RdRp's **18.29%** on
the 291,827-genome corpus, 11.16% on the 2,765-genome one, which fall to 0.00018 and
0.0004 when clause 2 requires the two mismatches the bench measured — and reframe the
category as **surveillance**: those genomes are one mutation away from a configuration
that does have wet-lab evidence behind it (§4, §4.5).

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
overstate risk [@wetlab2025]. A 2021 study re-evaluated eleven primer/probe sets against
84,305 genomes collected between June 2020 and January 2021 and reported the China CDC
N-gene assay's inclusivity continuing to fall as new mismatches appeared, the others
holding above 99% [@genes2021]. BioFire runs the only ongoing programme of this kind we
found: monthly in-silico screening that had covered 12,744,136 sequences by March 2023
[@biofire2023].

**What is new here, stated narrowly.** Time series of this kind exist — [@genes2021] is
one, over seven months and eleven assays, and [@biofire2023] is a continuing one. Neither
reports an **independence structure**: both count sequences, and sequences deposited in
one submission are not independent observations of a virus population. Neither publishes
a per-period effective sample size, and neither states a floor a period must clear before
it may carry an endpoint. Those are this paper's methodological contribution, over a
span — 2020 to 2026 — roughly ten times the longer of the two. The claim that no time
series existed appeared in an earlier draft and was false; it is the fourth absence-based
claim in this author's related work to fail on a second literature pass, and the pattern
is recorded in the citation note rather than the single correction.

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
slow erosion. In 2023-Q3 it was 5.1% exact and 94.9% likely failing over 469 usable
sequences from 57 submitter groups ($n_\text{eff} = 14.2$), so not an artefact of one
submission. Delta's NSP12
G671S sits at the 3′ terminus of its forward primer; the rate recovered when Omicron
displaced Delta. **An assay can be fine at both ends of a series and have been unusable in
the middle**, which is what an endpoint-only tool cannot report and why the worst
qualifying quarter is published beside the endpoints.

The same four oligos, at the same base of each, appear independently in a corpus 105 times
smaller, with shares differing because the year mixes differ — not because the finding
moved.

# 4. How much of `likely_failing` is the rule

This section was rewritten after its two load-bearing references were read in full rather
than quoted from search summaries. Doing so contradicted the tool's most alarming
verdict, and the contradiction is reported here rather than worked around.

`likely_failing` has three clauses: (1) an oligo with no binding site, (2) a mismatch
within 5 bases of a primer's 3′ end, (3) ≥3 mismatches across the assay. Clause 1 is the
only one a bench would not dispute.

## 4.1 What the wet-lab literature actually licenses

Lefever et al. [@lefever2013] is the origin of the 5-base window. Their conclusion is
that *"single mismatches located >5 bp from the 3′ end have a moderate effect on qPCR
amplification and can be tolerated."* That is a statement about mismatches **outside** the
window. It does not say that a single mismatch inside the window causes failure, and
clause 2 asserts the converse of what the source supports.

Their blocking threshold is also far above clause 3: *"≥4 mismatches in a single primer"*
blocks a reaction completely, as does three in one primer with two in the other. Clause 3
fires at three mismatches **anywhere across the assay**, which is below the published
threshold for a single primer.

Otwell et al. [@wetlab2025] then tested 16 assays against over 200 synthetic templates in
vitro, and their results separate position from count in a way clause 2 cannot:

| configuration tested | observed effect |
|---|---|
| single reverse-primer mismatch ~1 nt from the 3′ end (BVP 501Y) | ≈1 cycle Ct shift, **no false negative** |
| two forward-primer mismatches, same assay | 9–12 cycle shift, **genuine false negative** |
| three mismatches within the last four 3′ nucleotides (C3 ORF3a) | failure below 5,000 copies, Ct shift > 15 |
| up to four reverse-primer mismatches, not 3′-clustered (C4 ORF8) | 1–2 cycle shift |
| three 5′-end forward-primer mismatches (China_N) | comparable to control despite ≈7 °C predicted Tm drop |

A single 3′-proximal mismatch is tolerated; a cluster is not. Their overall conclusion is
also more balanced than this paper's earlier draft reported: they write that assays are
*"more robust to genetic drift than anticipated"* — but add *"or at least can be more
robust than predicted, depending on the specific wet lab conditions used"*, and they
record the **opposite** error too, templates not flagged in silico that performed poorly.
The earlier draft cited them for "in-silico predictions overstate risk" alone, which took
one half of a two-sided finding.

## 4.2 Re-scoring under a rule the evidence supports

A fourth rule shape was added: clause 2 requires **two or more** mismatches in the 3′
window rather than one. Re-scoring the 2,765-genome corpus, 3 windows × 4 loads × 4
shapes:

| Assay | published | **clustered 3′ (≥2)** | without clause 2 | clause 1 alone |
|---|---:|---:|---:|---:|
| CDC N1 | 0.0239 | 0.0235 | 0.0235 | **0.0000** |
| CDC N2 | 0.0245 | **0.0000** | 0.0000 | **0.0000** |
| CDC N3 | 0.0058 | **0.0000** | 0.0000 | **0.0000** |
| Charité E | 0.0011 | **0.0000** | 0.0000 | **0.0000** |
| Charité RdRp | 0.1116 | **0.0004** | 0.0004 | **0.0000** |

The clustered column is identical to dropping clause 2 altogether, which is not a
coincidence and not an artefact — at this corpus size. §4.5 repeats the re-scoring over
291,827 genomes, where the two columns part company. The distribution of 3′-proximal mismatch counts over all
2,765 genomes is:

| Assay | 0 in window | 1 | ≥2 |
|---|---:|---:|---:|
| CDC N1 | 2,718 | 1 | **0** |
| CDC N2 | 2,663 | 67 | **0** |
| CDC N3 | 2,745 | 16 | **0** |
| Charité E | 2,744 | 3 | **0** |
| Charité RdRp | 2,451 | 308 | **0** |

**No genome in this corpus carries two mismatches in a 3′ window.** Every one of the 395
3′-proximal hits — including all 308 that produce Charité RdRp's 11.16% — is a *single*
mismatch near a 3′ terminus. That is precisely the configuration Otwell et al. measured
directly and found to cost about one cycle.

At clause 3's load raised to 4, Lefever's single-primer blocking threshold, Charité
RdRp's rate is **0.0000 at every window value** (n = 2,759, 95% CI [0, 0.0014]).

## 4.3 What follows

The headline `likely_failing` rate for Charité RdRp is composed **entirely** of a
configuration for which the only direct in-vitro measurement available reports no false
negative. We therefore do not claim that 11.16% of genomes defeat the RdRp assay, and the
tool's own most alarming number is withdrawn as a failure prediction.

The category is kept rather than deleted, because the positional observation is real and
is the thing worth surveilling: G15451A one base from a 3′ terminus did occur, it did
spread, and a *second* mutation arriving in the same window would move those 308 genomes
into the configuration that Otwell et al. found does cause false negatives. The honest
framing is **surveillance, not prediction** — `likely_failing` marks assay–genome pairs
one mutation away from a configuration with wet-lab evidence behind it, and that is a
useful thing to watch and a different claim from the one the name makes.

Accordingly the tool now reports the clustered rate beside the published one, and this
paper's abstract and §3 are to be read with 4.3 in force. §4.5 carries the whole argument
to the 291,827-genome corpus, where the rate being withdrawn is 0.1829 rather than 0.1116
and the conclusion is the same.

## 4.4 Figures

All four rebuild via `python research/sensitivity.py && python research/run_all.py`.
Figures 1 to 3 are over the 2,765-genome corpus, because the sensitivity grid is twelve
three-hour passes; Figure 4 is the published setting over all 291,827 genomes, which
needs no grid.

![Figure 1](../research/figures/fig1_clause_decomposition.png)

**Figure 1.** Which clause of the `likely_failing` rule produces the published rate, at
the published settings (3′ window 5, mismatch load 3), about 2,719 usable sequences per
assay. The grey floor is the only clause a bench would not dispute — an oligo with no
binding site cannot prime — and it is **zero for all five assays**. For Charité RdRp the
3′-window clause alone accounts for 0.1112 of a 0.1116 rate.

![Figure 2](../research/figures/fig2_sensitivity_heatmap.png)

**Figure 2.** `likely_failing` for Charité RdRp across both thresholds, 2,759 usable
sequences. The blue box is the published setting. The rate ranges from 0.1116 to 0.1174
across choices that are all defensible, and neither threshold was calibrated against a
bench.

![Figure 3](../research/figures/fig3_three_prime_distribution.png)

**Figure 3.** Mismatches falling inside a primer's 3′ window, per assay, over 2,765
genomes: **395 assay–genome pairs carry exactly one and none carries two or more.** The
rule fires on one. The only direct in-vitro measurement of one reports about a single
cycle of Ct shift and no false negative; two produce a 9–12 cycle shift and a real one.
Charité RdRp's 11.16% is composed entirely of the left-hand bars.

![Figure 4](../research/figures/fig4_big_corpus_decomposition.png)

**Figure 4.** The re-scoring of section 4.5 over all 291,827 genomes, at the published settings, 281,456 to 289,945 usable sequences per assay. Clause 1 - an oligo with no binding site, the only clause a bench would not dispute - is zero for all five and so does not appear. Requiring **two** 3'-proximal mismatches instead of one divides Charité RdRp's rate by 990 and leaves CDC N1's essentially unchanged (0.0062 to 0.0050), because N1's rate is mismatch load and RdRp's is single 3'-proximal mismatches. The pair does occur at this corpus size - 52 times in 54,755 3'-proximal hits - which is what makes the surveillance framing a measured claim rather than an absence.

## 4.5 The same re-scoring over 291,827 genomes

§4.2 is scoped to the 2,765-genome corpus, and its load-bearing sentence — no genome
carries two mismatches in a 3′ window — is the kind of claim a hundredfold larger corpus
is apt to overturn. It does, and the overturning is the more useful result.

Re-scoring all 291,827 genomes at the published settings (3′ window 5, mismatch load 3):

| Assay | usable | published | **clustered 3′ (≥2)** | without clause 2 | clause 1 alone |
|---|---:|---:|---:|---:|---:|
| CDC N1 | 289,287 | 0.0062 | 0.0050 | 0.0050 | **0.0000** |
| CDC N2 | 288,108 | 0.0073 | 0.000014 | 0.000014 | **0.0000** |
| CDC N3 | 289,945 | 0.0023 | 0.000014 | 0.000014 | **0.0000** |
| Charité E | 285,298 | 0.0004 | 0.000070 | 0.000070 | **0.0000** |
| Charité RdRp | 281,456 | 0.1829 | **0.00018** | 0.000014 | **0.0000** |

Three things hold across the two corpora and one does not.

**Clause 1 still fires for nothing.** Zero for all five assays over 281,456 to 289,945
usable sequences each — the one clause with uncontested wet-lab backing contributes
nothing to the published rate at either scale. The single record that would have fired it,
a 2,859-base partial deposit, is now counted as `region_absent` rather than as a lost
oligo (§5), because a fragment that never contained the target region cannot be evidence
that the region was lost.

**The published rate is still almost entirely clause 2 at one mismatch.** Charité RdRp's
0.1829 is 51,465 genomes, and all 51,465 of them are 3′-proximal hits. Re-scoring clause 2
at two mismatches takes that rate to 0.00018, a factor of a thousand. §4.3's withdrawal
therefore does not rest on the small corpus; it is the same withdrawal, over a hundred times as
many genomes.

**Pairs exist, and they are rare in a way worth quantifying.** The column that was all
zeroes in §4.2 is not:

| Assay | 3′-proximal hits | of which a pair | share |
|---|---:|---:|---:|
| CDC N1 | 423 | 4 | 0.95% |
| CDC N2 | 2,097 | 0 | 0% |
| CDC N3 | 672 | 0 | 0% |
| Charité E | 98 | 0 | 0% |
| Charité RdRp | 51,465 | 48 | 0.093% |

52 pairs in 54,755 3′-proximal hits, and 52 in 1,434,094 assay–genome evaluations —
0.0036% of evaluations.

**An external number, and why it is not a validation.** BioFire's ongoing surveillance
reports 481 of 12,744,136 sequences, **0.0038%**, carrying "paired mismatches of concern"
[@biofire2023]. The two figures agree to one significant figure and they are not the same
quantity. Theirs is a 3′-half mismatch in *both* assays of a two-assay panel; ours is two
mismatches inside one primer's 3′ window. An agreement in magnitude between two different
definitions is a coincidence until someone computes one of them on the other's corpus,
and this paper does not. What does carry over is their design argument: a panel that calls
positive on either of two independent targets tolerates a mismatch that would be
diagnostic failure for a single-target test. All five assays studied here are
single-target laboratory-developed tests, so the mitigation that makes 0.0038% tolerable
for BioFire does not apply to them — which is a reason to watch the configuration, and
the same conclusion §4.3 reaches from the other direction.

This strengthens §4.3 rather than weakening it: the configuration Otwell et al. found does
cause false negatives is **reachable** — it was reached, by two assays, in six quarters —
and its base rate among the genomes `likely_failing` flags is about one in a thousand.
That is what "one mutation away" is worth as a number, and it is a figure only a corpus of
this size can supply. The RdRp pairs are concentrated in 2021-Q3 and 2021-Q4, 46 of the
48, the Delta quarters in which G15451A was spreading — which is the mechanism §4.3
describes, arriving at the 3′ terminus a second time.

**What changes is the arithmetic of the clustered column.** In §4.2 the clustered rate
equalled the rate with clause 2 dropped, because there were no pairs for clause 2 to
contribute. Here the two part for Charité RdRp — 0.00018 against 0.000014, 48 pairs
against 4 severe-load genomes — which is why the tool now carries `three_prime_pair` and
`severe_load` as separate counters rather than one `failing_clustered` total. A union
count cannot say which clause fired, and the published 0.1829 was taken apart only
because each clause is counted on its own.

# 5. Threats to validity

**The failure rule is uncalibrated** (§4). This is the dominant threat and it is the
paper's own section rather than a line here.

**The corpus is GenBank, not a clinical sample set.** Submission is uneven by country and
institution; the Kish correction addresses the arithmetic of that, not the selection.

**Short and partial records.** A record shorter than an oligo has no binding site by
construction, which is not the same fact as a site having been lost; those are `unknown`.
That sentence stood in an earlier draft of this section while the code behind it checked
only oligo length, and it was too weak by a factor of a hundred. LR877184 is 2,859 bases,
collected 2020-04-26 in Stockholm, with no unknown bases anywhere in it: a partial deposit
of a 29,900-base genome. It is longer than every 24-mer, so the length check passed it,
and it then aligned none of the fifteen oligos of the five assays — so it was reported, for
all five, as an assay with no binding site left. That is clause 1 of `likely_failing`, the
one clause with no judgement call in it, and in 291,827 genomes its only instance was this
record. The predicate is now whether the assay's target *region* is present: every oligo
absent **and** the record too short to be a whole genome (25,000 bases). The length half is
deliberate — a full-length genome that has genuinely lost an entire target region is still
reported as a lost oligo, because that is the finding the rule must not be able to hide.
A threat bullet that describes a guard is worth only as much as the guard, and this one
was describing a narrower check than it claimed.

**One alignment implementation.** The matcher is compared against Biopython's
`PairwiseAligner` on 300 oligo–genome pairs: 270 agree on the binding site, and every
full-length disagreement is a published primer's deliberate degenerate base that a
substitution matrix counts as a mismatch. That validates the aligner, not the rule.

**The sensitivity grid is the small corpus.** One pass over 291,827 genomes is about three
hours and the grid is twelve passes. The direction is reported; the magnitudes are the
2,765-genome subset's. The published setting alone was re-scored over all 291,827
genomes (§4.5) and the conclusion does not move, but the grid around it has not been.

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
- [@wetlab2025] Otwell, T., Knight, B., Coryell, M., Stone, J., Davis, P., Necciai, B., Carlson, P. and Sozhamannan, S. *Reality check: testing the in silico predictions of false negative results due to mutations in SARS-CoV-2 PCR assays using templates with mismatches in vitro.* Front. Cell. Infect. Microbiol. 15, 2025. <https://doi.org/10.3389/fcimb.2025.1524025>
- [@genes2021] *Deepening of In Silico Evaluation of SARS-CoV-2 Detection RT-qPCR Assays
  in the Context of New Variants.* Genes 12(4):565, 2021. Eleven primer/probe sets against
  84,305 genomes collected June 2020 to January 2021.
  <https://doi.org/10.3390/genes12040565>
- [@biofire2023] *Ongoing Surveillance of SARS-CoV-2 for In Silico Analysis to Monitor the
  Performance of SARS-CoV-2 Assays in BIOFIRE Panels.* Open Forum Infect. Dis. 10(Suppl 2),
  abstract 452, 2023. 12,744,136 sequences analysed as of 21 March 2023; 481 (0.0038%)
  carrying paired mismatches of concern, defined as a mismatch in the 3′ half of a primer
  in **both** assays of a two-assay panel.
  <https://doi.org/10.1093/ofid/ofad500.522>
- [@lefever2013] Lefever, S., Pattyn, F., Hellemans, J. and Vandesompele, J. *Single-Nucleotide Polymorphisms and Other Mismatches Reduce Performance of Quantitative PCR Assays.* Clin. Chem. 59(10):1470–1480, 2013. <https://doi.org/10.1373/clinchem.2013.203653>

> **Citation status.** @wetlab2025 **has now been read** — the article text, not a search
> summary — and doing so changed §4 substantially: it supplied the position-versus-count
> separation, the per-configuration effect sizes in 4.1, and the correction that the
> paper's conclusion is two-sided rather than the one-sided "in-silico overstates risk"
> an earlier draft attributed to it. About 100,000 of its 127,000 characters were read;
> the unread remainder is the latter part of the discussion, and no quantitative claim
> here is drawn from it.
>
> @lefever2013 was added as the actual origin of the 5-base window and the ≥4-mismatch
> blocking threshold. **Its abstract and conclusions were read; its full results section
> is paywalled.** Both quoted sentences are from the abstract and conclusions, which is
> enough for the two claims made of it, but the full text should be obtained before
> submission — in particular to check whether its results section reports anything about
> single mismatches at 1–2 bases from the 3′ terminus, which the abstract does not and
> which is the exact configuration §4.2 turns on.
>
> @biofire2023 **was read** — it is a conference abstract, so there is little to read,
> and both numbers quoted here (12,744,136 and 481 / 0.0038%) and its definition of a
> paired mismatch are taken from its own sentences. The definition is quoted because it
> differs from ours, which is the whole point of the comparison in §4.5.
>
> @genes2021 is **quoted from a search summary**: the publisher returns 403 and the PMC
> mirror is behind a browser check. The three facts used of it — eleven primer/probe
> sets, 84,305 genomes, June 2020 to January 2021, the China CDC N assay's inclusivity
> falling while the others held above 99% — must be confirmed against the article before
> submission, and the §1 contribution claim rests on a *negative* about it (that it
> reports no independence structure) which a summary cannot establish. That negative is
> the weakest sentence in §1 and is flagged as such here rather than in the text.
>
> The remaining four references are still quoted from search summaries and must be read
> before submission.
