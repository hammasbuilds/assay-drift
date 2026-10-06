# Changelog

Dates are when the work landed, not when a version was tagged. Nothing is on PyPI yet.

## Unreleased

### Two corrections to published claims

- **"both stayed in clinical use" is retracted.** It was the only external fact in the
  repository, it carried the whole argument, and it had no citation — and it is false as
  written: CDC withdrew its EUA request for the 2019-nCoV Real-Time RT-PCR Diagnostic
  Panel, which N1 and N2 belong to, effective after
  [2021-12-31](https://www.aacc.org/cln/articles/2021/september/cdc-planning-to-withdraw-request-for-eua-of-sars-cov-2-pcr-test),
  while the window here runs to 2026-Q3. The stated reason was that hundreds of other
  tests had been authorised, **not** that it had stopped detecting its target, which is
  the opposite of what a reader would infer — so that is now said, with the source.
- **"is the one that broke" is softened to "most likely to have failed".** It contradicted
  the repository's own Scope — *"It predicts, it does not test… `likely_failing` is a
  hypothesis, not a result"* — and it welded two quarters into one sentence: 90.8% is
  RdRp's 2026-Q3 exact match, which is its *healthiest* recent quarter, while the evidence
  behind "broke" is 44.5% likely-failing in 2021-Q4.

### Matcher

- **A coincidence is no longer reported as a binding site.** The mismatch cap was a flat 8;
  the catalogue's oligos are 18–26 bases, so that is 31% of a 26-mer and 44% of an 18-mer.
  Brute force over the reference puts the best *wrong* site for all fifteen oligos at
  27%–42% divergence, so the flat cap admitted a coincidence for ten of them. Masking all
  18 bases of CDC N2's reverse primer with `N` made `find` report a site 8,436 bases away
  with 7 mismatches and `usable` true. The cap is now a fifth of the oligo's length.
- **A drifted site is no longer missed.** Seeds were fixed 8-mers at every fourth offset,
  giving an 18-mer three seeds over bases 0–15 and leaving two bases in no seed at all, so
  two substitutions could break every seed. Enumerated over every pair of positions:
  **31.4% of two-mismatch sites in an 18-mer came back NOT_FOUND** — read as "this assay
  has no binding site" — and 60.8% of three-mismatch sites. Now 0.0% at every catalogue
  length, and 0.3%–2.2% for three. This changed no number in the committed run, because
  the mutations that circulate in these sites are single substitutions; it matters for
  `check`, where the primers are yours.

### Honesty about what the numbers rest on

- **Every rate now reports how independent it is.** `MIN_PERIOD_SAMPLES` counts sequences,
  and GenBank hands out consecutive accessions to one submission: 2024-Q1 holds 159
  sequences from 4 submitter groups, an effective n of **1.05**, and passes a 25-sequence
  floor. Every period row carries `clusters` and `effective_n`; the median across the
  quarters that pass the floor is about 4.
- **The corpus is reproducible.** `data/accessions.tsv` commits the identity of all 2,765
  records with a hash, and `scripts/fetch.py --from-manifest` fetches exactly those. The
  previous instruction — "`scripts/fetch.py` rebuilds it" — rebuilt *a* corpus: the
  deposit-year windows for 2025 and 2026 are still filling, so 28% of the records came
  back different.
- **`results/drift.json` records its own provenance** in a `_run` block: when, over which
  corpus, which manifest, and the command to rebuild it.
- **The stratification is described accurately.** "Stratified by year" is the *deposit*
  year, and the newest-deposited records in a year were mostly collected in its final
  quarter — so 2021-Q4/2022-Q4/2023-Q4/2025-Q4 hold 330/421/423/327 while 2022-Q1,
  2024-Q2 and 2024-Q3 hold nothing. The series is seven annual clumps, not a quarterly
  time series, and the README now says so and points at the year view.

### Fixed

- `check` reported a four-base FASTA record as `likely_failing` with every oligo "lost".
  A record too short to hold an oligo is not evidence that the assay fails; it is
  `unknown`.
- `--json` dropped `Hit.usable` — whose docstring is *"whether this measurement can go
  into a rate"* — because it is a property. A consumer building a rate from that file
  counted exactly the rows it exists to exclude. It now also carries
  `three_prime_mismatches` as a count and `indels`.

### Packaging

- PEP 639 licence metadata (`License-Expression: MIT`), replacing the deprecated table
  form plus classifier.
- README links are absolute, so they resolve on PyPI rather than 404.
- A release workflow using PyPI Trusted Publishing, matching the sibling projects.
