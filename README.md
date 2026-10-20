# assay-drift

A tool that checks whether published PCR diagnostic assays still match the virus that is circulating. We built it
for diagnostic labs, public-health analysts and assay designers: a PCR test fails silently when the target mutates
under a primer, because a failed amplification looks exactly like a negative sample. assay-drift aligns every oligo
of an assay against dated GenBank genomes, separates "no longer an exact match" from "would plausibly no longer
amplify", and names the mutation responsible. The committed run covers 2,765 SARS-CoV-2 genomes from 2020 to 2026
against five assays clinical laboratories ran. Standard-library Python, zero dependencies.

**Headline:** two of five published assays lost their exact match entirely (99.3% to 0.0%) and neither mutation is
at a position that stops the reaction; Charité RdRp never lost its exact match, yet in 2021-Q4 44.5% of sequences
carried a mutation one base from its primer's 3′ end. Where the mutation lands matters more than how much of it
there is.

## What it does

- **Downloads dated genomes from GenBank** through NCBI E-utilities, rate-limited and cached, stratified by deposit
  year. A record with no usable collection date is dropped and counted (35 of 2,800), never defaulted to today.
  `--from-manifest` rebuilds the exact corpus from the committed accession list.
- **Aligns each oligo** with IUPAC codes read as sets on both sides: primer `R` matches `A` and `G`; a subject `R`
  where the primer wants `A` is ambiguous, not a mismatch. A site must be within a fifth of the oligo's length in
  mismatches to count as the binding site at all, and the seed search is sized so two substitutions anywhere cannot
  hide a site.
- **Treats an unknown base (`N`) under an oligo as missing data**, excluded from the rates and reported, so a change
  in sequencing quality does not show up as viral drift.
- **Scores each sequence** as `perfect`, `drifted` (mismatched but expected to amplify), `likely_failing` (an oligo
  lost, a mismatch in a primer's last five bases, or three or more mismatches across the assay) or `unknown`. The 3′
  window counts only on primers: a hydrolysis probe is blocked at its 3′ end and never extended.
- **Groups by collection quarter or year** and publishes, beside every rate, the number of submitter groups and the
  Kish effective sample size, because sequences from one deposit batch are not independent.
- **Names the mutations** behind each change from raw sequence, with genome coordinate, position in the oligo and
  the places they were collected.
- **Checks your own primers against your own FASTA** with `assay-drift check`, offline.

`likely_failing` is a prediction from sequence complementarity: no PCR was run, and there is no melting-temperature
model (the design side, with nearest-neighbour thermodynamics, is
[primer-designer](https://github.com/hammasbuilds/primer-designer)).

## Screenshots

A small browser front end drives the package for demonstrations. It is a local demo and is **not part of
the package** — the package has no web surface and no runtime dependencies. Every result on these screens
is a live call into the library.

![Overview](https://raw.githubusercontent.com/hammasbuilds/assay-drift/main/docs/screenshots/00-overview.png)
*The front end as it opens: every feature is a tab, and the page states that nothing on it is precomputed.*

![Break An Assay](https://raw.githubusercontent.com/hammasbuilds/assay-drift/main/docs/screenshots/01-break-an-assay.png)
*A substitution placed on the 3′ base of primer N1-F: mismatches 0 → 1, and **blocks extension 0 → 1**. The same change in the middle of the oligo would usually be tolerated.*

![Score Reference](https://raw.githubusercontent.com/hammasbuilds/assay-drift/main/docs/screenshots/02-score-reference.png)
*An assay scored against the reference genome, oligo by oligo.*

![Six Years Drift](https://raw.githubusercontent.com/hammasbuilds/assay-drift/main/docs/screenshots/03-six-years-drift.png)
*Exact match and likely-failing share per quarter, 2020 to 2026, for all five assays.*

![Your Own Primers](https://raw.githubusercontent.com/hammasbuilds/assay-drift/main/docs/screenshots/04-your-own-primers.png)
*Your own primers and sequences scored — typed, pasted, or dropped in as files.*


## Sample input / output

`assay-drift check` on the shipped example: the CDC N2 primers and two genomes, the 2019 reference and the same
genome carrying C29215T.

```
$ assay-drift check --primers examples/cdc_n2.tsv --sequences examples/two_genomes.fasta
sequence                verdict                 N2-F        N2-R        N2-P
NC_045512.2             perfect                  0mm         0mm         0mm
mutant_C29215T          likely_failing           0mm      1mm/3'         0mm

perfect 1, drifted 0, likely_failing 1, unknown 0  (of 2)
```

`python demo.py` first checks every published oligo against the 2019 reference (a mistyped base would read as drift
that is not there), then reprints the finding from `results/drift.json` (trimmed):

```
  oligo     role      len  position  mismatch   note
  N1-F      forward    20     28286         0
  N1-R      reverse    24     28334         0
  N1-P      probe      24     28308         0
  N2-F      forward    20     29163         0
  N2-R      reverse    18     29212         0
  ...
  RdRp-R    reverse    26     15504         1   known: S vs T, 14 bases from the 3' end
  RdRp-P    probe      25     15469         0

  14 of 15 match perfectly.

  assay                            exact match        likely failing
  CDC N1            2020-Q1  99.3%  ->  2026-Q3   0.0%       10.9%   worst 10.9% (2026-Q3)
  CDC N2            2020-Q1  98.7%  ->  2026-Q3  87.3%        4.2%   worst 34.5% (2026-Q2)
  CDC N3 (retired)  2020-Q1  98.7%  ->  2026-Q3  86.7%        0.5%   worst 3.4% (2021-Q1)
  Charite E         2020-Q1  99.3%  ->  2026-Q3   0.0%        0.0%   worst 1.8% (2026-Q1)
  Charite RdRp *    2020-Q1  99.3%  ->  2026-Q3  90.8%        9.2%   worst 44.5% (2021-Q4)
```

## Results

Every number below comes from [`results/drift.json`](results/drift.json) and
[`results/mutations.json`](results/mutations.json), and `tests/test_published_numbers.py` checks each one against
them. Assays: CDC-006-00019 rev.06 (N1, N2, N3) and Corman et al. 2020, *Euro Surveill* 25(3) (Charité E, RdRp).

&#128202; **The same figures as a dataset:
[PCR Assay Drift — SARS-CoV-2](https://www.kaggle.com/datasets/muhammadhammas13/pcr-assay-drift-sars-cov-2)**
on Kaggle — one row per assay per quarter as CSV, the accession list for all 2,765 records,
and both result files. CC BY-SA 4.0; please cite it if you use it.

**Exact match and likely failing, by collection quarter** (first and last quarters with at least 25 usable
sequences):

| Assay | Exact match 2020-Q1 | Exact match 2026-Q3 | Change | Likely failing, worst quarter |
|---|---:|---:|---:|---|
| **CDC N1** | 99.3% | **0.0%** | **−99.3** | 10.9% (2026-Q3) |
| **Charité E** | 99.3% | **0.0%** | **−99.3** | 1.8% (2026-Q1) |
| **Charité RdRp** | 99.3% | 90.8% | −8.6 | **44.5% (2021-Q4)** |
| CDC N3 *(retired)* | 98.7% | 86.7% | −12.0 | 3.4% (2021-Q1) |
| **CDC N2** | 98.7% | 87.3% | −11.4 | **34.5% (2026-Q2)** |

**How independent those quarters are.** Sequences from one submission are not independent observations, so each
quarter is grouped by (country, accession block) and given a Kish effective sample size:

| quarter | sequences | submitter groups | effective n |
|---|---:|---:|---:|
| 2020-Q1, the baseline of every row above | 153 | 11 | **2.4** |
| 2021-Q4, Charité RdRp's worst quarter | 330 | 20 | 5.2 |
| 2024-Q1 | 159 | **4** | **1.1** |
| 2026-Q3, the endpoint of every row above | 196 | 33 | 4.4 |

The median across quarters that pass the floor is about 4. The 0.0% endpoints are large enough to survive that;
quarter-to-quarter shape and worst-quarter figures are not, and no confidence interval is quoted. Downloads are
stratified by deposit year, so collection quarters are uneven (2022-Q1, 2024-Q2 and 2024-Q3 hold no sequences);
the year view below is the steadier one.

### The same data by year

`scripts/analyse.py --granularity year` regroups the same 2,765 sequences:

| Assay | Exact match 2020 | Exact match 2026 | Likely failing 2026 | Likely failing, worst year |
|---|---:|---:|---:|---|
| **CDC N1** | 99.0% | **0.0%** | 7.0% | 7.0% (2026) |
| **Charité E** | 98.8% | **0.0%** | 0.9% | 0.9% (2026) |
| **Charité RdRp** | 98.3% | 87.0% | 13.0% | **38.3% (2021)** |
| CDC N3 *(retired)* | 99.3% | 90.3% | 0.9% | 2.0% (2023) |
| **CDC N2** | 98.0% | 68.8% | **15.8%** | 15.8% (2026) |

### The mutation behind each change

| Assay | Position in oligo | Genome coordinate | Change | Share of all sequences |
|---|---|---|---|---:|
| CDC N1 probe | base 3 of 24 (21 from the 3′ end) | **C28311T** | C→T | 76.5% |
| Charité E forward | base 2 of 26 (24 from the 3′ end) | **C26270T** | C→T | 76.7% |
| Charité RdRp forward | **1 base from the 3′ end** | **G15451A** | G→A | 11.1% |
| CDC N2 reverse | **2 bases from the 3′ end** | **C29215T** | C→T | 2.3% |

The first two are Omicron substitutions, mid-oligo, which weaken binding without blocking extension: that is why
CDC N1 and Charité E lost their exact match and are still not predicted to fail (N1's residual likely-failing rate
comes from sequences with three or more mismatches across the assay). The third is Delta's NSP12 G671S, at the 3′
terminus where polymerase extends: Charité RdRp's likely-failing rate hit 44.5% in 2021-Q4, fell to 3.1% in 2022-Q4
when Omicron displaced Delta, and runs at 9-24% across the 2025-2026 quarters. The share column is pooled over all
years, so a recent mutation looks small there and large in its own quarter.

### CDC N2 in 2026

The fourth mutation sits two bases from the 3′ end of the N2 reverse primer. Its share of sequences where the site
is readable:

| Quarter | carrying C29215T | n |
|---|---:|---:|
| 2020-Q1 … 2024-Q4 | 0.0% | 153–420 per quarter |
| 2025-Q3 | 1.6% | 62 |
| 2025-Q4 | 3.8% | 316 |
| 2026-Q1 | 28.7% | 115 |
| 2026-Q2 | 34.5% | 29 |
| 2026-Q3 | 4.2% | 189 |

It is absent from 2,000+ sequences before 2025 and the base is unambiguous in every sequence counted. Within 2026 the
mutation is 51 of 333 readable sequences (15.3%). The per-quarter swing is mostly where the sequencing was done: the
hits arrive in single-submitter batches, and the places stored for it, pooled over all years, are
`USA: Wisconsin 24/74, USA: California 10/186, USA: Michigan 7/27, USA: Washington 4/163, USA 3/170`. 2026-Q2 rests
on 29 sequences. GenBank is not a random sample of infections, so these rates describe deposited sequences.

CDC withdrew its EUA request for the 2019-nCoV Real-Time RT-PCR Diagnostic Panel (N1, N2) effective after
[2021-12-31](https://www.aacc.org/cln/articles/2021/september/cdc-planning-to-withdraw-request-for-eua-of-sars-cov-2-pcr-test),
because other tests had been authorised, not because it stopped detecting its target.

## Tech stack

- Python 3.10+, standard library only: `urllib` for NCBI E-utilities, `re`, `json`, `gzip`, `argparse`
- The 2019 reference genome (NC_045512.2) vendored as 9 KB of gzipped JSON with a pinned sha256, so tests never
  touch the network
- pytest and ruff (pinned 0.16) for development; CI on Ubuntu and Windows, Python 3.10-3.13

## Run it

```bash
git clone https://github.com/hammasbuilds/assay-drift
cd assay-drift
python demo.py          # the setup check and the finding, no network, ~2 s
pytest -q               # 224 tests, no network, no install step
```

To rebuild the data from scratch:

```bash
python scripts/fetch.py --per-year 400    # ~10 min, polite to NCBI, cached for a month
python scripts/analyse.py                 # writes results/drift.json and results/mutations.json
```

Or with `make`: `make demo`, `make test`, `make lint`, `make fetch`, `make analyse`.

Your own primers against your own sequences:

```bash
pip install -e .        # only needed for the CLI; demo and tests need no install
assay-drift check --primers my_assay.tsv --sequences genomes.fasta --json out.json
```

`my_assay.tsv` is one oligo per line, `name<TAB>role<TAB>sequence`, role being `forward`, `reverse` or `probe`; `#`
comments and blank lines are skipped. `--json` writes per-oligo positions, mismatch counts and unknown-base counts.

## Licence

MIT, see [LICENSE](LICENSE). The assay sequences are published in the cited sources; the genomes are GenBank's
(NCBI).
