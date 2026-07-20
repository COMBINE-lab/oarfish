# Degradation-aware coverage evaluation: 2026-07-19

> Follow-up: the intact/degraded kernel evaluated here has been superseded by
> the validated [three-component competing-risk model](competing-risk-degradation-evaluation-2026-07-20.md),
> which separates length-invariant technical truncation from degradation.

## Model

`--coverage-model degradation --seq-tech ont-drna` models the global
sample-level process that generates 3'-anchored fragments before interpreting
endpoint patterns as transcript evidence. For transcript length `L` and
sequenced fraction `u`, the degraded component is a truncated exponential
survival distribution with rate `hazard_per_kb * L`. An intact point component
accounts for alignments in the first 5' gap bin.

The hazard is selected on a small regularized grid and the intact fraction is
estimated continuously by MAP-EM using
only single-candidate, 3'-anchored reads. Every read is scored with parameters
trained without its fold. Candidate endpoint likelihood is divided by this
global nuisance likelihood and renormalized.

The unified `auto` mode learns three levels of weighting. A sample-level null
sets degradation correction to zero when a minimum-hazard technical process
fits as well as degradation. Otherwise, each read is weighted by its posterior
probability of degradation. Finally, normalized mixture variance
`4 * intact_fraction * (1 - intact_fraction)` controls whether both components
are represented well enough to identify the correction. The same learned
quantity attenuates endpoint reliability before the existing agreement gate
and Bayes-factor cap.

This is deliberately an ONT direct-RNA kernel. The command requires an explicit
`--seq-tech ont-drna`, including for BAM input, so it cannot silently apply the
wrong directional model to cDNA or PacBio data.

## Complete degradation trajectories

All 15 direct-RNA samples with public FASTQ in
[ENA PRJEB53210](https://www.ebi.ac.uk/ena/browser/view/PRJEB53210) were evaluated
as fixed 50,000-read prefixes. Parameters were frozen before this expansion.

| Trajectory | RIN sequence | Learned hazard sequence (/kb) | RIN/hazard Pearson |
|---|---|---|---:|
| TS10 | 9.8, 9.6, 9.3, 8.4, 7.7 | 0.20, 0.20, 0.20, 0.36, 0.40 | -0.969 |
| TS11 | 9.7, 9.6, 9.3, 8.9, 8.8, 8.7 | 0.10, 0.10, 0.40, 0.70, 0.40, 0.40 | -0.767 |
| TS12 | 9.9, 8.2, 7.3, 7.2 | 0.20, 0.40, 0.70, 0.58 | -0.948 |
| pooled | 15 samples | 0.356 mean | -0.764 |

The hazard tracks lower RIN strongly in every trajectory. It is not perfectly
monotone—most visibly around TS11 RIN 8.9—which is expected from 50k-read
subsets and from RIN measuring total RNA rather than only sequenceable poly(A)
molecules. The fitted intact fraction is less informative: it is 0.5 for 14/15
samples and 0.25 for TS11 RIN 8.9. The hazard, rather than mixture weight, is
currently the identifiable degradation diagnostic.

### High-to-low RIN robustness

| Trajectory | Model | Pearson | Jensen-Shannon | L1 distance |
|---|---|---:|---:|---:|
| TS10 9.8→7.7 | none | **0.9662** | **0.1415** | **0.6463** |
|  | adaptive | 0.9542 | 0.1594 | 0.6825 |
|  | degradation | 0.9605 | 0.1585 | 0.6780 |
|  | **auto** | 0.9655 | 0.1515 | 0.6651 |
| TS12 9.9→7.2 | none | **0.9468** | **0.1298** | **0.6313** |
|  | adaptive | 0.9196 | 0.1538 | 0.6824 |
|  | degradation | 0.9266 | 0.1526 | 0.6795 |
|  | **auto** | **0.9479** | 0.1435 | 0.6541 |
| TS12 9.9→7.3 | none | **0.9378** | **0.1283** | **0.6155** |
|  | adaptive | 0.9118 | 0.1537 | 0.6698 |
|  | degradation | 0.9200 | 0.1513 | 0.6616 |
|  | **auto** | **0.9392** | 0.1423 | 0.6365 |

Automatic degradation-aware correction improves over adaptive on every
trajectory and metric shown. On both TS12 comparisons its Pearson correlation
slightly exceeds no coverage, while Jensen-Shannon and L1 remain between no
coverage and adaptive. Since
degradation genuinely changes observed abundance, this is a robustness test,
not absolute transcript truth; nevertheless, the remaining gap is a reason to
keep the model opt-in.

## Frozen truth benchmarks

| Dataset | Model | CCC | MARD | RMSE |
|---|---|---:|---:|---:|
| large synthetic dRNA | none | 0.956952 | 0.169178 | 89.4590 |
|  | adaptive | **0.996935** | 0.143259 | **24.0550** |
|  | degradation | 0.991278 | **0.143207** | 40.3891 |
|  | **auto** | **0.996935** | 0.143259 | **24.0550** |
| SIRV E0 dRNA 50k | none | 0.639192 | 0.166302 | 0.529723 |
|  | adaptive | **0.671348** | **0.154340** | **0.493311** |
|  | degradation | 0.641909 | 0.165376 | 0.526606 |
|  | **auto** | **0.671323** | **0.154345** | **0.493338** |
| source-labelled dRNA | none | 0.981515 | 0.0001168 | 0.012606 |
|  | degradation | 0.982339 | 0.0001127 | 0.012317 |
|  | adaptive | **0.983387** | **0.0001062** | **0.011939** |
|  | **auto** | 0.983374 | 0.0001063 | 0.011944 |

The final automatic model retains adaptive's SIRV and large-synthetic accuracy:
both select the minimum hazard and correction weight zero. On the nearly
full-length source-labelled set, `auto` is numerically indistinguishable from
adaptive because intact fraction is 0.998 and normalized mixture variance is
near zero. Thus the degradation robustness no longer costs accuracy on any of
the frozen truth benchmarks.

Across the 15 degradation samples, coverage construction averages 0.157 s and
EM averages 1.235 s. This is about 0.12 s more coverage work than adaptive on
the earlier six-sample panel. The large 1.36-million-read synthetic run spends
4.49 s on coverage and 16.95 s in EM.

## Decision

- Prefer experimental `auto` over asking users to choose `adaptive` versus
  `degradation`. It preserves adaptive behavior when degradation is unsupported
  and continuously increases degradation correction when it is identifiable.
- Keep the overall default at `none` until the full multi-sample LongBench truth
  evaluation is complete; automatic model selection solves the interface and
  major accuracy tradeoff, but the external evidence base is still modest.
- The next model iteration should replace the coarse intact-fraction grid with
  a continuous or beta-binomial estimate and allow a third component for
  technical sequencing termination. Those components need sample-level outer
  validation before changing defaults.
- Technology separation is warranted. ONT cDNA and PacBio should receive their
  own observation kernels rather than reuse this 3'-anchored survival process.

Artifacts are under
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/degradation` and
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/degradation-model-eval`.
