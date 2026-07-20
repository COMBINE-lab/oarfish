# Competing-risk degradation evaluation: 2026-07-20

This stage tested two frozen alternatives to the original ONT direct-RNA
degradation kernel: internally learned scalar correction weights and an
observation model that separates degradation from technical truncation.

## Rejected: alignment-predictive scalar calibration

The original intact/degraded model was evaluated at correction weights
`0, 0.1, 0.25, 0.5, 0.75, 1`. Cross-fitted ambiguous reads were scored by the
agreement of corrected endpoint probabilities with their independent alignment
score probabilities. A conservative prior favored zero and likelihood evidence
was capped at 50 effective reads.

This successfully detected LongBench overcorrection. Learned weights were
0.03--0.07 in the four active 50,000-read libraries, and HCC827 at 250,000
reads improved from CCC 0.558229 to 0.562535, essentially equal to adaptive
(0.562512). It failed the degradation-trajectory gate, however:

| Trajectory | Old auto Pearson | Predictive-weight Pearson |
|---|---:|---:|
| TS10 9.8→7.7 | 0.96549 | 0.95439 |
| TS12 9.9→7.2 | 0.94790 | 0.91991 |
| TS12 9.9→7.3 | 0.93919 | 0.91205 |

The selector assigned the degradation series the same small weights as
LongBench, collapsing robustness almost to adaptive. Alignment agreement is
therefore not an adequate single-sample proxy for quantification benefit. The
selector was removed rather than retained as unused runtime overhead.

## Accepted: three-component competing-risk model

The accepted model represents a 3'-anchored read as one of:

1. intact, concentrated in the full-length 5' endpoint bin;
2. technical truncation, a length-invariant Beta(2,1) retained fraction; or
3. degradation, a truncated exponential whose rate scales with physical
   transcript length.

Mixture weights are estimated continuously by cross-fitted MAP-EM with a weak
Dirichlet(5,2,2) prior. Hazard is selected on the frozen regularized grid. The
correction divides only the degradation likelihood ratio relative to the
intact/technical mixture, rather than removing the entire observed endpoint
profile. Its strength remains learned from degradation-versus-null evidence,
the degraded mixture fraction, and each read's degradation posterior. Existing
adaptive reliability and Bayes-factor safeguards remain in force.

Across LongBench dRNA, the model estimates 7.6--20.1% degraded molecules and
20.0--26.5% technical truncation at 50,000 reads. These components are separated
by their different transcript-length dependence rather than endpoint shape
alone.

## Results

### Degradation robustness

| Trajectory | None Pearson | Old auto Pearson | Competing-risk Pearson | Competing-risk JSD | Competing-risk L1 |
|---|---:|---:|---:|---:|---:|
| TS10 9.8→7.7 | **0.96623** | 0.96549 | 0.96500 | 0.15519 | 0.66866 |
| TS12 9.9→7.2 | 0.94682 | 0.94790 | **0.94876** | 0.14290 | 0.64831 |
| TS12 9.9→7.3 | 0.93781 | 0.93919 | **0.94013** | 0.14257 | 0.63269 |

The new model preserves the original auto robustness and slightly improves both
TS12 comparisons.

### LongBench matched-Illumina comparator

At 50,000 reads, competing-risk auto is within 0.00027 CCC of adaptive or
better in all eight cell lines, while every cell line remains above `none`.
It improves over the old auto in each library where old degradation correction
was active. At 250,000 HCC827 reads:

| Model | CCC | Pearson | RMSE |
|---|---:|---:|---:|
| adaptive | 0.562512 | 0.739364 | 6869.94 |
| old auto | 0.558229 | 0.735039 | 6919.07 |
| **competing-risk auto** | **0.562581** | **0.739437** | **6869.19** |

### Frozen truth

| Dataset | CCC | MARD | RMSE |
|---|---:|---:|---:|
| large synthetic dRNA | 0.996926 | 0.143267 | 24.0920 |
| SIRV E0 dRNA | 0.671364 | 0.154332 | 0.493293 |
| source-labelled dRNA | 0.983387 | 0.0001062 | 0.011939 |

These are effectively equal to the frozen adaptive optima. The model assigns
only 0.69% degraded fraction to SIRV and 0.22% to the nearly full-length
source-labelled set, making correction negligible without a manually selected
mode.

## Runtime and decision

Caching component probabilities during mixture EM reduced coverage construction
on the 1.36-million-read synthetic benchmark from 11.3 s during development to
2.66 s; EM took 17.45 s. This is also below the previous two-component auto
measurement of 4.49 s. The cached and uncached implementations produce the same
fitted parameters and estimates.

Retain the three-component constant-hazard model. Bounded beta and piecewise
challengers were subsequently tested; neither displaced the incumbent. See the
[challenger study](degradation-kernel-challengers-2026-07-20.md). Gene-random
effects and joint abundance/degradation inference remain unjustified.

Artifacts are under `oarfish-evaluation-data/degradation/three-component*` and
`oarfish-evaluation-data/longbench/three-component*`.
