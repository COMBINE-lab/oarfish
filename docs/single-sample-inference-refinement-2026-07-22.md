# Single-sample inference refinement (2026-07-22)

This round tested three annotation-free improvements: learned alignment-score
temperature, read-level alignment censoring, and Pareto-dominance candidate
pruning. Accuracy, wall time, and peak RSS were measured against the unchanged
`dev` baseline. The primary panel contains synthetic direct RNA and 24 LongBench
50k technology/sample combinations plus three PacBio 250k confirmations.

## Alignment-score calibration: rejected learner

A fixed sweep of score denominators `D={1,2,3,5,8,12}` established that score
strength matters. Sharper scores generally improved real ONT and PacBio rank and
error metrics, but the synthetic data preferred a flatter score for Pearson and
showed no stable rank optimum.

The first single-sample learner used two-fold held-out predictive likelihood to
select `D`, with abundance fitted only on the complementary reads. It selected
the flattest candidate (`D=12`) for every sample. This objective rewards placing
mass across the observed candidate set rather than identifying the generating
transcript. It regressed real-data rank/MARD and added 3.3--3.8 seconds to 50k
human runs and 34 seconds to the synthetic run. The implementation was removed.
The fixed `--score-prob-denom` remains available, but its default is unchanged.

## Adaptive read-level censoring

Each alignment now records terminal query clipping. For candidate transcript
`t`, unexplained clipping is

`min(left_clip, transcript_left_gap) + min(right_clip, transcript_right_gap)`.

Clipping beyond a transcript end is therefore neutral, whereas clipped query
sequence that could have aligned to remaining transcript sequence is evidence
against that candidate. The exponential censoring scale is learned from
single-candidate reads with 100 pseudo-observations at 50 nt, clamped to
25--500 nt. Sample support shrinks its exponent toward zero, and its per-read
Bayes factor is capped at four. This is selected by `--censoring-model adaptive`.

## Pareto-dominance pruning

Candidate B is pruned only if candidate A is no worse in alignment likelihood,
coverage/censoring likelihood, and aligned transcript span, and A's joint
likelihood is at least twice B's. At least one candidate is retained. The
implementation reuses a scratch mask across reads and zeros the losing
likelihood without reallocating the packed equivalence-class store. It is
selected by `--candidate-pruning dominance --dominance-bayes-factor 2`.

## Full-panel result

Dominance alone improved Spearman and MARD on all 28 panel samples, Pearson on
26/28, and every metric on all 8 cDNA and all 11 PacBio samples. Mean Pearson
changes were +0.00109 cDNA, +0.00075 dRNA, and +0.00729 PacBio. It also reduced
mean EM time by 0.076 seconds.

Combining censoring and dominance was stronger:

| Technology | Samples | Pearson wins | Spearman wins | Mean Pearson delta | Mean Spearman delta | Mean MARD delta |
|---|---:|---:|---:|---:|---:|---:|
| ONT cDNA | 8 | 8 | 8 | +0.00408 | +0.00250 | -0.000183 |
| ONT dRNA | 9 | 8 | 9 | +0.00190 | +0.00264 | -0.001037 |
| PacBio | 11 | 11 | 11 | +0.00837 | +0.00279 | -0.000284 |
| All | 28 | 27 | 28 | +0.00506 | +0.00266 | -0.000497 |

The sole Pearson regression was synthetic dRNA (-0.00177), accompanied by a
Spearman gain of +0.00374 and MARD improvement of -0.00697. One PacBio sample
had a negligible MARD regression (+0.000013) while its Pearson and Spearman
improved. Median peak-RSS change was only +56 KiB; mean wall time fell by 0.079
seconds because inference converged faster.

## Independent outer validation

The frozen combination was then tested on public SIRV and separately generated
100k ONT/PacBio simulations:

| Dataset | Pearson delta | Spearman delta | MARD delta |
|---|---:|---:|---:|
| SIRV E0 dRNA 50k | 0 | 0 | -0.04983 |
| SIRV E2 dRNA 50k | +0.03102 | +0.11898 | -0.06403 |
| SIRV E0 cDNA 50k | 0 | 0 | -0.01247 |
| Independent ONT 100k | +0.01103 | +0.06552 | -0.00524 |
| Independent PacBio 100k | +0.01197 | +0.05510 | -0.00549 |

The combined model is therefore part of the new automatic transcriptome inference stack through
`--censoring-model auto` and `--candidate-pruning auto` when transcriptome input
is combined with `--coverage-model auto`. Both automatic modes abstain for
other coverage selections and genome-projection input, which were not part of this validation.
Users can reproduce the previous behavior with `--censoring-model none
--candidate-pruning none`.

Benchmark outputs are under `oarfish-evaluation-data` in directories beginning
`alignment-temperature-screen`, `dominance-full`, `censoring-dominance-full`,
and `censoring-dominance-outer` dated 20260722.
