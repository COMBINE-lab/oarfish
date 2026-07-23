# Single-sample accuracy refinement (2026-07-23)

This round evaluated two annotation-free refinements in sequence against the
adaptive rank-preserving model. Both accuracy and wall-time/peak-RSS costs were
gated. Only the second refinement was retained.

## Equivalence-component mass conservation: rejected

The first experiment built transcript components from candidates co-occurring
on ambiguous reads. After abundance blending, each component was rescaled to
preserve its pre-blend corrected abundance. A broad implementation incorrectly
rescaled singleton components and visibly lost rank accuracy. A restricted
version conserved only multi-transcript components representing at least 0.25%
of the library.

The restricted model improved synthetic Pearson by 0.00015, CCC by 0.00016,
and RMSE by 0.91 with essentially unchanged Spearman. It also improved H69 cDNA
RMSE by 9.4 and H69 dRNA RMSE by 5.1. However, it reduced H69 cDNA Pearson by
0.00141 and CCC by 0.00040, and worsened both PacBio test cases. Building and
rescaling the components added approximately 0.5--1.1 seconds and 9--17 MiB RSS
on the human screens. The accuracy/cost trade was not acceptable, so the
implementation and CLI variants were removed.

Artifacts:

- `oarfish-evaluation-data/rank-blend-family-screen-20260723`
- `oarfish-evaluation-data/rank-blend-family-multitarget-screen-20260723`

## Agreement-calibrated alignment likelihood: retained

The retained model sharpens a read's alignment-score likelihood only when the
score-winning alignment also has both the least unexplained terminal clipping
and the longest aligned transcript span. The exponent is strongest for score
ties and is regularized toward one as candidate multiplicity grows:

`alpha = 1 + 4 * (second_score_probability / top_score_probability)^2
             / sqrt(candidate_count - 1)`.

Probabilities are updated in place as `p <- p^alpha`. This combines alignment
score, read-level censoring evidence, aligned span, and candidate multiplicity
without an annotation or truth labels. It avoids the failure of the earlier
held-out-likelihood temperature learner, which selected the flattest candidate
for every sample.

On the 28-case primary panel, relative to adaptive rank preservation alone:

| Metric | Mean delta | Wins |
|---|---:|---:|
| Pearson | +0.000018 | 17/28 |
| Spearman | +0.000252 | 22/28 |
| CCC | +0.000052 | 21/28 |
| RMSE | -0.795 | 21/28 |
| MARD | -0.000021 | 23/28 |

The six Spearman regressions were all small (largest absolute regression
0.00050). PacBio improved most consistently: mean Spearman +0.000316 (10/11),
Pearson +0.000061, CCC +0.000072, RMSE -0.432, and MARD -0.000039.

Outer validation also passed:

| Dataset | Pearson delta | Spearman delta | CCC delta | MARD delta |
|---|---:|---:|---:|---:|
| SIRV E2 dRNA | +0.000414 | 0 | +0.000380 | -0.000787 |
| Independent ONT | +0.000166 | +0.000768 | +0.000171 | -0.000058 |
| Independent PacBio | -0.000081 | +0.000430 | -0.000075 | -0.000065 |

SIRV E0 rank correlations were unchanged. A stronger exponent coefficient of
eight improved PacBio further but reduced synthetic and dRNA Spearman; it was
rejected in favor of the conservative coefficient of four.

Two same-binary paired repeats on synthetic dRNA, H69 cDNA, and H69 PacBio
measured a mean 3.1 ms increase in the coverage/calibration phase. Mean peak-RSS
changed by -27 KiB (measurement noise), confirming that the in-place scan adds
no material memory. Whole-run and EM timing varied in both directions between
repeats and showed no reproducible overhead.

`--alignment-calibration auto` is the default. It activates only for automatic
transcriptome coverage inference and abstains for genome projection, explicit
non-auto coverage modes, and an explicitly supplied `--score-prob-denom`.
`agreement` forces the model and `none` disables it.

Artifacts:

- `oarfish-evaluation-data/alignment-agreement-rank-full-20260723`
- `oarfish-evaluation-data/alignment-agreement-outer-20260723`
- `oarfish-evaluation-data/alignment-agreement-cost-{none,on}-20260723`
- `oarfish-evaluation-data/alignment-agreement-strength8-screen-20260723`
