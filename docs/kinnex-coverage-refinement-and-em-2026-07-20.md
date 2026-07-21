# Kinnex coverage refinement and EM optimizer study (2026-07-20)

## Questions

1. Can a sample-adaptive transcript correction recover the Figure 2D loss
   without discarding the Figure 2E benefit?
2. Do ordinary EM, SQUAREM, and DAAREM change accuracy, convergence, or cost?

All comparisons use the six author-provided one-million-read Kinnex SIRV BAMs.
Matched Day 0/Day 5 truth was used only for evaluation, never for fitting.

## Rejected stability-displacement refinement

We tested a transcript-specific stability gate. It compared the coverage-aware
estimate with a coverage-free warm start and increasingly shrank unusually
large transcript displacements. Its scale was the sample's median relative
displacement, so no truth or annotation was needed.

With SQUAREM, the candidate produced Figure 2D Pearson 0.7993 and Figure 2E
Pearson 0.6849. This improved D relative to the corrected abundance blend
(0.7713) but reduced E (0.6974). Multipliers from 0.5 to 4 produced essentially
the same point, indicating saturation rather than useful sample adaptation.
The candidate therefore did not dominate the simpler model and was removed
from production code after evaluation.

## Logistic-strength sweep

Changing coverage strength inside inference also traces a sharp Pareto front:

| Logistic weight | Figure 2D Pearson | D RMSE | Figure 2E Pearson |
|---:|---:|---:|---:|
| 0.10 | 0.8272 | 2.0967 | 0.6396 |
| 0.25 | 0.8271 | 2.0973 | 0.6396 |
| 0.35 | 0.8271 | 2.0977 | 0.6396 |
| 0.40 | 0.8076 | 2.2343 | 0.6573 |
| 0.45 | 0.8044 | 2.2588 | 0.6582 |
| 0.50 | 0.7763 | 2.4524 | 0.6919 |
| 0.75 | 0.7556 | 2.6126 | 0.6997 |
| 1.00 | 0.7526 | 2.6412 | 0.7013 |

No scalar weight recovers both objectives. The discontinuity near 0.4--0.5
also argues against trying to learn one global weight from a noisy single
sample diagnostic.

## Optimizer results

For the standard models, optimizer choice has negligible accuracy impact:
pooled Figure 2 correlations differ by less than 0.001 between ordinary EM,
SQUAREM, and corrected DAAREM. It has a large runtime impact.

Mean final-EM performance across six samples:

| Model | Optimizer | Evaluations | EM seconds | Peak RSS MiB |
|---|---|---:|---:|---:|
| None | ordinary | 611.5 | 24.87 | 196.2 |
| None | SQUAREM | 89.0 | 3.28 | 199.4 |
| None | DAAREM | 310.0 | 12.62 | 188.4 |
| Logistic | ordinary | 499.7 | 22.66 | 226.7 |
| Logistic | SQUAREM | 79.5 | 3.19 | 227.4 |
| Logistic | DAAREM | 181.7 | 7.96 | 233.2 |
| Auto | ordinary | 499.3 | 23.43 | 241.4 |
| Auto | SQUAREM | 75.5 | 2.82 | 241.3 |
| Auto | DAAREM | 61.0 | 2.32 | 241.5 |

SQUAREM is the most consistent choice: it reproduces ordinary-EM accuracy and
is about 7--8 times faster here. DAAREM is fastest for `auto`, but less
consistent across the other likelihoods.

Abundance blending is different because it initializes from a coverage-free
solution in a highly non-identifiable transcript equivalence class. After
making the warm-up use SQUAREM and preserving tiny positive support for every
transcript, ordinary EM and SQUAREM agree:

| Final optimizer | D Pearson | E Pearson | Evaluations | EM seconds |
|---|---:|---:|---:|---:|
| Ordinary | 0.7714 | 0.6974 | 202.7 | 6.68 |
| SQUAREM | 0.7713 | 0.6974 | 67.8 | 2.09 |
| DAAREM | 0.8268 | 0.6395 | 103.7 | 3.43 |

DAAREM selects a different point on a non-identifiable solution face and
effectively returns toward the coverage-free result. This is not a useful
accuracy gain because it discards the absolute-abundance coverage signal.

## Implementation changes retained

- DAAREM convergence is now checked against the fixed-point residual rather
  than the size of a damped accepted step. A tiny damped step is not evidence
  of convergence.
- Warm-started EM replaces exact zeros with a negligible positive floor.
  Exact zero is an absorbing state in multiplicative EM and can prevent a new
  likelihood from reactivating a transcript.
- The abundance warm-up now uses SQUAREM rather than DAAREM because it is
  faster and less path-dependent on this non-identifiable problem.
- The Kinnex runner accepts `--em-accel` and `--logistic-weight` for factorial
  reproduction.

## Decision

Do not add the stability-displacement gate or change PacBio `auto` based on
this dataset alone. Retain the EM correctness fixes and use SQUAREM for the
internal abundance warm-up. For user-selected final optimization, SQUAREM is
the current recommendation. A genuinely better coverage refinement needs
additional information beyond a scalar strength or final-count displacement,
most likely read-level censoring/quality evidence or cross-equivalence-class
predictive validation.
