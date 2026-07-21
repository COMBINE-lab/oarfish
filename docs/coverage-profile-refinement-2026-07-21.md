# Coverage-profile refinement study (2026-07-21)

## Goal and evaluation policy

This study followed the Kinnex Figure 2 result by testing richer, annotation-free
ways to prevent coverage overcorrection. Candidate parameters were learned only
from each sample. Kinnex Day 0/Day 5 truth was used for evaluation, followed by
an outer panel consisting of the synthetic dRNA truth set, 24 LongBench
50,000-read libraries, and three PacBio 250,000-read confirmations. Candidates
were removed unless the gain generalized across objectives and datasets.

All new runs used SQUAREM for the final EM. This optimizer reproduces ordinary
EM on the standard models while greatly reducing runtime.

## Explicit PacBio censoring

Three increasingly explicit interpretations were tested:

1. shrink all coverage evidence when no candidate alignment was terminally
   anchored;
2. neutralize coverage evidence candidate-by-candidate according to total 5'
   and 3' terminal gap; and
3. learn an intact-read fraction from the sample, then score candidates under
   an intact-versus-censored terminal-gap mixture.

The read-level gate was almost always one when defined by either end (mean
0.9975). Requiring both ends produced a mean gate near 0.896 but changed Kinnex
only slightly (D Pearson 0.7549, E Pearson 0.7002). Candidate neutralization
collapsed to the no-coverage result.

The generative mixture learned intact fractions of 0.872--0.899 and slightly
improved Kinnex Figure 2E to 0.7051, but Figure 2D remained 0.7516. On the eight
LongBench PacBio libraries it reduced CCC by 0.0171 on average and increased
RMSE by 113.2; all three 250k confirmations regressed. Endpoint completeness
therefore captures plausible full-length isoform fit but not unbiased abundance.
All censoring variants were removed.

## Equivalence-class-aware coverage profiles

The current logistic profile construction counts every retained alignment with
weight one, so an ambiguous read contributes to every candidate transcript.
Three alternatives were evaluated.

### Coverage-free abundance responsibilities

A 100-evaluation SQUAREM coverage-free fit supplied transcript abundances.
Coverage bins were rebuilt from alignment-score times abundance posterior
responsibilities before coverage-aware inference.

On Kinnex this recovered much of the fold-change loss (D Pearson 0.8025) while
retaining some absolute signal (E Pearson 0.6599). It failed every outer
technology group:

| Outer group | Mean CCC change | Mean Spearman change | Mean RMSE change |
|---|---:|---:|---:|
| ONT cDNA | -0.00332 | -0.00251 | +42.0 |
| ONT dRNA | -0.00267 | -0.00256 | +18.6 |
| PacBio | -0.02493 | -0.00190 | +158.3 |

All three 250k PacBio confirmations also regressed. Feeding a biased
coverage-free abundance estimate back into profile construction is circular
and amplifies its assignment errors. The candidate was removed.

### Fractional alignment-score profiles

Each read contributed total profile mass one across its candidates, divided by
alignment-score probabilities, without abundance feedback. It made a small
Kinnex Pareto improvement over current auto (D 0.7565 versus 0.7526; E 0.7020
versus 0.7013). This did not survive the outer panel: Spearman fell in 23/24
LongBench libraries, synthetic RMSE rose from 40.83 to 45.05, and two of three
250k confirmations regressed. The candidate was removed.

### Unique-read-only profiles

Logistic profiles were rebuilt only from reads with one retained transcript
alignment. This strongly improved Kinnex Figure 2D (Pearson 0.8439, RMSE 1.955)
and slightly exceeded no coverage on Figure 2E (0.6425 versus 0.6396).

It also improved Spearman and MARD in every LongBench library, demonstrating
that ambiguous-read duplication contributes to rank distortion. However,
calibration was unacceptable:

- synthetic CCC fell from 0.99704 to 0.94994;
- synthetic RMSE rose from 40.83 to 167.61;
- mean PacBio LongBench CCC fell by 0.02235;
- all three 250k PacBio confirmations lost CCC and gained RMSE.

This is a rank-oriented shrinkage model, and the existing abundance blend
achieves that trade-off more safely. The candidate was removed.

## Runtime and memory

The profile-only candidates added little direct coverage-construction time at
50,000 reads (roughly 0.06--0.14 seconds depending on technology) and no clear
RSS increase beyond process noise. The responsibility candidate additionally
cost approximately 3.15 seconds per one-million-read Kinnex sample for its
warm-up. Accuracy, not computation, rejected these models.

## Conclusions

No new coverage model from this round is promoted into `auto`. The negative
results establish three useful constraints:

1. Terminal completeness is not a sufficient proxy for trustworthy abundance
   correction.
2. Ambiguous reads should not train coverage profiles through in-sample
   abundance feedback.
3. Unique reads expose a real rank/calibration trade-off, but optimizing rank
   alone can severely damage abundance calibration.

The next model should keep cross-fitted unique-read training but learn a
calibrated correction mapping on a multi-sample truth panel with entire-sample
outer holdouts. Candidate features can include unique support, transcript
length, ambiguity degree, endpoint completeness, correction magnitude, and
technology. Without that external calibration, the current `auto` plus the
experimental abundance blend remain safer than an increasingly elaborate
single-sample heuristic.

Raw outputs are under `oarfish-evaluation-data` in directories prefixed
`censoring-*`, `responsibility-coverage-*`, `fractional-coverage-*`, and
`unique-coverage-*`.
