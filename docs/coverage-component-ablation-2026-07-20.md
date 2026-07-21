# Coverage-component ablation: 2026-07-20

## Question and frozen design

This study asks which parts of the adaptive coverage likelihood are useful,
whether three annotation-free extensions improve it, and what each part costs.
No transcript-to-gene annotation was used. The control was the `auto` model at
commit `312da81`: adaptive logistic/endpoint evidence for ONT cDNA and PacBio,
and the same evidence after competing-risk correction for ONT dRNA.

The screen used one fixed synthetic dRNA BAM with transcript truth and 24
LongBench libraries (eight cell lines each for ONT cDNA, ONT dRNA, and PacBio;
50,000-read prefixes). LongBench matched Illumina estimates are an external
biological comparator, not molecule-level truth. Three available PacBio
250,000-read prefixes were a confirmation tier. Every arm reused the identical
BAM, filters, EM settings, and four threads.

Accuracy includes CCC, Pearson, Spearman, RMSE, and MARD. `/usr/bin/time -v`
recorded whole-process wall time and peak RSS; oarfish metadata recorded
alignment, coverage-model, and EM time separately. The reproducible driver is
`scripts/run_coverage_ablation.py`, with inputs in
`benchmarks/coverage_ablation_manifest.tsv`.

## Arms

The leave-one-component-out arms removed logistic evidence, endpoint evidence,
the empirical-support gate, logistic/endpoint agreement gate, or per-read
Bayes-factor cap. `none` supplied the no-coverage reference.

Three annotation-free additions were tested independently and together:

- **quality gate:** shrink coverage evidence according to the normalized
  separation of the best alignment-score probability from the alternatives;
- **uncertainty gate:** shrink according to endpoint-probability variation
  among the five cross-fit models; and
- **equivalence-class training:** add ambiguous reads to endpoint training,
  fractionally weighted by their existing alignment-score probabilities.

These experimental arms remain opt-in and never affect `auto` unless explicitly
selected with `--coverage-ablation`.

## Large synthetic truth result

| Arm | CCC | RMSE | MARD | Coverage time | Peak RSS |
|---|---:|---:|---:|---:|---:|
| full | **0.997038** | **40.83** | 0.33991 | 2.709 s | 359.1 MB |
| no logistic | 0.972790 | 122.78 | 0.32889 | 2.734 s | 359.3 MB |
| no endpoint | 0.994975 | 53.12 | 0.34125 | 2.718 s | 359.2 MB |
| no support gate | 0.997076 | 40.57 | 0.33995 | 2.718 s | 359.4 MB |
| no agreement gate | 0.997054 | 40.72 | **0.33983** | 2.745 s | 359.4 MB |
| no Bayes cap | 0.997038 | 40.83 | 0.33984 | 2.728 s | 359.1 MB |
| quality gate | 0.957863 | 152.65 | 0.30243 | 2.720 s | 359.1 MB |
| uncertainty gate | 0.997038 | 40.83 | 0.33991 | 2.759 s | 385.9 MB |
| equivalence-class training | 0.996458 | 44.63 | 0.33983 | 2.737 s | 363.0 MB |
| all additions | 0.957361 | 153.56 | **0.30234** | 2.772 s | 389.9 MB |

The quality gate is a decisive failure for abundance calibration even though it
improves MARD. Alignment ambiguity is precisely where coverage evidence is
needed; treating ambiguity itself as evidence that coverage is unreliable
removes useful correction. This is also a warning against selecting a model on
one error statistic.

Fold uncertainty changes no displayed truth metric and costs 26.7 MB, exactly
the expected extra eight-byte value per 3.21 million retained alignments, plus
about 50 ms. Equivalence-class training is inexpensive but slightly worsens
synthetic CCC and RMSE.

## LongBench screen

The table reports mean CCC change from the full control within each technology.
Ranges are across eight cell lines.

| Arm | ONT cDNA | ONT dRNA | PacBio |
|---|---:|---:|---:|
| no coverage | -0.006702 | -0.006968 | -0.008326 |
| no logistic | -0.000590 | -0.001069 | +0.006139 |
| no endpoint | -0.002366 | -0.001272 | **+0.012212** |
| no support gate | -0.000255 | +0.000011 | +0.001359 |
| no agreement gate | -0.000123 | -0.000156 | +0.000362 |
| no Bayes cap | +0.000092 | -0.000043 | -0.001934 |
| quality gate | -0.000049 | -0.002218 | -0.002305 |
| uncertainty gate | +0.000003 | -0.000004 | -0.000003 |
| equivalence-class training | -0.000216 | +0.000066 | +0.001850 |
| all additions | -0.000096 | -0.002137 | -0.001228 |

The full model beats no coverage in every ONT cDNA and dRNA library. Logistic
and endpoint evidence are both useful for ONT; removing either lowers mean CCC.
The support, agreement, and Bayes-factor safeguards have small mean effects but
prevent sample-specific losses, so their negligible cost argues for retaining
them.

PacBio behaves differently. Equivalence-class training improves all eight
50,000-read libraries, but removing endpoint improves seven of eight and has a
much larger mean gain. This suggests a technology mismatch, not insufficient
endpoint training.

## PacBio 250,000-read confirmation

| Library | none | full endpoint | equivalence-class | logistic without endpoint |
|---|---:|---:|---:|---:|
| H69 | 0.512762 | 0.564665 | 0.564990 | **0.589353** |
| H146 | 0.313482 | 0.322583 | 0.322731 | **0.323708** |
| H526 | 0.285937 | 0.319979 | 0.320096 | **0.334463** |

The simpler PacBio kernel wins all three confirmations. Its coverage-model time
is 0.136--0.162 s, indistinguishable from the old full model's 0.136--0.165 s;
peak RSS is likewise indistinguishable at about 494--496 MB. Whole-process wall
times vary by more than the coverage difference because all arms run the same
non-converged 1,000-iteration EM, so component time is the cleaner cost measure.

Across the 24-library screen, full coverage construction averaged 60.0 ms.
Quality gating and equivalence-class training added less than 1 ms on average;
uncertainty gating added 3.0 ms. Median RSS deltas among 50,000-read libraries
were below 0.3 MB and dominated by process noise. The large synthetic allocation
provides the reliable memory measurement for uncertainty gating.

## Decisions

1. Keep the current ONT cDNA adaptive endpoint/logistic model and ONT dRNA
   competing-risk model, including support, agreement, and Bayes-factor gates.
2. Change `auto` for PacBio and PacBio HiFi to the adaptively gated logistic
   kernel without endpoint evidence. The global CLI default remains `none`.
3. Do not enable the quality or uncertainty gates. Quality is harmful;
   uncertainty has measurable memory cost and no accuracy benefit.
4. Do not enable equivalence-class endpoint training. Its small PacBio gain is
   superseded by removing the mismatched endpoint term, and it does not
   generalize to synthetic dRNA or ONT cDNA.
5. Do not yet implement sequence-derived neighborhoods or a fully alternating
   abundance/coverage EM. The lightweight equivalence-class result sets a weak
   upper-priority signal; a future joint model should first have independent
   PacBio molecular truth and an explicit monotonic held-out likelihood gate.

The new PacBio `auto` output was checked bit-for-bit against the corresponding
`no-endpoint` ablation (maximum count difference `0.0` over 385,659 targets),
and its metadata records `pacbio_adaptive_logistic` and the effective ablation.

## Artifacts

Screen outputs are in
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/coverage-ablation-screen-20260720`;
the 250,000-read confirmation is in
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/coverage-ablation-confirm-250k-20260720`.

## Allocation optimization and rank follow-up

The adaptive implementation originally held a cloned full-length logistic
vector, full-length endpoint probabilities, per-alignment duplicated support
gates, and a second full-length combined output at the same time. The optimized
implementation now:

- combines into the alignment store's existing coverage vector in place;
- stores support and optional uncertainty gates once per read as `f32`, rather
  than once per alignment as `f64`; and
- releases cross-validation observations before allocating per-alignment
  endpoint output.

On the 3,214,282-alignment synthetic benchmark, peak RSS fell from 359,148 KB
to 303,932 and 308,132 KB in two optimized runs. The optimized mean is about
53 MB (14.8%) below the original peak, and reduces the incremental memory over
the 216,592 KB no-coverage run by approximately 37%. Coverage construction fell
from 2.709 seconds to 2.573 seconds. The largest target-count difference was
`9.34e-6`, total absolute difference was `3.73e-4`, and all accuracy metrics
were unchanged at reported precision.

The simulated Spearman decrease is not an unavoidable consequence of improved
calibration. As a diagnostic only, convexly mixing no-coverage and full-model
counts raised Spearman from 0.84319 to 0.84764 at 20% full correction, although
CCC was then only 0.97012. More promisingly, retaining full correction for
well-supported transcripts while applying 25% correction below a 50-read
preliminary-abundance threshold produced CCC 0.99706, Spearman 0.84640, RMSE
40.68, and MARD 0.29837. This is post-hoc feasibility evidence, not a deployable
estimator.

A cheap inference-level substitute based only on unambiguous endpoint-training
support was implemented and rejected: it improved Spearman only from 0.82380
to 0.82597 while lowering CCC to 0.99596 and increasing RMSE to 47.57. The code
was removed. Recovering the full Pareto improvement therefore appears to require
an actual preliminary abundance estimate. The next candidate should be a short,
accelerated coverage-free warm-up followed by abundance-gated coverage EM,
with the warm-up time included in the benchmark and sample-level validation
used to learn the shrinkage curve.

### Two-stage abundance blend

The proposed two-stage estimator was subsequently implemented. It runs 100
DAAREM-accelerated coverage-free evaluations, uses that estimate to initialize
the ordinary coverage EM, then shrinks each final target count toward its
warm-up count. Shrinkage has a 25% floor and a smooth fourth-power transition
centered at 37 warm-up counts per million aligned reads. Counts are renormalized
to the aligned library size. This is available only through the experimental
`--coverage-ablation abundance-blend`; it is not selected by `auto`.

On synthetic dRNA it recovers more than the lost rank performance while also
slightly improving calibration:

| Model | CCC | Spearman | RMSE | MARD |
|---|---:|---:|---:|---:|
| no coverage | 0.956779 | 0.843190 | 154.61 | 0.29892 |
| full coverage | 0.997038 | 0.823800 | 40.83 | 0.33991 |
| abundance blend | **0.997070** | **0.859898** | **40.59** | **0.29448** |

The warm-up adds 1.78 seconds and approximately 4 MB on the large simulation.
Warm-up initialization alone lowered Spearman to 0.82327. Applying preliminary
abundance inside the read likelihood raised it only to 0.82458 at the frozen
midpoint, and stronger likelihood shrinkage reached only 0.82534 while worsening
RMSE. Those two rejected paths were removed.

Across the 24 LongBench 50,000-read libraries, the final count blend improved
Spearman in every sample relative to the current technology-specific `auto`
model. Mean Spearman changes were +0.00922 for ONT cDNA, +0.00732 for ONT dRNA,
and +0.00680 for PacBio. It also slightly exceeded no coverage in all three
technologies (mean +0.00079, +0.00074, and +0.00075 respectively).

This rank recovery gives back some of the full model's calibration. Relative to
current `auto`, mean CCC changed by -0.00194 for cDNA, -0.00339 for dRNA, and
-0.00108 for PacBio; mean RMSE increased by 54.8, 42.5, and 13.2 in the
LongBench comparator scale. Nevertheless, compared with no coverage it retained
positive mean CCC in all technologies and improved RMSE in most libraries.
Whole-process wall time increased by about 0.52--0.56 seconds at 50,000 reads,
with median peak-RSS increases of roughly 1.5--3.3 MB.

Decision: retain the abundance blend as an experimental accuracy/rank tradeoff,
but do not make it the default. The synthetic result is unusually strong, while
the biological comparator shows that a single frozen shrinkage curve can
over-regularize calibration. A default would require sample-level learning of
blend strength against a truth-independent objective and validation on
molecular mixtures with transcript truth.

### Attempted learned and transcript-dynamic weights

Two empirical-Bayes versions were tested to replace the fixed blend curve.
Both used transcript-specific log corrections between warm-up and coverage
counts. A spike-and-slab prior learned the active-correction fraction and slab
scale within each sample; posterior weights therefore varied dynamically with
transcript abundance, correction magnitude, and estimated uncertainty.

The count-noise approximation recovered part of the simulated rank loss
(CCC 0.997056, Spearman 0.836214, RMSE 40.71), but depended on an assumed
paired-noise scale. Across LongBench it improved Spearman in all libraries yet
gave back nearly the same CCC as the fixed blend. Learning that noise scale by
in-sample marginal likelihood selected the smallest candidate and collapsed
simulated Spearman to 0.825172.

A second implementation estimated transcript-specific correction variance from
five deterministic read folds in a linear responsibility pass. It dynamically
incorporated abundance, ambiguity, coverage likelihood, and fold stability.
The fitted model assigned mean posterior weight 0.483 because the coverage
correction was highly stable, but simulated Spearman was only 0.823595. This
demonstrates the central identifiability problem: a reproducible correction can
still be biased relative to abundance truth. Neither held-out likelihood nor
within-sample stability can discover that bias by itself.

Both empirical-Bayes implementations were removed after evaluation. A more
adaptive blend remains promising, but its mapping from observable features to
weights must be calibrated on truth-bearing samples with entire-sample outer
holdouts. Candidate features should include preliminary abundance, ambiguous
and unique support, coverage reliability, warm-up/full disagreement,
degradation posterior, and fold variance. Until such a training panel is large
enough, the fixed abundance blend is the only retained experimental rank-aware
option and `auto` remains unchanged.
