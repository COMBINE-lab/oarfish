# Coverage-signal review (2026-07-21)

## Question

Which signals observable from the alignments and an initial abundance fit tell
Oarfish whether a coverage correction is needed, in which direction, and with
what strength?  The immediate constraints are to improve both Kinnex Figure 2D
(fold changes) and Figure 2E (absolute concentration), and to avoid disturbing
the rank order of low-abundance transcripts.

## What the Kinnex data say

The reusable analysis is `scripts/analyze_coverage_signals.py`; its row-level
output is `oarfish-evaluation-data/kinnex-wtc11/coverage-signal-analysis-20260721.tsv`.
Predictors use only ordinary Oarfish outputs: transcript length, no-coverage
abundance, total alignment support, unique-support fraction, and the proposed
auto correction. Truth is used only to define the desired correction.

The most important interaction is transcript length:

| Length | Fig. 2E auto win fraction | Fig. 2E correction alignment | Fig. 2D effect |
|---|---:|---:|---|
| <750 nt | 0.354 | -0.195 | essentially no change |
| 750--1249 nt | 1.000 | 0.400 | essentially no change |
| 1250--1999 nt | 0.641 | 0.469 | slight improvement |
| >=2000 nt | 0.500 | 0.510 | large regression |

“Correction alignment” is Spearman correlation between the correction proposed
by auto and the correction required by truth.  Thus, for absolute abundance,
the current model has real directional signal above 750 nt and its strongest
signal above 2 kb.  But the same correction is not stable between conditions:
for >=2 kb transcripts, mean absolute Figure 2D error rises from 1.682 (none) to
2.709 (auto).  A single per-sample correction can therefore improve absolute
calibration while introducing a condition-specific error that destroys fold
changes.

Support and ambiguity summaries do not solve this alone. Unique-support
fraction has near-zero association with correction benefit in both panels.
Total support and inferred abundance weakly predict Figure 2E benefit, largely
because both correlate with length and estimability. The magnitude of the
proposed correction is not a reliable confidence score.

These observations explain the earlier unique-read ablation: unique reads
carry useful ordering information, but unique support is not evidence that a
transcript-specific positional shape is estimated without bias.

## Signals supported by prior work

1. **Length-conditioned positional shape.** Existing positional-bias models
   pool or tie distributions by transcript length because the same normalized
   endpoint has different censoring implications for short and long molecules.
   Mix2 showed that multiple positional modes exist and that mixtures can be
   learned jointly with abundance. It also found that simple length-only tying
   can fail on real data, motivating partial rather than complete pooling.
2. **Two-sided endpoint completeness, not one coverage curve.** Long-read
   terminal errors are asymmetric. Direct RNA is sequenced 3'-to-5' and is
   particularly prone to 5' truncation; cDNA protocols add reverse-transcription
   and template-switching failure modes; PacBio and ONT also favor shorter
   molecules during sequencing. The 5' and 3' endpoint gaps should therefore be
   modeled separately and conditional on technology and molecule length.
3. **A survival/censoring process.** A read ending internally is not merely a
   low-density observation. It may represent RNA degradation, incomplete RT,
   pore termination, or an alignment error. A discrete competing-risk label
   (intact, 5'-censored, 3'-censored, internally broken) is more faithful than
   treating every end under one smooth positional density.
4. **Shared shape plus transcript deviation.** Per-transcript profiles are too
   noisy at low abundance; a single global profile erases real heterogeneous
   biases. Hierarchical shrinkage should learn technology/length-class shapes
   from high-information reads and allow transcript deviations only when the
   data support them.
5. **Equivalence-class contrast.** Coverage is useful only where candidates in
   the same read's equivalence class predict different observations. Absolute
   profile likelihood can penalize every candidate together and alter abundance
   without resolving ambiguity. Coverage evidence should be centered within
   each equivalence class and gated by its information gain (for example,
   entropy reduction or variance of candidate log likelihoods).
6. **Sequence and alignment quality.** Sequence context can explain substantial
   positional sampling variation in short-read models, and long-read terminal
   failures can be caused by structure, modifications, motor stalls, and
   basecalling/alignment errors. Available proxies include alignment-score
   margin, edit rate, soft clipping, homopolymer/GC context near observed ends,
   and read quality. These should primarily control confidence, not directly
   dictate abundance.
7. **Replicate/condition invariance.** A biochemical sampling bias should be
   considerably more stable across technical replicates than true abundance.
   The Figure 2D failure shows that independently learned sample profiles are
   not stable enough. When multiple samples are supplied, shared nuisance
   parameters should be learned jointly. In a single sample, an empirical prior
   trained across reference datasets can provide the same regularization.

## Recommended model

The next candidate should be a hierarchical competing-risks model, used only
for *relative evidence inside an equivalence class*:

1. Fit coverage-free abundances and retain them as the calibration anchor.
2. For each alignment, compute separate normalized 5' gap, 3' gap, aligned
   fraction, internal clipping/gaps, score margin, and read quality.
3. Estimate technology-by-length shared distributions for intact and censored
   molecules using high-information reads, with empirical-Bayes shrinkage.
4. Compute candidate log Bayes factors, center them within each equivalence
   class, and cap them according to posterior uncertainty.
5. Apply almost no coverage evidence when the class has low contrast, the
   transcript is low-support, or the proposed correction is not reproducible
   under read splitting.
6. Blend toward coverage-free abundance transcript-by-transcript. The blend
   depends on cross-fit stability and information gain, not abundance alone.

This construction protects low-abundance ranks in two ways: low-support shapes
are strongly pooled, and coverage cannot move a read when it supplies no
candidate-discriminating information. It also separates the useful Figure 2E
signal for long molecules from unstable per-condition profile estimates.

## Evaluation sequence

Before implementing the full likelihood, export alignment-level diagnostics
and run three falsifiable analyses:

1. Does endpoint/coverage evidence distinguish the true SIRV from alternatives
   within ambiguous equivalence classes more often than alignment score alone?
2. Are length-conditioned endpoint shapes reproducible across the three
   technical replicates and between day 0/day 5 after controlling for mixture?
3. Does split-half stability predict transcript-level oracle correction benefit
   on Kinnex, synthetic dRNA, and whole-sample-held-out LongBench data?

Promote a model only if it improves Figure 2D and 2E over no coverage, preserves
low-abundance Spearman, and improves held-out CCC/RMSE without a technology-wide
regression.

## Alignment-level diagnostic results

Oarfish now has an opt-in `--write-coverage-signals` diagnostic export. It
writes one deterministically sampled read in every
`--coverage-signal-sample-rate` reads, retaining the complete equivalence class:
candidate names, lengths, alignment intervals, strands, score probabilities,
and final coverage probabilities. The export does not alter inference and is
disabled by default. `scripts/analyze_alignment_coverage_signals.py` analyzes
the tables.

Six Kinnex samples contributed approximately 100,000 sampled reads each,
including 73,462--78,804 ambiguous reads. The central findings are:

| Signal | Six-sample mean |
|---|---:|
| Coverage/score winner agreement | 0.527 |
| Coverage agreement with high-confidence score winner | 0.579 |
| Combined likelihood agreement with high-confidence score winner | 1.000 |
| Mean coverage log Bayes-factor range | 0.092 |
| Coverage probability vs transcript length, Spearman | -0.615 |
| Coverage probability vs aligned fraction, Spearman | 0.313 |

PacBio auto coverage is consequently weak within individual equivalence
classes: it never changes the winner selected by a high-confidence alignment
score in this sample. Its dominant systematic signal is instead transcript
length. Repeating that small tilt over millions of reads can change abundance
without resolving read origin. This is a direct mechanism for improving an
absolute length bias while damaging fold changes.

Endpoint profiles are highly reproducible between random read halves of one
sample (length-class Jensen--Shannon divergence 0.0004--0.0035), so sampling
noise is not the main issue. They are less reproducible between samples,
especially for transcripts >=2 kb:

| Length | Mean between-sample JS | Maximum between-sample JS |
|---|---:|---:|
| <750 nt | 0.0366 | 0.1077 |
| 750--1249 nt | 0.0120 | 0.0227 |
| 1250--1999 nt | 0.0160 | 0.0281 |
| >=2000 nt | 0.0925 | 0.1964 |

Thus each sample has enough observations to estimate a profile precisely, but
the long-transcript nuisance profile varies across samples. Independently
applying those precise but different profiles explains the large Figure 2D
regression for >=2 kb SIRVs. The next candidate should share the nuisance shape
across related samples when possible and otherwise shrink toward a
technology-by-length empirical prior. It should also gate coverage by
candidate-relative information gain; evidence too weak to distinguish
candidates should be neutral.

At a 1-in-10 sampling rate, each million-read Kinnex table occupies about 20 MB.
The diagnostic name buffer raised peak RSS by roughly 60--80 MB and table/name
handling added approximately 2--6 seconds in these short runs. These costs are
strictly opt-in and are not paid during normal quantification.

## Candidate-relative contrast gate

The most direct implication was implemented and tested as an experimental
candidate: multiply coverage reliability by `r / (r + s)`, where `r` is the
within-equivalence-class coverage log Bayes-factor range and `s` is a tunable
half-strength scale. This makes zero-contrast evidence exactly neutral.

The Kinnex scale sweep demonstrated a smooth trade-off rather than a Pareto
improvement:

| Scale | Figure 2D Pearson | Figure 2E Pearson |
|---:|---:|---:|
| current auto | 0.7530 | 0.7015 |
| 0.10 | 0.7555 | 0.7000 |
| 0.25 | 0.7634 | 0.6971 |
| 0.35 | 0.7888 | 0.6756 |
| 0.40 | 0.8042 | 0.6584 |
| 0.50 | 0.8072 | 0.6575 |
| 1.00 | 0.8273 | 0.6396 |
| no coverage | 0.8273 | 0.6396 |

Increasing the gate scale simply interpolates from auto toward no coverage.
Scale 0.35 was carried to the outer panel as a compromise. Relative to current
auto, its mean changes were:

| Technology | CCC | Spearman | RMSE |
|---|---:|---:|---:|
| ONT cDNA | -0.00033 | -0.00044 | +4.96 |
| ONT dRNA | +0.00115 | -0.00154 | -12.02 |
| PacBio 50k | +0.00993 | -0.00205 | -75.10 |

Synthetic dRNA CCC fell from 0.99704 to 0.99232 and RMSE rose from 40.83 to
65.38, although Spearman rose from 0.82380 to 0.82667. Most decisively, on all
three PacBio 250k confirmations the gate was slightly worse than the current
PacBio no-endpoint auto kernel for CCC, Spearman, and RMSE. For example, H526
CCC changed from 0.33446 to 0.33195 and RMSE from 4510.27 to 4517.12.

The candidate was therefore removed. Contrast is a real measure of evidence
strength, but a per-read scalar cannot distinguish a reproducible biochemical
length correction from a sample-specific length tilt. The next experiment must
regularize the nuisance profile across samples or against a frozen empirical
technology prior.

## Literature

- Oarfish coverage motivation and observed dependence on fidelity/ambiguity:
  https://doi.org/10.1093/bioinformatics/btaf240
- Mix2 hierarchical positional mixtures and joint abundance/bias estimation:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5448817/
- Long-read terminal truncation mechanisms and technology asymmetry:
  https://genome.cshlp.org/content/34/11/1719
- Kinnex length bias and transcript-quantification benchmark:
  https://www.biorxiv.org/content/10.1101/2025.05.30.656561v3
- LRGASP protocol/depth/full-length benchmark:
  https://www.nature.com/articles/s41592-024-02298-3
- CapTrap-seq evidence that full-length enrichment is protocol dependent:
  https://www.nature.com/articles/s41467-024-49523-3
- Fragment GC/sequence bias and sample-specific correction:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5143225/
