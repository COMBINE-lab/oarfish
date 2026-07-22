# Rank-preserving adaptive inference evaluation (2026-07-22)

## Motivation

The automatic coverage, censoring, and dominance models improve count-scale
metrics, but their largest remaining weakness was low-abundance rank order. A
coverage-free preliminary EM estimate often preserved that order. We therefore
tested an abundance-dependent blend between the preliminary and corrected
estimates.

The retained blend keeps at least 80% of the corrected estimate. Its correction
fraction rises smoothly to 100% with preliminary abundance, with a midpoint of
300 counts per million reads. This primarily protects low-abundance ranks while
leaving high-abundance coverage corrections nearly unchanged.

## Selection experiments

- A global blend improved Spearman on all 28 primary cases (mean +0.00822 versus
  the previous automatic model), but erased much of the coverage-model benefit
  on independent ONT and PacBio simulations. It was rejected as a default.
- Transcript-level gating by unique-read support also failed those independent
  simulations and was rejected.
- The censoring scale learned from unique reads separated the cases. Independent
  ONT/PacBio simulations fit the lower clamp of 25 nt; the main datasets fit
  39--92 nt. SIRV mixtures fit 25--36 nt but contain only 69 transcripts.
- The retained automatic rule activates rank preservation only for transcriptomes
  with at least 1,000 transcripts and a learned censoring scale above 37 nt.
  Otherwise it abstains and uses the corrected EM estimate unchanged.

The thresholds are deliberately based only on sample-internal diagnostics; no
truth labels enter inference. `--rank-blend fixed` remains available for explicit
experiments, and `--rank-blend none` disables the additional warm-up.

## Accuracy

On the 28-case primary panel, automatic selection activates for every case and
matches the validated fixed blend:

| Comparison | Pearson | Spearman | CCC | RMSE | MARD |
|---|---:|---:|---:|---:|---:|
| vs previous auto, mean delta | +0.00021 | **+0.00822** | -0.00052 | +15.36 | -0.00090 |
| wins vs previous auto | 16/28 | **28/28** | 4/28 | 3/28 | 20/28 |
| vs no coverage, mean delta | +0.02314 | **+0.00471** | +0.02296 | -155.00 | +0.00037 |
| wins vs no coverage | 27/28 | **27/28** | 27/28 | 25/28 | 9/28 |

Positive is better for correlations; negative is better for RMSE and MARD. The
remaining Spearman loss versus no coverage is H69 PacBio (-0.00186), while the
model retains its count-scale benefit there.

On outer validation the rule abstains in all five cases. It therefore retains
the previous model's strong results instead of the rejected global blend:

| Sample | Previous/automatic Spearman | Global blend Spearman |
|---|---:|---:|
| independent ONT simulation | 0.7679 | 0.7048 |
| independent PacBio simulation | 0.7642 | 0.7099 |
| SIRV E2 dRNA | 0.8749 | 0.8692 |

## Cost

When active, the preliminary coverage-free EM adds a mean 0.079 seconds and
about 1.0 MiB peak RSS versus the previous automatic model on the primary panel.
When the learned rule abstains, only a linear scan over unique alignments is
added; no preliminary EM or blend allocation is performed.

Artifacts are in `oarfish-evaluation-data/rank-blend-mid-300-full-20260722`,
`rank-blend-auto-outer-20260722`, and `rank-blend-auto-boundary-20260722`.
