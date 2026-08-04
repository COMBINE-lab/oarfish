# Recalibrating `--score-prob-denom` (2026-08-03)

## Summary

The default denominator in the transcriptome score→probability conversion
`exp((score - best)/D)` changes from **5 to 3**. Measured over 34 samples across
both evaluation panels: **+0.0015 Spearman** on exact simulated truth (5/6
samples improving, MARD better on 6/6) and **+0.0007** on 28 real LongBench
samples (**27/28** improving), uniformly across ONT-cDNA, ONT-dRNA and PacBio.

This is a one-value change. No new code, flag, input or runtime cost.

## Why this parameter

The error decomposition in
`docs/coverage-auto-rebenchmark-results-2026-08-03.md` found that roughly a
quarter of misassigned reads have the **true source already scoring higher**
than the transcript the EM picks, and lose to the abundance prior anyway. That
is a calibration problem between alignment likelihood and abundance, not a
missing-feature problem, and `D` is exactly the knob controlling it: smaller `D`
sharpens the weighting toward the best-scoring alignment.

The remaining ~73% of misassignments are between transcripts whose alignments
score *identically*, and are unreachable by any per-alignment signal.

## Sweep

NanoSim NA12878 1D-cDNA, `--model-coverage`:

| `D` | Spearman | Δ vs 5 | MARD | RMSE |
|---|---|---|---|---|
| 1 | 0.879447 | -0.007271 | 0.060715 | 68.61 |
| 2 | 0.884568 | -0.002149 | 0.058781 | 68.58 |
| **3** | **0.887130** | **+0.000413** | **0.057870** | 68.59 |
| 5 (previous default) | 0.886717 | - | 0.058173 | 69.38 |
| 8 | 0.882268 | -0.004450 | 0.059805 | 71.47 |
| 12 | 0.875640 | -0.011078 | 0.062116 | 75.35 |
| 20 | 0.863546 | -0.023172 | 0.066372 | 83.71 |

The optimum is broad: 3 and 5 are both far better than 1 or 8. This moves toward
the centre of a flat region rather than to a sharp peak, so the change is not
fragile. Sharpening much further (`D=1`) is clearly harmful, which is why the
~27% bucket cannot simply be handed to the alignment score outright.

## Panel B - exact read-level truth

| sample | D=5 | D=3 | Δ Spearman | Δ MARD |
|---|---|---|---|---|
| nanosim-NA12878-cdna | 0.886717 | 0.887130 | +0.000413 | -0.000303 |
| nanosim-NA12878-drna | 0.919238 | 0.920082 | +0.000844 | -0.000295 |
| nanosim-H9-cdna | 0.908404 | 0.910829 | +0.002425 | -0.001085 |
| nanosim-H9-drna | 0.936348 | 0.936289 | -0.000058 | -0.000256 |
| tksm-RSII | 0.935371 | 0.938823 | +0.003452 | -0.002224 |
| tksm-SQ2 | 0.945946 | 0.947869 | +0.001924 | -0.002390 |
| **mean** | | | **+0.001500** (5/6) | **-0.001092** (6/6) |

## Panel A - 28 real LongBench samples (held out)

`D=3` was chosen by sweeping a Panel B sample, so all 28 experimental samples are
held out.

| group | n | Δ Spearman | improving | Δ MARD | improving |
|---|---|---|---|---|---|
| all | 28 | **+0.000721** | **27/28** | -0.000065 | 25/28 |
| ont-cdna | 8 | +0.000875 | 8/8 | -0.000055 | 8/8 |
| ont-drna | 9 | +0.000623 | 8/9 | -0.000047 | 7/9 |
| pac-bio | 11 | +0.000689 | 11/11 | -0.000086 | 10/11 |

Combined across both panels: **32 of 34 samples improve**.

## Why this cleared the bar when other candidates did not

Three features evaluated in the same series were rejected or removed:
the splice-junction endpoint term (+0.000257 on Panel B, +0.000067 experimental),
alignment calibration (-0.000143 on Panel B) and censoring (-0.000233). Two
differences matter here:

1. **It generalises.** The junction term's experimental effect was 4x smaller
   than its simulated one and its consistency dropped to 20/28. This change
   halves in magnitude but *gains* consistency, to 27/28 on fully held-out real
   data.
2. **It costs nothing.** Every other rejection turned on whether a feature was
   worth carrying. Nothing is carried here; a single default value changes.

## Caveats

- The effect is small in absolute terms (+0.0007 to +0.0015 Spearman). No
  biological conclusion changes.
- Changing a default shifts results slightly for every existing user, which
  matters for reproducing previously published numbers. `--score-prob-denom 5`
  restores the old behaviour exactly.
- Both panels use `--model-coverage`. The parameter is transcriptome-mode only;
  genome mode weights projected alignments by bramble similarity instead.
