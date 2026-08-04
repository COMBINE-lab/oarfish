# Retiring the non-logistic coverage models (2026-08-03)

## Summary

Every coverage model oarfish offers was run head-to-head on Panel B — six
NanoSim/TKSM simulations with exact read-level truth — using the current binary,
serially on an idle node. **`--model-coverage` (logistic) beats all of them.**

`--coverage-model auto|adaptive|endpoint|hybrid|degradation` and their four
supporting modules are removed. `--coverage-model none|logistic` remain. The
retired stack is preserved on `archive/coverage-kernels-2026-08-03` together
with this evidence.

## Results

Spearman against exact truth, scored over an identical reference transcriptome:

| sample | none | logistic | auto | endpoint | hybrid | degradation |
|---|---|---|---|---|---|---|
| nanosim-NA12878-cdna | 0.840612 | **0.887130** | 0.885205 | 0.820957 | 0.877806 | n/a |
| nanosim-NA12878-drna | 0.873359 | **0.920082** | 0.917575 | 0.848241 | 0.911196 | 0.917575 |
| nanosim-H9-cdna | 0.876814 | **0.910829** | 0.910073 | 0.857114 | 0.905001 | n/a |
| nanosim-H9-drna | 0.911101 | **0.936289** | 0.934951 | 0.878742 | 0.931588 | 0.934951 |
| tksm-RSII | 0.932431 | 0.938823 | 0.940940 | 0.928092 | **0.941208** | n/a |
| tksm-SQ2 | 0.950101 | 0.947869 | 0.949717 | 0.951190 | **0.952083** | n/a |

Paired against `logistic`:

| arm | mean Spearman | Δ vs logistic | better on | mean wall |
|---|---|---|---|---|
| `none` | 0.897403 | -0.026101 | 1/6 | 111.3 s |
| **`logistic`** | **0.923504** | — | — | **115.0 s** |
| `auto` | 0.923077 | -0.000427 | 2/6 | 144.8 s (**+25.9%**) |
| `endpoint` | 0.880723 | **-0.042781** | 1/6 | 116.0 s |
| `hybrid` | 0.919814 | -0.003690 | 2/6 | 117.2 s |
| `degradation` | — | -0.001923 | 0/2 | 109.9 s |

Notes on reading this table:

- **`endpoint` is worse than no coverage model at all** (0.880723 vs `none`'s
  0.897403), losing to `none` on 5 of 6 samples.
- **`degradation` is byte-identical to `auto` on its two applicable samples.**
  It is what `auto` dispatches to for ONT direct-RNA. It requires
  `--seq-tech ont-drna` and errors on cDNA and PacBio, so it is only evaluable
  on 2 of 6 samples; its higher raw mean is an artifact of those two samples
  having the highest absolute Spearman in the panel.
- **Coverage modelling itself is clearly worthwhile**: `logistic` - `none` is
  **+0.026** on 5/6 samples. This retirement removes alternatives to logistic,
  not coverage modelling.

## The PacBio caveat, stated explicitly

On both TKSM PacBio samples, `hybrid` and `auto` beat `logistic`, and on
TKSM SQ2 even `none` beats it (0.950101 vs 0.947869). This is the same pattern
seen throughout this evaluation series: `logistic` underperforms on PacBio.

That is not sufficient to keep the stack — n=2, no kernel wins overall, and the
gains are 0.002-0.004 — but it is a real, repeatedly-observed signal and the
most likely place a future coverage model would find headroom. It is recorded
here rather than dropped.

## What was removed and what was kept

Dependencies were verified rather than assumed. `binomial_probability` initially
appeared orphaned to a module-path search, but is reached from `single_cell.rs`
via a crate-root re-export; deleting it would have broken single-cell
quantification. `single_cell.rs` itself does not touch the adaptive stack.

| | lines | |
|---|---|---|
| removed | 388 | `util/endpoint_probability.rs` |
| removed | 271 | `util/hybrid_probability.rs` |
| removed | 545 | `util/degradation_probability.rs` |
| removed | 274 | `util/pacbio_endpoint_probability.rs` |
| **removed total** | **1,478** | |
| kept | 81 | `util/logistic_probability.rs` (logistic path) |
| kept | 76 | `util/normalize_probability.rs` (logistic path) |
| kept | 224 | `util/binomial_probability.rs` (single-cell path) |

## Relationship to the rest of this series

This completes a sequence in which every addition to the coverage model made
after `logistic` was measured against it on exact truth and failed:

| change | Δ vs logistic on exact truth | outcome |
|---|---|---|
| rank blending | -0.0878 worst family | removed |
| dominance pruning | -0.0437 worst family | removed |
| alignment calibration | -0.000143 | removed |
| censoring | -0.000233 | removed |
| splice-junction endpoint | +0.000257 | archived, not adopted |
| `endpoint` kernel | -0.042781 | removed here |
| `hybrid` kernel | -0.003690 | removed here |
| `auto` / `adaptive` kernel | -0.000427, +25.9% cost | removed here |
| `degradation` kernel | -0.001923 | removed here |

The one change that *did* survive was not a feature at all: recalibrating
`--score-prob-denom` from 5 to 3, which improved 32 of 34 samples at zero
structural cost (`docs/score-prob-denom-recalibration-2026-08-03.md`).

The root cause is documented in
`docs/coverage-auto-rebenchmark-results-2026-08-03.md`: these features were
selected against `none`/`endpoint`/`hybrid`/`adaptive` on a panel whose truth is
matched-Illumina, and on which `none` outranks `logistic`. The strongest
baseline was never in the comparison set.
