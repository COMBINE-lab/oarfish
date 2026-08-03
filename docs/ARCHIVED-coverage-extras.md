# Archived: the four `auto` coverage extras (2026-08-03)

This branch preserves four experimental coverage features removed from mainline.
It is the last state in which they exist and in which the original selection
panel is reproducible. **Do not merge without new evidence.**

| feature | flag | status when archived | evidence |
|---|---|---|---|
| rank blending | `--rank-blend`, `--rank-blend-floor` | opt-in, **harmful** | worst family **-0.0878** Spearman |
| dominance pruning | `--candidate-pruning`, `--dominance-bayes-factor` | opt-in, **harmful** | worst family **-0.0437** |
| alignment calibration | `--alignment-calibration` | **default on**, unsupported | Panel B **-0.000143**, 2/6 samples |
| censoring | `--censoring-model` | **default on**, unsupported | Panel B **-0.000233**, 3/6 samples |

Full evidence: `docs/coverage-auto-rebenchmark-results-2026-08-03.md` and the
per-sample table `benchmarks/coverage-rebenchmark-2026-08-03.tsv` (448 runs).

## Why all four went at once

Neither of the two that were *enabled by default* can be supported. On Panel B —
the only panel with exact read-level truth at realistic depth, and the one
decision rule 5 prefers — both are net negative against the bare adaptive kernel.
Their apparent benefit on the union of panels (+0.0037 to +0.0043) comes almost
entirely from SIRV E2: four samples of 69 transcripts with four concentration
levels. Excluding SIRV they are +0.00017 and +0.00048.

For calibration: a separately evaluated feature, the splice-junction endpoint
term (`archive/junction-endpoint-2026-08-03`), was archived for being too weak at
**+0.000257** Spearman on Panel B with 6/6 samples improving. Both features that
shipped enabled perform *worse than that* on exact truth.

The two harmful features were already demoted to opt-in in `af8f3eb`. They are
removed here rather than left as unused surface area.

`censoring_probability::estimate_censor_scale` also fed rank blending's
activation gate, so the four were coupled: removing the two unsupported defaults
cleanly required removing the two harmful opt-ins as well.

## What mainline keeps

`--model-coverage` (logistic), the adaptive coverage kernel and its
technology dispatch, and the rest of `--coverage-model`. Note that the
re-benchmark also concluded (decision rule 4) that `auto` should not be the
default coverage model; that is a separate change from this removal.
