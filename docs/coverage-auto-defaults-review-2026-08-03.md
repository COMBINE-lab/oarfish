# Review of `--coverage-model auto` defaults (2026-08-03)

## Summary

An independent evaluation on NanoSim (ONT) and TKSM (PacBio) simulations found
that `--coverage-model auto` regressed accuracy relative to the published
`--model-coverage` (logistic) model. Ablation attributes essentially all of the
regression to two features that `auto` enabled implicitly — abundance rank
blending and dominance pruning — and not to the adaptive coverage kernel, which
matches or beats the logistic baseline once those two are disabled.

Both features are demoted to opt-in. They remain reachable
(`--rank-blend auto|fixed`, `--candidate-pruning auto|dominance`) so the panel
that originally selected them can be reproduced.

These defaults are **provisional**. The decisive experiment — re-running the
original selection panel with `logistic` reinstated as a comparator — has not
been run. See "Open question" below.

## Evaluation setup

Simulated data with exact read-level ground truth, generated independently of
oarfish development:

| Family | Datasets | Truth |
|---|---|---|
| NanoSim (ONT) | H9 and NA12878, each 1D-cDNA and direct-RNA | simulated read counts per transcript |
| TKSM (PacBio) | SQ2-HiFi, RSII | simulated read counts per transcript |

All comparisons use one fixed transcriptome BAM per dataset, so mapping
variation cannot contribute. `--filter-group no-filters` throughout. Metrics
computed against truth after accession-version stripping, outer-joined over
177,816 transcript identifiers.

## Accuracy

Spearman / MARD, alignment mode:

| Dataset | main `--model-coverage` | dev `auto` (previous defaults) | dev `auto` (kernel only) |
|---|---|---|---|
| NanoSim NA12878 1D-cDNA | 0.8850 / 0.0590 | 0.7677 / 0.1576 | 0.8846 / 0.0593 |
| NanoSim H9 1D-cDNA | 0.9061 / 0.0589 | 0.8099 / 0.1403 | **0.9077 / 0.0575** |
| NanoSim NA12878 dRNA | 0.9139 / 0.0396 | 0.8724 / 0.0538 | **0.9168 / 0.0382** |

"kernel only" disables all four auto-enabled extras. The adaptive kernel is not
the problem: it wins on two of three datasets and ties on the third, and
improves RMSE on NA12878 1D-cDNA (68.54 vs 72.09).

Two controls confirm the comparison is sound:

- `--coverage-model none` produces **identical** metrics on both branches
  (Spearman 0.8391, RMSE 227.4457, MARD 0.0873), in both alignment and read
  mode. The base pipeline and the read-mode mapper are unchanged.
- dev `--model-coverage` reproduces main `--model-coverage` (0.8867 vs 0.8850).
  The logistic path did not regress.

## Attribution

Add-one ablation on NanoSim NA12878 1D-cDNA, starting from the bare adaptive
kernel (Spearman 0.8846, MARD 0.0593), enabling exactly one feature at a time:

| Feature enabled | Spearman | Δ | MARD | wall |
|---|---:|---:|---:|---:|
| `--alignment-calibration agreement` | 0.8842 | -0.0004 | 0.0594 | 1:56.8 |
| `--censoring-model adaptive` | 0.8827 | -0.0019 | 0.0599 | 1:52.4 |
| `--candidate-pruning auto` | 0.8564 | **-0.0282** | 0.0713 | 2:20.5 |
| `--rank-blend auto` | 0.7879 | **-0.0967** | 0.1426 | 1:52.4 |

Leave-one-out ablation agrees: removing rank blending alone recovers the most
(0.7677 to 0.8382), and removing all four recovers fully (0.8846).

Alignment calibration and censoring are within noise on this data and are left
enabled; they are in scope for the re-benchmark below.

## Defects identified

**Rank blending.** `bulk.rs` blends each estimate toward a coverage-free warm-up
as `count = baseline + gate * (count - baseline)`, with `gate` floored at 0.8.
On NanoSim NA12878 1D-cDNA the midpoint lands at 2,870.87 counts, so `mean_gate`
is 0.8006 — pinned at the floor for essentially every transcript, i.e. a flat
20% pull toward the baseline rather than the intended abundance-dependent
protection. The warm-up itself is capped at `--coverage-warmup-iterations`
(default 100) and reports `"converged": false`, so 20% of every published
estimate came from a deliberately truncated EM.

Its activation gate requires a learned censoring scale above 37 nt.
`docs/rank-preserving-inference-2026-07-22.md` justifies that threshold on the
grounds that independent ONT/PacBio simulations fit the 25 nt lower clamp and
the rule would therefore abstain — the same document records that a global blend
"erased much of the coverage-model benefit" on exactly those samples. On NanoSim
direct-RNA the gate does abstain (25.0 nt). On NanoSim 1D-cDNA it learns
**87.85 nt** and fires. The premise behind the threshold does not hold on this
data.

**Dominance pruning.** `dominance_pruning.rs` sets `as_probabilities` to 0.0 for
any candidate Pareto-dominated on (score, coverage, alignment span) with a joint
margin of `--dominance-bayes-factor` (default 2.0). On NanoSim NA12878 1D-cDNA
this removes 7,592,011 of 34,899,437 candidates (21.8%) affecting 2,502,120
reads; on direct-RNA, 8,114,953 of 30,402,807 (26.7%). Requiring the winner's
span to be no shorter than the loser's systematically favors longer isoforms,
and the decision is made before the EM, removing ambiguity that abundance
evidence should resolve. It also costs ~28 s of wall time.

## Effect of the change

`--coverage-model auto` with the new defaults, NanoSim NA12878 1D-cDNA:

| | Spearman | MARD | RMSE | wall | peak RSS |
|---|---:|---:|---:|---:|---:|
| main `--model-coverage` | 0.8850 | 0.0590 | 72.09 | 3:23.9 | 2.51 GB |
| dev `auto`, previous defaults | 0.7677 | 0.1576 | 97.28 | 2:06.2 | 3.01 GB |
| dev `auto`, new defaults | 0.8824 | 0.0600 | **70.72** | **1:54.1** | 3.01 GB |

Parity on rank metrics, better RMSE, 1.8x faster.

## Performance note (no regression)

A reported runtime/memory regression did not reproduce. Controlled back-to-back
runs on an idle 32-core node, same inputs and flags:

| Mode | Config | main | dev | Δ time | Δ mem |
|---|---|---:|---:|---:|---:|
| alignment | `--model-coverage` | 3:23.9 / 2.51 GB | 1:52.4 / 3.01 GB | **-45%** | +20% |
| read | no coverage | 14:24.1 / 5.54 GB | 12:13.6 / 5.95 GB | **-15.1%** | +7.4% |
| read | `--model-coverage` | 14:52.3 / 5.65 GB | 13:10.9 / 6.11 GB | **-11.4%** | +8.3% |
| read | dev `auto` vs main logistic | 14:52.3 / 5.65 GB | 13:47.7 / 6.04 GB | **-7.2%** | +7.0% |

The memory increase is accounted for by two new `u32` clip fields per alignment
record (8 B x 34.9 M = 0.28 GB predicted, 0.23 GB observed), which feed the
packed parallel EM. Read mode gains less because it is dominated by mapping,
which is byte-identical between the branches.

Two cost items in the PacBio `auto` path remain unaddressed and are the largest
single outliers in the evaluation (`bulk.rs`, physical-endpoint branch): it
clones `store.coverage_probabilities` (`Vec<f64>`; 0.59 GB at TKSM SQ2's
73.5 M alignments) and holds it across a **full, unaccelerated** extra EM
(`max_iter: args.max_em_iter`, `EmAccel::None`). This is the only configuration
in the whole evaluation where dev is slower than main (+4%) and it carries the
largest memory delta (+15%).

## Open question

`logistic` was benchmarked once, in `docs/coverage-evaluation-2026-07-19.md`,
where it won decisively on the large synthetic benchmark (CCC 0.9950 vs endpoint
0.9552; RMSE 30.7 vs 91.5). Every subsequent comparison in this series —
hybrid, adaptive, degradation, censoring, rank blending, auto — reports only
`none`, `endpoint`, `hybrid`, and `adaptive`. The strongest baseline, and the
shipping default, was dropped from the comparison set and never reinstated.
Reported wins such as "27/28 versus no coverage" are therefore measured against
an opponent the features were never at risk of losing to.

Until the original panel is re-run with `logistic` as a comparator column, it is
not established that any of the `auto` extras — including the two left enabled
here — improve on the shipping default. The plan for that re-benchmark is in
`docs/coverage-auto-rebenchmark-plan-2026-08-03.md`.
