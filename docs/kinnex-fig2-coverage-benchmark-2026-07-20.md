# Kinnex Figure 2 coverage-model benchmark (2026-07-20)

## Data and design

The author-provided `kinnex_wtc11_oarfish.tar.gz` archive contains 1-million-read
transcriptome-aligned BAMs for the three Day 0 (SIRV E1) and three Day 5 (SIRV
E2) replicates, the SIRV Set 4 reference, and archived Oarfish counts. We ran
the current `dev` release build with eight threads and no alignment filters.

The compared modes were:

- `none`: no positional coverage likelihood;
- `logistic`: the current repaired logistic implementation;
- `auto`: the current PacBio adaptive-logistic kernel;
- `abundance-blend`: `auto` followed by transcript-abundance-dependent blending
  with a coverage-free warm-up estimate.

Panel D is pooled across 69 SIRVs and three matched replicate pairs (207
observations): observed Day0-Day5 logCPM difference versus expected
`log2(E1/E2)`. Panel E pools the three Day 0 replicates: logCPM versus
`log(molecular_weight * E1_molarity)`. CPM values use edgeR's default prior
count convention. Pearson is the paper's displayed statistic; Spearman and
Panel-D RMSE are additional diagnostics.

## Results

| Model | D Pearson | D Spearman | D RMSE | E Pearson | E Spearman |
|---|---:|---:|---:|---:|---:|
| Archived paper Oarfish | 0.8055 | 0.8286 | 2.2472 | 0.7044 | 0.7354 |
| None | **0.8273** | **0.8578** | **2.0962** | 0.6396 | 0.6706 |
| Current logistic | 0.7515 | 0.8192 | 2.6568 | **0.7029** | **0.7320** |
| Current auto | 0.7530 | 0.8202 | 2.6383 | 0.7015 | 0.7303 |
| Abundance blend | 0.7886 | 0.8326 | 2.3797 | 0.6776 | 0.7027 |

The effects are concentrated in long SIRVs (at least 1250 nt). Short-SIRV
Panel-D Pearson is approximately 0.887 for every mode. For long SIRVs, Panel-D
Pearson is 0.778 with no model, 0.651 with auto, and 0.711 with abundance
blending. Conversely, long-SIRV Panel-E Pearson is 0.281 with no model, 0.669
with auto, and 0.449 with abundance blending.

Thus this benchmark exposes a real objective trade-off. Coverage correction
recovers the length/mass-associated absolute-abundance relationship in Panel E,
but overcorrects relative E1/E2 changes in Panel D. The adaptive reliability
gate in current PacBio `auto` barely changes the current logistic result. The
abundance blend moves toward no coverage and substantially recovers Panel D,
but gives up too much Panel-E accuracy to dominate either endpoint.

The archived counts reproduce the paper's Panel-D value (0.8055, displayed as
0.81) and approximately reproduce Panel E (0.7044 versus displayed 0.71). A
fresh run of current `logistic` does not exactly reproduce those archived
counts. This is expected because the current implementation is repaired and
other Oarfish internals have evolved; the archived arm is therefore retained
as a historical reference rather than treated as an executable mode.

## Runtime and memory

Means across six samples, using Oarfish's phase timers:

| Model | Alignment + coverage + EM (s) | Coverage phase (s) | Mean peak RSS (MiB) | Median EM evaluations |
|---|---:|---:|---:|---:|
| None | 27.26 | ~0 | 196.2 | 588.5 |
| Current logistic | 25.17 | 0.051 | 226.7 | 464.5 |
| Current auto | 25.96 | 0.237 | 241.4 | 464.5 |
| Abundance blend | 11.76 | 3.635 | 239.5 | 164.0 |

The blend is faster here despite its warm-up because that warm-up uses DAAREM
and supplies a strong initialization to the final EM. Three Day-5 runs reached
the 1000-evaluation limit in several non-blend modes, so convergence behavior
is a material part of the timing result.

## Reproduction

Quantifications, logs, metadata, resource measurements, and scored tables are
under `oarfish-evaluation-data/kinnex-wtc11/coverage-model-comparison-20260720`.
The reusable drivers are `scripts/run_kinnex_fig2_benchmark.py` and
`scripts/evaluate_kinnex_fig2.py`.

## Conclusion

Current `auto` does not help materially over current logistic on this data and
is not the best overall choice: no coverage wins relative-change recovery,
while logistic wins the paper's absolute-abundance metric. Abundance blending
is a useful compromise but is not yet the desired best-of-both-worlds rule.
This dataset should become a held-out calibration target for learning a
sample/transcript-specific coverage strength, with special attention to long
transcripts and paired-condition stability.
