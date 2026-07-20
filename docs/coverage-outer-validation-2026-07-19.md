# Coverage-model outer validation: 2026-07-19

This evaluation freezes all adaptive parameters from the preceding development
stage. No model choice, weight, prior grid, support scale, fold count, or Bayes
factor was changed after observing these samples.

The machine-readable sample manifest is
`benchmarks/coverage_outer_holdouts.tsv`, and
`scripts/run_coverage_panel.sh` runs the five coverage modes on a reusable
name-sorted BAM while reporting coverage and EM wall times.

## Degradation series

The public [PRJEB53210 degradation experiment](https://www.ebi.ac.uk/ena/browser/view/PRJEB53210)
contains matched ONT direct-RNA degradation trajectories with measured RNA
integrity. The study reports pervasive length- and isoform-dependent effects of
degradation, so the low-RIN sample is not assumed to have identical biological
truth to its high-RIN starting point. Instead, this test asks whether a coverage
model introduces an unusually large additional distribution shift.

Fixed 50,000-read prefixes were evaluated within TS10 (RIN 9.8 versus 7.7) and
TS12 (RIN 9.9 versus 7.2/7.3). Estimates were normalized to proportions and
compared with `scripts/evaluate_pair.py`.

| Trajectory | Model | Pearson | Jensen-Shannon divergence | L1 distance |
|---|---|---:|---:|---:|
| TS10 | none | **0.9662** | **0.1415** | **0.6463** |
| TS10 | endpoint | 0.9610 | 0.1572 | 0.6892 |
| TS10 | hybrid | 0.9598 | 0.1618 | 0.6986 |
| TS10 | adaptive | 0.9542 | 0.1594 | 0.6825 |
| TS12, RIN 7.2 | none | **0.9468** | **0.1298** | **0.6313** |
| TS12, RIN 7.2 | endpoint | 0.9098 | 0.1662 | 0.7354 |
| TS12, RIN 7.2 | hybrid | 0.9147 | 0.1662 | 0.7356 |
| TS12, RIN 7.2 | adaptive | 0.9196 | 0.1538 | 0.6824 |
| TS12, RIN 7.3 | none | **0.9378** | **0.1283** | **0.6155** |
| TS12, RIN 7.3 | endpoint | 0.9022 | 0.1625 | 0.7118 |
| TS12, RIN 7.3 | hybrid | 0.9058 | 0.1651 | 0.7170 |
| TS12, RIN 7.3 | adaptive | 0.9118 | 0.1537 | 0.6698 |

Adaptive consistently moderates the shift caused by endpoint/hybrid evidence
in TS12, but it does not recover the neutral model's cross-RIN stability. TS10
is less favorable. Because degradation genuinely changes observed molecule
abundance, these distances are a robustness diagnostic rather than accuracy
scores. They nevertheless argue against enabling coverage correction by
default on low-RIN samples.

Across six degradation samples, mean coverage/EM times were:

| Model | Coverage | EM |
|---|---:|---:|
| none | <0.001 ms | 1.194 s |
| endpoint | 7.43 ms | 1.212 s |
| hybrid | 34.99 ms | 1.193 s |
| adaptive | 41.83 ms | 1.205 s |

Adaptive mean reliability was 0.604 and all samples again selected prior mass
10.

## PacBio and non-SIRV biological expression

[LongBench](https://registry.opendata.aws/longbench/) supplies matched PacBio
Kinnex, ONT, and Illumina data from lung-cancer cell lines, with Sequins and
SIRV controls. We evaluated a fixed 50,000-read prefix of H69 PacBio Kinnex
against its matched Illumina transcript estimates. This is a deliberately
imperfect external comparator, not molecular truth. Metrics below use the
251,488 transcript identifiers shared after version removal between the public
matrix and the GENCODE v47 reference.

| Model | CCC | Pearson | RMSE | Coverage time | EM time |
|---|---:|---:|---:|---:|---:|
| none | 0.5152 | 0.5680 | 3218.7 | <0.001 ms | 1.063 s |
| endpoint | **0.5772** | **0.6201** | **2871.2** | 6.09 ms | 1.097 s |
| hybrid | 0.5393 | 0.5907 | 3102.4 | 34.24 ms | 1.106 s |
| adaptive | 0.5528 | 0.6012 | 3019.0 | 39.97 ms | 1.095 s |

Endpoint and adaptive both improve agreement over no coverage on this first
PacBio holdout. Adaptive is again intermediate between neutral and endpoint,
with mean reliability 0.706 and 7,610 reads reaching the Bayes-factor cap.

## Decision

- Keep adaptive opt-in. It transfers positively to PacBio, but degradation
  remains a credible failure mode for all coverage corrections.
- Surface the adaptive diagnostics in downstream QC. A low-RIN warning cannot
  be inferred reliably from alignments alone yet, so oarfish should not silently
  disable correction without sample metadata.
- Expand LongBench to all eight cell lines and use matched Illumina as an
  external endpoint. This expansion is now complete in the
  [LongBench coverage evaluation](longbench-coverage-evaluation-2026-07-19.md);
  dedicated Sequin/SIRV truth remains future work.
- For parameter learning, use leave-one-sample-out or leave-one-cell-line-out
  outer folds. Never split reads from the same library across the outer train
  and test partitions.

Artifacts are under
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/degradation` and
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/longbench`.
