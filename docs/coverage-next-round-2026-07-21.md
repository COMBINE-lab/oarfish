# Coverage and inference refinement round (2026-07-21)

## Evaluation policy and recent PacBio suite

The guarded PacBio physical `auto` model at commit `30a4056` is the frozen
accuracy baseline. Candidates must report transcript accuracy, wall time, and
peak RSS and are retained only when their benefit generalizes. The primary
PacBio panel is restricted to current high-throughput data: Revio Kinnex
LongBench at 50k and 250k depths and the Revio Kinnex WTC11/SIRV time course.
The independent HiFi simulation is a mechanism control. The existing public
FLNC SIRV sample is retained only as a secondary protocol control; no further
legacy Sequel datasets will be downloaded. The registry is
`benchmarks/pacbio_recent_suite.tsv`.

PacBio's current public catalog confirms that Kinnex is the MAS-Seq-derived
high-throughput protocol and that Revio is a current supported platform. A
2024 Revio SPRQ UHRR release is available, but each FLNC BAM is 14--16 GB and
UHRR lacks direct transcript-level ground truth. Downloading it would add
replicate consistency rather than accuracy evidence, so it is not presently
worth the space or processing time.

## Continuous nested-isoform guard: rejected

The first candidate replaced the hard extreme-creation clamp with a continuous
abundance gate only on transcripts that co-occur with a candidate at least 5 kb
longer. For physical increases, the gate was `b / (b + m)`, where `b` is the
coverage-free anchor count and `m` is the existing 37-counts-per-million
midpoint. It blended 883--1,156 transcripts per 250k sample instead of clamping
zero or one.

The candidate had effectively no runtime or memory cost, but worsened Pearson,
Spearman, CCC, RMSE, and MARD in all three 250k samples relative to the guarded
default. For example, H69 CCC fell from 0.590006 to 0.589540 and RMSE rose from
2841.29 to 2846.47. The implementation was removed. This confirms that
continuous shrinkage should not be applied to every geometrically nested
candidate; the narrow extreme guard remains better calibrated.

The earlier generic candidate-relative contrast gate is also not being
reintroduced. Its completed outer ablation merely interpolated toward no
coverage and regressed all three 250k confirmations.

## EM active-set stopping: rejected

The convergence calculation was swept over active-count thresholds of 0.01,
0.1, 1, 2, 5, and 10 reads. Thresholds through one read did not reduce the
1,000 evaluations. A two-read threshold saved only 109 final H69 evaluations
without an end-to-end gain. Five- and ten-read thresholds reduced 250k wall
time to roughly 3.4--4.1 and 3.0--3.6 seconds, respectively, but consistently
worsened MARD and perturbed low-abundance ranks. The candidate was removed.

The slow mode therefore involves supported competing isoforms, not merely
numerical dust. Future convergence work should use a likelihood or assignment
mass bound rather than excluding abundance strata.

## Packed bootstrap inference: retained

Bootstrap inference now packs invariant alignment/coverage/KDE weights once,
resamples read indices into equivalence-class multiplicities, and reuses the
packed representation across replicates. Replicates remain parallelized at the
outer level, avoiding nested parallelism.

On H69 PacBio 50k with ten bootstraps and four threads, wall time fell from
4.37 to **3.24 seconds** (26%). Peak RSS changed from 494,872 to 495,092 KB
(0.04%), and the primary quantification file is byte-identical. Bootstrap
draws are intentionally random, so replicate files are compared statistically
rather than byte-for-byte. Raw results are
`oarfish-evaluation-data/bootstrap-packed-*-h69-50k-20260721*`.

## Exact weighted equivalence-class collapsing: retained for serial EM

On H69 Kinnex/Revio 50k, grouping only by transcript identifiers could remove
47.9% of read records but would be invalid because coverage evidence is
read-specific. Grouping by the exact ordered transcript IDs plus the exact
floating-point score/coverage weights can safely remove 30.7% overall and
14.2% among ambiguous reads.

Adding exact multiplicities to serial packed EM reduces the guarded physical
50k run from 2.26 to **2.09 seconds** (about 8%) with unchanged Pearson, CCC,
RMSE, and MARD. Counts differ by at most `2e-13`; the displayed Spearman change
comes only from rank tie ordering. Peak RSS rises by about 1 MB. Applying the
same hashing to parallel 250k inference made
wall time slightly worse (4.98 to 5.27 seconds), so parallel EM deliberately
keeps its uncollapsed packed shards. Raw results are under
`oarfish-evaluation-data/exact-weighted-collapse-*` and
`weighted-eq-collapse-h69-50k-20260721*`.

## Joint abundance/coverage feedback: not retained

The existing responsibility-profile ablation already tests the critical
single-sample feedback mechanism: a coverage-free abundance fit assigns
ambiguous-read responsibilities, which then retrain coverage profiles before
coverage-aware inference. It improved one Kinnex fold-change objective but
regressed mean CCC, Spearman, and RMSE for ONT cDNA, ONT dRNA, and PacBio and
regressed all three PacBio 250k confirmations. It also added about 3.15 seconds
per million-read Kinnex sample. Repeating this circular update for more outer
iterations has no empirical justification, so it is not promoted.

Any future joint model must use held-out reads or a hierarchical external
technology prior; an in-sample abundance-to-coverage feedback loop is excluded.
