# Leg 1 — designed TKSM simulation: pre-registered specification (2026-08-17)

Purpose: a conditions × replicates dataset with exact per-read truth and
known differential structure, to measure end-to-end power/FDR/calibration of
the oarfish 2 inference stack (presence model, presence-aware inferential
replicates, identifiability groups) under a downstream swish-style test —
including the demonstration that within-experiment prior sharing is
anti-conservative.

## Design (fixed before generation)

- **2 conditions (A, B) × 4 replicates**, 10,000,000 molecules per replicate.
- Base profile: `SQ2_quant.tsv` (TKSM SQ2 trained abundances; md5
  `f943b848…`), restricted to transcripts present in the RefSeq GTF.
- Technology: PacBio HiFi-like (SQ2): truncate KDE `SQ2_kde.json` (md5
  `1a07252c…`), badread `pacbio2016` error+qscore models, identity
  `99.19,99.99,2.09` (matched to real SQ2 header identities), polyA
  `normal(30,7), min 10` (pilot-tuned toward the real SQ2 detected-tail rate
  of 7.61%).
- Perturbations in B relative to A (disjoint sets, sampled with master seed
  20260817 from eligible pools):
  - **DTE**: 2,000 transcripts with tpm ≥ 2 — 1,000 up, 1,000 down;
    |log2FC| ∈ {0.5, 1, 2} in equal thirds.
  - **Presence-off**: 400 transcripts with 2 ≤ tpm ≤ 50 set to 0 in B.
  - **Presence-on**: 400 GTF transcripts absent from the base profile,
    switched on in B with tpm resampled from the presence-off set's values.
  - **DTU**: 400 genes with ≥ 2 isoforms at tpm ≥ 2 (untouched by the sets
    above): the top two isoforms' tpms are swapped in B (gene-level
    expression conserved; pure usage switch).
  - Everything else identical in expectation (the null set, for FDR).
- **Replicate variability**: per-transcript, per-replicate multiplicative
  lognormal noise, σ_log = 0.2 (same process in both conditions), seeded.
- Seeds: sample seed = 20260817 + sample index (A1..A4 = 1..4, B1..B4 =
  5..8); module seeds = sample seed × 10 + {1 transcribe, 2 polyA,
  3 truncate, 4 shuffle, 5 sequence}.

## Generation chain (per sample)

`tksm transcribe` (RefSeq GTF + sample abundance TSV, `--use-whole-id
--non-coding --molecule-count 10000000`) → `polyA` → `truncate
--kde-model` → `shuffle` → `sequence` (GRCh38 references, badread output,
FASTQ). Per-read truth = FASTQ `molecule_id=` comment joined to the
transcribe MDF `tid=` field (the recovery chain validated on the real TKSM
samples). Intermediate MDFs are deleted after truth extraction; retained per
sample: `fastq.gz`, `read_truth.tsv.gz`, abundance TSV.

## Pilot gate (before full generation)

One sample (A1) at 1M molecules must satisfy: read-identity mean within 1pp
of 99.19; detected-tail rate within ±3pp of 7.61% (else adjust polyA
parameters only); read-length distribution qualitatively matching real SQ2.

## Pre-registered evaluation

1. Per-sample accuracy: full stack vs baseline (Spearman/MARD vs per-sample
   truth), as in Panel B.
2. Inference: swish-style test (median-ratio scaled, infRV-moderated) over
   presence-aware inferential replicates (30/sample), A vs B:
   - power and FDR at nominal 0.05 for DTE / presence-switch / DTU sets;
   - the same with plain (presence-off) replicates — expected: worse
     calibration at presence-boundary transcripts;
   - the same with within-condition presence-prior sharing — expected:
     anti-conservative (inflated FDR), the post-selection demonstration.
3. Calibration: replicate-interval coverage of truth; infRV vs
   cross-replicate variance.
4. Identifiability: group-level (sidecar components) vs transcript-level
   testing on the unidentifiable subset.

## Storage

`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/leg1-dtu-sim/`
(regenerable from this directory's scripts + seeds + models; ~100–150 GB).
