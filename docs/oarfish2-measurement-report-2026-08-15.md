# oarfish 2 — measurement-phase report (2026-08-15)

Execution of the measurement-first phases of the oarfish 2 accuracy program
(plan: `nifty-wibbling-panda`; approved scope: Foundations + Measurement, then
re-prioritize). Starting state snapshotted as `oarfish2-baseline-2026-08-15`
in COMBINE-lab/oarfish (cf415a4), zrudnick/bramble (372a1e6 — commits the
previously-uncommitted projection-API working tree), and rob-p/juncprobe (new
repo). Work branches: `oarfish2-dev` (oarfish 74a87b3+, bramble).

## Phase 0 — EM convergence

**Fix** (`src/em.rs convergence_distance`): symmetric `max(prev, curr)`
denominator (dividing by *current* made geometric decay `c = r·p` report the
constant `(1-r)/r` forever) and active-set gate raised from 1e-5 to
`MIN_ACTIVE_COUNT = 1e-2` reads. A global relative-L1 companion is logged.
New `--count-floor` flag (default 1e-5 = historical) at the post-EM zeroing
step.

**What the fix revealed**: the criterion was not the whole story. With
`--em-accel none`, the residual max-relative change is pinned by transcripts
*genuinely* decaying at ~0.65%/iteration (rel diff ≈ 6.5e-3 > 1e-3 threshold)
while global relative-L1 sits at ~5e-7 — slow geometric decay of dying
transcripts, which a per-transcript relative criterion correctly reports as
still-moving. Consequently the criterion fix alone is **byte-identical** to
the baseline at the 1000-eval cap (validated old-vs-new on Panel B: identical
Spearman to 16 digits on every completed arm) — a safe no-op.

**True convergence needs a mass-weighted criterion, not just acceleration**:
`--em-accel squarem` alone still hits the 1000-eval cap on the 177k-transcript
panel (probed to 5000 evals: still `converged:false`) — the per-transcript
L∞-at-1e-3 bar is unreachable at this scale. A secondary criterion was added:
`--convergence-l1-thresh` declares convergence when `Σ|Δcount|/Σcount`
between iterations falls below the threshold. **Default 0 (disabled) — the
shipped default remains byte-identical** — because early-stopping at 1e-6
moves Spearman by up to ±0.0018 (mixed sign), exceeding the pre-registered
1e-4 no-op gate. Whether to flip defaults is an F-phase decision under the
normal accuracy gates.

### Validation (Panel B, 6 samples × {none, logistic})

| configuration | vs capped baseline | converged | evaluations |
|---|---|---|---|
| criterion fix only (`accel none`, cap 1000) | **byte-identical** (12/12 quant md5s equal) | 0/12 | 1000 |
| + `--em-accel squarem` (cap 1000) | mean Δ +0.0009 (none) / +0.000005 (logistic); max +0.0031 | 0/12 | 1000 |
| + `--convergence-l1-thresh 1e-6` | mean Δ +0.0006 (none) / −0.0005 (logistic); max &#124;Δ&#124; 0.0018 | 9/12 | 424–1000 (median ~730) |
| accel none + L1 1e-6 (single-sample probe) | — | yes | 766 |

The residual motion at the cap is genuine slow geometric decay of dying
low-count transcripts (rel diff ≈ 6.5e-3/iter while global rel-L1 ≈ 5e-7) —
i.e. the historical `converged:false` was a truthful statement about an
unreachable criterion, not (only) the divide-by-current bug. The bug is fixed;
the semantics are now controllable.

## M1 — genome-mode truth panel expansion

All six Panel B samples' simulated-read FASTQs were located on NFS
(`nanosim_data/dataset/*/…fastq.gz`, `tksm_new_model/dataset/*_shuffled.fastq.gz`)
and aligned to GRCh38 with the exact command that produced
`eval/parity/genome.bam` (`minimap2 -ax splice -N 100 --junc-bed …`; dRNA
`-uf -k14`; SQ2 `splice:hq`): `eval/genome-panel/align_panel.sh` →
`eval/genome-panel/<sample>.genome.bam` (unsorted, query-grouped, matching
the parity BAM). This takes genome-mode exact-truth samples from 1 (a 2M-read
NA12878-cDNA subset) to 6 full samples, including both PacBio chemistries.

**TKSM new-sim generation is blocked locally**: only the quantification
snakemake and generation *outputs* (kde.json, mdf chains, FASTQs) exist on
NFS; the generation pipeline/configs do not. Faithful new replicates need the
tksm tool plus the original model configs (ask Zahra, or install vpc-ccg/tksm
and drive from the per-sample kde.json + quant.tsv). Remains the
highest-value data acquisition for the PacBio anomaly (F3).

## M2 — projection information ceiling/oracle (the decisive result)

New instrumentation, all committed and pushed:

- **bramble-rs (`oarfish2-dev`)**: opt-in diagnostics
  (`ProjectionContext::enable_diagnostics` / `take_diagnostics`) recording
  every previously-silent kill site — `ElimReason::{SegmentOutsideTranscript,
  ExonSkip, DuplicateExon, NoCigar, LowSimilarity(sim), BeyondTranscriptEnd}`
  per transcript and `StrandFailure` (a read segment overlapping no annotated
  exon) per strand pass — plus `clip_score`, `total_coverage`,
  `total_operations` now exposed on `ProjectedAlignment`. Zero-cost when off.
- **oarfish (`oarfish2-dev`)**: hidden `--projection-dump <path>` (genome-BAM
  mode) writes one row per (read × candidate) outcome with all signals +
  genome AS, plus a `.tnames` sidecar.
- **juncprobe `proj-oracle`**: census, offline EM, signal-separation table,
  truth-boost/truth-exclusion oracles. 15 s for the full 7.3M-row dump.

### Complete annotation (eval/parity/genome.bam, 71,220 txps, 1.37M reads)

Offline baseline reconciles with the real oarfish run (0.8125 vs 0.8212
Spearman, same regime — deltas below are within-harness).

| measurement | value |
|---|---|
| misassigned reads (truth among candidates) | 191,501 / 1,355,059 (14.1%) |
| truth vs winner **similarity-TIED** | **96.70%** (worse than the 73–77% transcriptome-mode AS-tie rate) |
| truth strictly better — genome AS | 0.85% |
| — aligned length | 0.52% |
| — raw similarity | 0.18% |
| — ksw2 clip rescue score | 0.10% |
| — junction misses | 0.01% |
| **ANY discarded signal favors truth** | **1.10%** |
| Oracle A (perfect candidate down-weight) | +0.0534 Spearman — but ≤ ~1% addressable ⇒ realistic ≈ +0.0006 |
| Oracle B (perfect exclusion of external-origin reads) | **+0.0003** |
| beta sweep 5→40 (weight calibration) | flat (±0.0003) |

Projection's hard filters already exclude nearly all reads that "match the
genome better outside the transcript": only 0.11% of kept reads have their
truth outside the candidate set, and of 10,236 fully-failed reads only ~256
are false kills of the true transcript.

**Conclusion (pre-registered gate: ≥ +0.005 Spearman or ≥ +0.5pp read
accuracy):** per-alignment splice-agreement weighting (Arm S) and
external-read exclusion (Arm X) are **dead under a complete annotation** —
the gate is missed by an order of magnitude. This confirms and sharpens the
transcriptome-mode information ceiling: candidates surviving projection are
indistinguishable over the read's span; residual misassignment is an
abundance-prior phenomenon, not a missing per-alignment feature.

### Incomplete annotation (20% expressed-transcript holdout) — the live lead

`make_holdout.py --fraction 0.2 --seed 20260815` (6,049 expressed transcripts
deleted), dump re-run, reads cross-tabbed by truth status:

| read's true isoform | failed: strand-failure | failed: all-eliminated | kept, all candidates junc-miss | kept clean |
|---|---|---|---|---|
| **held out** (263,582 reads) | **50.37%** | **5.00%** | 2.42% | 42.21% |
| present (1,103,146 reads) | 0.79% | 0.03% | 0.59% | 98.59% |

- The projection-failure signal identifies **55.4%** of reads from missing
  isoforms at a **0.82%** background rate (67:1 per read, before locus-level
  aggregation) — and these reads currently vanish *before any quantification
  signal exists*.
- The kept-with-junction-misses signal that `--model-unannotated-isoforms`
  relies on today covers only **2.4%** of held-out-isoform reads — this *is*
  the explanation of its measured ~1% recall.
- 42.2% of held-out-isoform reads are absorbed cleanly by sibling isoforms:
  the per-read detection ceiling is ~58%.

**Conclusion:** the projection track's value concentrates in **Arm N** —
routing projection-*failed* reads (with reason codes) into the novel-state
machinery — with a measured recall pool ~23× the current signal. Arms S and X
are cut.

## M3 — converged baselines and headroom re-measurement

Presence/absence headroom re-measured on SQUAREM logistic quants
(`p0-validation/headroom.py`):

| sample | baseline | perfect detection | perfect abundance | FP txps | FP mass % |
|---|---|---|---|---|---|
| nanosim-NA12878-cdna | 0.8871 | 0.9361 | 0.9562 | 3,076 | 0.194 |
| nanosim-NA12878-drna | 0.9201 | 0.9513 | 0.9714 | 1,674 | 0.102 |
| nanosim-H9-cdna | 0.9108 | 0.9433 | 0.9708 | 2,332 | 0.103 |
| nanosim-H9-drna | 0.9363 | 0.9574 | 0.9804 | 1,088 | 0.068 |
| tksm-RSII | 0.9388 | 0.9541 | 0.9857 | 1,301 | 0.239 |
| tksm-SQ2 | 0.9479 | 0.9625 | 0.9863 | 1,960 | 0.279 |
| **mean** | **0.9235** | **0.9508** | **0.9752** | | |

**Perfect-detection headroom = +0.0273 mean Spearman** (was +0.0261 in the
handoff) — Lead A confirmed and slightly larger under converged EM, on every
sample including both PacBio chemistries.

**`--count-floor 0.5` stopgap measured** (squarem, logistic): mean **+0.0029**
(NanoSim samples +0.0041…+0.0063; TKSM PacBio −0.0005/−0.0012 — no family
regresses beyond the 0.01 rule). Captures ~11% of the detection headroom at
zero cost.

Genome-mode pipeline verified end-to-end on the new panel: full NA12878-cdna
genome BAM quantifies at Spearman 0.9022 over its 71,220-transcript universe,
converged in 511 evaluations.

## Re-prioritized feature phases

1. **F2 presence/absence (unchanged, top priority)** — the +0.026 detection
   headroom is ~40× anything the projection signals can deliver and is
   confirmed on converged baselines (see M3, +0.0273). Spike-and-slab wrapper
   around the EM as planned; `--count-floor 0.5` stopgap measured at +0.0029
   mean (safe on every family).
2. **F1′ novel-from-failures (Arm N only; S and X cut)** — plumb
   `ReadFailure`/`ElimReason` reads into `--model-unannotated-isoforms`
   (locus by genomic-span overlap, per-reason odds), genome-BAM path parity
   included. Pre-registered target: recall 1% → ≥10% at precision ≥ 0.98 on
   the holdout harness built here; also raises the `.unexplained.tsv`
   reporting value.
3. **F3 PacBio endpoint + assigned coverage** — unchanged, gated on data:
   TKSM generation configs (external ask) or an independent PacBio simulator
   family. The new RSII/SQ2 genome BAMs double the PacBio evaluation surface
   meanwhile.
4. **EM defaults** — the shipped default remains byte-identical to baseline.
   The converged configuration (`--em-accel squarem --convergence-l1-thresh
   1e-6`) is validated and opt-in; flipping it to default moves Spearman by up
   to ±0.0018 (mixed sign, within the historical cross-branch noise floor),
   so it goes through the normal F-phase accuracy gates together with the
   presence/absence work, which touches the same loop.

## Artifacts

- `eval/genome-panel/` — 6 genome-mode truth BAMs + `align_panel.sh` + logs
- `eval/proj-oracle/` — dumps (complete + holdout20), holdout GTF/manifest,
  `proj_oracle_report.txt`, baseline quants
- `rebench-2026-08-03/p0-validation/` — old/new/squarem/floor arm results,
  `headroom.py`
- juncprobe `proj-oracle` (rob-p/juncprobe), bramble diagnostics API,
  oarfish `--projection-dump` (both `oarfish2-dev`)
