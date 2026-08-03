# Plan: re-benchmark the `auto` coverage stack against `logistic` (2026-08-03)

## Why

Every feature in the `auto` coverage stack was selected against `none`,
`endpoint`, `hybrid`, and `adaptive`. `logistic` — the shipping default and the
model used in the paper — was benchmarked once
(`docs/coverage-evaluation-2026-07-19.md`, where it won decisively) and then
dropped from every subsequent comparison table.

An independent evaluation on NanoSim/TKSM simulations then found that `auto`
regressed against `logistic` by 0.04-0.12 Spearman
(`docs/coverage-auto-defaults-review-2026-08-03.md`). Two features were demoted
to opt-in as a result. Two others (`--alignment-calibration`,
`--censoring-model`) were left enabled because they are within noise on that
data — but they have still never been compared against `logistic` either.

This plan runs the comparison that was never run.

## Already done (branch `fix/demote-rank-blend-and-dominance-pruning`)

- `--rank-blend` default `auto` -> `none`
- `--candidate-pruning` default `auto` -> `none`
- Both remain explicitly reachable, so the original panel is reproducible
- Regression test pinning both defaults and their opt-in reachability
- `CHANGELOG.md` and `docs/coverage-auto-defaults-review-2026-08-03.md`

Effect on NanoSim NA12878 1D-cDNA: Spearman 0.7677 -> 0.8824 (vs 0.8850 for
main's `--model-coverage`), RMSE 97.28 -> 70.72 (vs 72.09), wall 2:06 -> 1:54
(vs 3:24). Parity on ranks, better RMSE, 1.8x faster.

**These defaults are provisional and this benchmark decides whether they stand.**

## Two panels

**Panel A — original selection panel.** Already on the benchmark machine under
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/`. Manifests are in
the repo: `benchmarks/coverage_ablation_manifest.tsv` (29 rows),
`benchmarks/coverage_outer_holdouts.tsv` (17), `benchmarks/longbench_coverage_panel.tsv` (30),
`benchmarks/pacbio_250k_manifest.tsv`, `benchmarks/coverage_default_public_manifest.tsv`.
Truth here is mostly *proxy*: matched Illumina quantification (the outer-validation
doc calls it "a deliberately imperfect external comparator, not molecular truth")
and SIRV controls of 69-171 transcripts.

**Panel B — simulations with exact truth.** Currently only on
`nexuscbcb01.umiacs.umd.edu`. NanoSim ONT (4 datasets) and TKSM PacBio (2).
These are the only datasets in the whole evaluation with exact read-level ground
truth at realistic scale and transcriptome size (177,816 identifiers). Copy
instructions below.

Panel B must become a permanent part of the selection panel, not a
post-hoc check.

## Copying Panel B to the benchmark machine

Both machines share a network; `nexuscbcb01.umiacs.umd.edu` has the NFS mount,
so the source paths below are literal. Run **from the benchmark machine**.

Only the transcriptome BAMs and truth files are needed — approximately **88 GB**.
Do not copy the `_sorted.bam` variants (another 36 GB, unused: oarfish consumes
the name-grouped `_T.bam`), the `.mmi` indexes, or the FASTQs. Read mode is not
needed for this benchmark: the read-mode mapper is byte-identical between the
branches, confirmed by `--coverage-model none` producing identical metrics on
both in read mode.

```bash
SRC=nexuscbcb01.umiacs.umd.edu
DEST=/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/sim-panel   # adjust
mkdir -p "$DEST"/{nanosim,tksm}

# --- NanoSim ONT: 4 transcriptome BAMs, ~60 GB ---
for d in H9_1DcDNA:cdna H9_directRNA:drna NA12878_1DcDNA:cdna NA12878_directRNA:drna; do
  rep=${d%%:*}; samp=${d##*:}
  mkdir -p "$DEST/nanosim/$rep"
  rsync -avP \
    "$SRC:/fs/cbcb-lab/rob/students/zahra/nanosim_data/alignment/$rep/${samp}_T.bam" \
    "$DEST/nanosim/$rep/"
done

# --- NanoSim ground truth, ~2.6 MB ---
rsync -avP "$SRC:/fs/cbcb-lab/rob/students/zahra/nanosim_data/ground_truth/" \
  "$DEST/nanosim/ground_truth/"

# --- TKSM PacBio: 2 transcriptome BAMs, ~28 GB ---
# One call each: brace expansion inside a quoted remote path relies on remote
# shell expansion and silently does nothing on some setups.
for s in RSII SQ2; do
  rsync -avP \
    "$SRC:/fs/cbcb-lab/rob/students/zahra/tksm_new_model/alignment/${s}_T.bam" \
    "$DEST/tksm/"
done

# --- TKSM ground truth, ~2 MB ---
rsync -avP "$SRC:/fs/cbcb-lab/rob/students/zahra/tksm_new_model/ground_truth/" \
  "$DEST/tksm/ground_truth/"
```

Dataset-to-truth mapping (ground-truth files are headerless TSV,
`transcript_id<TAB>count`, unversioned accessions; oarfish `.quant` uses
versioned `tname`, so strip `\.\d+$` before joining, then outer-join and fill 0):

| BAM | truth file | `--seq-tech` | extra flags |
|---|---|---|---|
| `nanosim/H9_1DcDNA/cdna_T.bam` | `H9_1DcDNA_ground_truth.csv` | `ont-cdna` | |
| `nanosim/H9_directRNA/drna_T.bam` | `H9_directRNA_ground_truth.csv` | `ont-drna` | `-d fw` |
| `nanosim/NA12878_1DcDNA/cdna_T.bam` | `cdna_ground_truth.csv` | `ont-cdna` | |
| `nanosim/NA12878_directRNA/drna_T.bam` | `drna_ground_truth.csv` | `ont-drna` | `-d fw` |
| `tksm/RSII_T.bam` | `RSII_ground_truth.csv` | `pac-bio` | |
| `tksm/SQ2_T.bam` | `SQ2_ground_truth.csv` | `pac-bio-hifi` | |

Use `--filter-group no-filters` for every run, matching how both panels were
generated.

## Benchmark design

**Arms.** On every sample in both panels:

| Arm | Flags |
|---|---|
| `none` | (no coverage flags) |
| **`logistic`** | `--model-coverage` |
| `adaptive-bare` | `--coverage-model auto --alignment-calibration none --candidate-pruning none --censoring-model none --rank-blend none` |
| `auto-new` | `--coverage-model auto` (new defaults: calibration + censoring on) |
| `auto-old` | `--coverage-model auto --rank-blend auto --candidate-pruning auto` |
| `+calib` | `adaptive-bare` + `--alignment-calibration agreement` |
| `+censor` | `adaptive-bare` + `--censoring-model adaptive` |
| `+prune` | `adaptive-bare` + `--candidate-pruning auto` |
| `+rank` | `adaptive-bare` + `--rank-blend auto` |

`logistic` is the arm that matters. `auto-old` reproduces the pre-change
behavior so the regression is visible on Panel A, where it was never measured.
The four `+X` arms give per-feature attribution on both panels.

**Metrics.** Spearman, Kendall, Pearson on log1p, CCC on log1p, RMSE, MARD, plus
wall time and peak RSS from `/usr/bin/time -v`. Report per sample; do not report
only panel means — the regression here is dataset-family-specific and a mean
over 28 cases hides it.

**Cost control.** 9 arms x ~40 Panel A samples + 9 x 6 Panel B samples. Panel A
samples are 50k-250k read prefixes and cheap. Panel B samples are 11-25 GB BAMs
at roughly 2-8 minutes per run. Run Panel B serially on an otherwise idle node;
timing comparisons are only valid if nothing else is competing for cores.

## Decision rules

Fix these before looking at results.

1. **A feature is an `auto` default only if it beats `logistic`** — not `none` —
   on the union of both panels, with no dataset family regressing more than
   0.01 Spearman.
2. **Rank blending and dominance pruning stay demoted** unless they clear rule 1.
   If either does, its activation gate must also be re-derived: the current 37 nt
   censoring-scale threshold was calibrated on the premise that independent
   simulations fit the 25 nt clamp, and NanoSim 1D-cDNA fits 87.85 nt.
3. **Alignment calibration and censoring are on probation.** They survived the
   demotion because they are within noise on Panel B (-0.0004, -0.0019). If they
   do not beat `logistic` on Panel A either, they come out too.
4. **If `adaptive-bare` does not beat `logistic` on the union**, `auto` should
   not be a default at all; ship `--model-coverage` and keep `auto` opt-in.
5. **If Panel A and Panel B disagree**, prefer Panel B for accuracy claims — it
   has exact truth, and Panel A's primary comparator is matched Illumina, which
   is a cross-protocol proxy. Report the disagreement explicitly rather than
   averaging it away.

## Independent of the benchmark

Fix regardless of what the numbers say:

- **Rank blending blends toward a non-converged EM.** The warm-up terminates at
  the `--coverage-warmup-iterations` cap (100) with `"converged": false`. If the
  feature survives at all, the warm-up must converge first.
- **The PacBio `auto` path is the one real cost regression.** In `bulk.rs` it
  clones `store.coverage_probabilities` (`Vec<f64>`, 0.59 GB at TKSM SQ2's
  73.5 M alignments) and holds it live across a full, unaccelerated extra EM
  (`max_iter: args.max_em_iter`, `EmAccel::None`). This is the only run in the
  whole evaluation where dev is slower than main (+4%) and carries the largest
  memory delta (+15%). Recompute rather than retain the copy, and use the same
  capped/accelerated warm-up the non-PacBio path uses.
- **`auto` should log which sub-features it activated, at INFO.** A single flag
  silently enabling four independent behaviors is why localizing this took a
  full ablation.
- **Panel A lives on node-local scratch.** `/scratch1/rob/long-read-ecosystem/`
  is not shared storage and is not backed up. Move it somewhere durable, or the
  next person to ask "does this still hold?" cannot answer.

---

# Agent prompt

Give this verbatim to an agent on the benchmarking machine.

```
You are running a decisive accuracy benchmark for the oarfish `--coverage-model auto`
stack. Read this whole prompt before starting.

## Background

oarfish's `auto` coverage stack was developed against the baselines `none`,
`endpoint`, `hybrid`, and `adaptive`. The `logistic` model -- which is what
oarfish actually ships as `--model-coverage`, and what the paper used -- was
benchmarked once on 2026-07-19 (where it won decisively) and then dropped from
every subsequent comparison table.

An independent evaluation on NanoSim/TKSM simulations with exact read-level
ground truth then found `auto` regressed against `logistic` by 0.04-0.12
Spearman, with MARD up to 2.7x worse. Add-one ablation attributed it to two
features that `auto` enabled implicitly: rank blending (-0.0967 Spearman) and
dominance pruning (-0.0282). Alignment calibration (-0.0004) and censoring
(-0.0019) were within noise.

Those two features have already been demoted to opt-in on branch
`fix/demote-rank-blend-and-dominance-pruning`. Read
`docs/coverage-auto-defaults-review-2026-08-03.md` and
`docs/coverage-auto-rebenchmark-plan-2026-08-03.md` on that branch first --
they contain the full evidence, the exact numbers, and the decision rules.

Your job is to run the comparison that was never run: every arm, including
`logistic`, on BOTH the original selection panel and the simulation panel.
The demoted defaults are provisional and your results decide whether they stand.

## Setup

1. Clone/fetch COMBINE-lab/oarfish, check out
   `fix/demote-rank-blend-and-dominance-pruning`, and build with
   `cargo build --release`. Also build `origin/main` into a separate
   CARGO_TARGET_DIR -- you need it to confirm `--model-coverage` reproduces
   across branches (it should; dev's logistic path was verified unchanged).
2. Panel A (original selection panel) should already be on this machine at
   /scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/. Verify it exists
   before planning anything. Manifests are in the repo under benchmarks/*.tsv.
   If it is missing, STOP and report -- it is node-local scratch, not shared
   storage, and may have been purged. Do not silently substitute other data.
3. Panel B (simulations) must be copied from nexuscbcb01.umiacs.umd.edu. The
   exact rsync commands, the ~88 GB size, the dataset-to-truth mapping, and the
   per-dataset --seq-tech flags are in the "Copying Panel B" section of
   coverage-auto-rebenchmark-plan-2026-08-03.md. Only transcriptome BAMs and
   truth files are needed -- do NOT copy the _sorted.bam variants, the .mmi
   indexes, or the FASTQs. Put everything on local scratch, never on an NFS
   mount.

## What to run

Nine arms, defined in the "Benchmark design" section of the plan doc, on every
sample in both panels. The `logistic` arm (`--model-coverage`) is the one that
matters -- do not omit it. `--filter-group no-filters` on every run.

Metrics: Spearman, Kendall, Pearson on log1p, CCC on log1p, RMSE, MARD, plus
wall time and peak RSS via `/usr/bin/time -v`. Ground truth joins need
accession-version stripping (`\.\d+$`) before an outer join, filling missing
with 0 -- see the plan doc.

Report per-sample results. Do NOT report only panel means: the known regression
is dataset-family-specific and a mean over 28 cases hides it entirely. That is
precisely how this was missed the first time.

Run Panel B serially on an otherwise idle node; its timing numbers are only
valid without core contention.

## Decision rules

The plan doc fixes five decision rules before results are seen. Apply them as
written. Do not relax a threshold after seeing a number. In particular: a
feature earns an `auto` default only by beating `logistic` -- beating `none` is
not evidence of anything, and is the specific error that produced this
situation.

If Panel A and Panel B disagree, say so explicitly and prefer Panel B for
accuracy claims (exact simulated truth vs. cross-protocol matched-Illumina
proxy). Do not average the disagreement away.

## Deliverable

1. A machine-readable per-sample TSV of every arm x sample x metric.
2. A markdown summary in docs/, dated, following the style of the existing
   docs/coverage-*.md evaluations, with a comparison table that includes a
   `logistic` column for every panel.
3. An explicit recommendation for each of the four features (rank blending,
   dominance pruning, alignment calibration, censoring): default-on, opt-in, or
   remove -- each justified against the `logistic` baseline and tied to the
   numbered decision rule it satisfies.
4. If any result contradicts the demotion already committed, say so plainly and
   recommend reverting it. The demotion is a hypothesis, not a conclusion to
   defend.

Report honestly: if an arm fails to run, a panel is missing, or a result is
ambiguous, say that rather than filling the gap with a plausible number.
```
