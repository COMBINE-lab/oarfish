# Changelog

All notable changes to `oarfish` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> Release tooling note: `dist` (cargo-dist) extracts the section whose header
> matches the release tag and uses it as the GitHub Release body. When tagging,
> rename the `## Unreleased` heading below to the chosen version, e.g.
> `## 0.10.0 - 2026-06-06`, so the notes are picked up automatically.

## Unreleased

### Added

- Opt-in `--em-accel none|squarem|daarem` for bulk estimates and bootstraps,
  with shared convergence semantics and inference diagnostics.
- `--coverage-model none|logistic|endpoint|hybrid|adaptive|degradation|auto`, including a repaired
  historical logistic baseline, a regularized joint endpoint model, and a
  support-gated log-linear combination of both. Adaptive mode learns endpoint
  shrinkage by cross-validation, smooths neighboring endpoint cells, gates by
  support and model agreement, and caps per-read coverage Bayes factors.
- An ONT direct-RNA degradation-aware model that cross-fits a sample-level
  3'-anchored competing-risk mixture, separates technical termination from
  length-scaled degradation, and removes only the degradation likelihood ratio
  from transcript evidence.
- An `auto` coverage mode that dispatches by sequencing technology and learns
  continuous degradation and endpoint-evidence weights. A learned null makes
  ONT direct-RNA samples without degradation evidence reduce exactly to the
  adaptive model.
- A competing-risk ONT direct-RNA degradation kernel that separates intact,
  length-invariant technical truncation, and length-scaled degradation before
  correcting endpoint evidence.
- An experimental bounded `--degradation-kernel piecewise2` challenger, with
  the validated `constant` kernel remaining the default. More complex
  piecewise3 and beta challengers were evaluated and removed after failing
  accuracy/runtime gates.
- Coverage-model and EM wall-clock times plus adaptive-model diagnostics in
  `.meta_info.json` for reproducible accuracy/cost evaluation.
- A documented, gated evaluation process for future coverage-model stages.
- An eight-cell-line LongBench evaluation across ONT cDNA, ONT direct RNA, and
  PacBio Kinnex, including a selective 250,000-read depth confirmation tier and
  machine-readable accuracy/runtime results.

## 0.10.3 - 2026-07-16

Patch release improving the numeric precision of the `--write-assignment-probs`
output.

### Fixed

- **Assignment probabilities in the `.prob`/`.prob.lz4` file could be reported
  as `0.000`.** The values were always printed with 3 decimal places, but the
  default `--display-thresh` (`1e-6`) admits assignments far smaller than
  3-decimal resolution, so a probability that passed the reporting threshold
  could still round to `0.000`. The printed precision now tracks
  `--display-thresh` (`decimals = ceil(-log10(display_thresh))`, clamped to
  `[3, 9]`): the default now prints 6 decimals, and an assignment that passes
  the threshold is never shown as zero. Thresholds of `1e-3` or larger are
  unchanged (still 3 decimals), and the `.prob` file format is otherwise
  identical. See [#71](https://github.com/COMBINE-lab/oarfish/issues/71).

## 0.10.2 - 2026-07-16

Patch release fixing raw-read mapping of ONT direct-RNA data basecalled in RNA
(uracil) mode.

### Fixed

- **ONT direct-RNA reads containing `U` (uracil) mapped at a near-zero rate in
  raw-read mode (`--reads`).** Reads basecalled in RNA mode carry `U`/`u`
  instead of `T`/`t`. The `rammap` seed encoding recognizes only `A/C/G/T` (`U`
  encodes as ambiguous), so `U`-containing k-mers never matched the `T`-based
  transcriptome index and ~99% of reads produced no mapping — even though the
  same reads map normally through minimap2 or the `rammap` CLI, both of which
  normalize `U` internally. oarfish now converts `U`→`T` (and `u`→`t`) once at
  read ingestion, so direct-RNA reads map as expected. BAM-input modes
  (`--alignments`, `--genome-alignments`) were never affected. See
  [#70](https://github.com/COMBINE-lab/oarfish/issues/70).

### Changed

- The raw-read discard table now reports two additional per-read outcomes —
  **read had no mapping** (the mapper returned no alignment record) and **read
  had no valid alignment** (alignments were found but none had a positive best
  score). These reads were previously dropped without being counted, so the
  table did not reconcile to the total number of reads; the new rows make the
  cause of unmapped reads visible (and were what pinpointed #70).

## 0.10.1 - 2026-07-11

Patch release fixing a `cargo install oarfish` build failure. No behavior change.

### Fixed

- **`cargo install oarfish` failed to compile against `bramble-rs` 0.1.5.** A
  fresh install re-resolves dependencies and picked up `bramble-rs` 0.1.5, whose
  `FastaDb::from_seqs` had switched to an internal deterministic-hasher map,
  rejecting the default-`RandomState` `HashMap` oarfish passed (a type mismatch
  in `src/util/aligner.rs`). See [#68](https://github.com/COMBINE-lab/oarfish/issues/68).
  Fixed by requiring `bramble-rs` 0.1.6, which makes `from_seqs` hasher-agnostic
  and exports `SeqMap` + `FastaDb::from_seq_map`; oarfish now builds bramble's map
  type directly and hands it over with no rehash. (`bramble-rs` 0.1.5 was yanked.)

## 0.10.0 - 2026-06-06

This release replaces the read-mode aligner and adds genome-based quantification.
Most users of alignment/BAM input are unaffected; users of raw-read mapping mode
should read the **Breaking changes** below.

### Breaking changes

- **Read-mode mapper swapped from `minimap2-rs` to the pure-Rust `rammap`
  backend.** The `minimap2-rs` backend and the `rammap` Cargo feature that gated
  the two have been removed; `rammap` is now the sole raw-read mapper. Quantification
  from raw reads (`--reads`) may differ slightly from previous `minimap2`-based
  runs. BAM-input modes (`--alignments`, `--genome-alignments`) are **unaffected**
  and still accept output from minimap2, pbmm2, STAR, etc.
- **`--thread-buff-size` removed.** It only configured `minimap2`'s per-thread
  `cap_kalloc` arena. It is replaced by **`--dp-cache-cap-mb`**, which bounds the
  per-thread DP alignment-scratch buffer the rammap mapper retains (default
  128 MB; also settable via the `RAMMAP_DP_CACHE_CAP_MB` environment variable;
  `0` disables the cap).

### Added

- **Genome-based quantification modes** (projection onto annotated transcripts via
  [bramble](https://github.com/COMBINE-lab/bramble)):
  - *Genome read-projection* (`--reads --genome --annotation`): spliced-align raw
    reads to the genome and project their alignments onto the transcripts.
  - *Genome-alignment projection* (`--genome-alignments --annotation`): project an
    existing spliced **genome** BAM (from minimap2/pbmm2/STAR/…) onto the transcripts.
  - Soft-clip **rescue**, on by default, recovers discriminating sequence from
    soft-clipped read ends (`--no-rescue` to disable; `--genome-fasta` to supply the
    rescue reference in genome-BAM mode).
  - `--junctions` (BED12) to hint spliced alignment instead of deriving junctions
    from `--annotation`.
- **`--dp-cache-cap-mb`**: cap (MB) on the mapper's retained per-thread DP scratch
  buffer; bounds peak RSS at high thread counts (default 128 MB; `0` = unbounded).
- **`--score-prob-denom`**: tune the score→probability denominator used to weight a
  read's alignments in the EM (transcriptome mode only).
- Documentation: genome modes, soft-clip rescue, the rammap mapping backend, and a
  GitHub Actions workflow that deploys the docs to GitHub Pages.

### Changed

- Genome-mode peak RSS reduced (~20% at high thread counts) by capping the rammap
  per-thread DP scratch-cache; alignment output is bit-identical.
- A transcriptome `--alignments` BAM produced by an unrecognized aligner now yields
  a clear error message instead of a panic.
- `--single-cell` now explicitly requires `--alignments` and conflicts with
  `--reads`, `--genome`, and `--genome-alignments`. Single-cell quantification has
  always required a transcriptome-aligned BAM; previously combining it with a
  genome mode silently ran bulk projection instead. It now errors clearly.
- **Build no longer requires a C toolchain / clang / bindgen** (the `minimap2-sys`
  dependency is gone), simplifying installation and cross-compilation.

### Notes

- Prebuilt **minimap2 `.mmi` indices remain usable** with `--index`: rammap reads
  both its own RMMI format and legacy `.mmi` files, and oarfish recomputes the
  reference signature from the index when no oarfish digest footer is present.
- The `rammap` dependency is temporarily pinned to the published fork
  `rammap-core-temp` 1.2.0, which carries the DP-cache cap (upstream
  [jwanglab/rammap#9](https://github.com/jwanglab/rammap/pull/9)). It will be
  switched back to `rammap-core` once that change is released upstream.
