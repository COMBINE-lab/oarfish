# Changelog

All notable changes to `oarfish` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> Release tooling note: `dist` (cargo-dist) extracts the section whose header
> matches the release tag and uses it as the GitHub Release body. When tagging,
> rename the `## Unreleased

### Added

- **`--model-truncation`: a learned placement model.** Each alignment's score probability is
  multiplied by the probability of where it lands on the transcript under a mixture of
  3'-anchored, 5'-anchored, full-length and uniformly placed reads, whose weights, read-length
  laws and end-gap laws are learned in the EM; the fitted mixture is logged. It adapts to the
  protocol without settings. On the SIRV E2 direct-RNA spike-ins it raises isoform-level
  Spearman from 0.817 to 0.881, and on the over-annotated SIRV reference cuts the mass placed on
  fabricated isoforms from 9.4% to 0.5%; on TKSM SQ2 simulations Spearman goes 0.950 to 0.954;
  on LongBench (50k-read subsamples) it agrees better with Illumina TPM than the default model
  on all 8 direct-RNA samples, and on 6 of 8 cDNA samples (all 8 on expressed transcripts). Works with bootstraps and in genome mode; conflicts with `--model-coverage`.
  (Supersedes the 3'-anchored model proposed in #75, which mis-ranked isoforms whenever reads
  were not 3'-anchored.)
- **`--presence-model spike-slab|bayes` (experimental).** A per-transcript presence posterior
  that discounts low-evidence transcripts during the EM; reduces false positives at some cost to
  true low-abundance transcripts. Off by default.

### Changed

- **`--score-prob-denom` defaults to 3 (was 5).** Sharper weighting of a read's alignments toward
  the best-scoring one. On the simulations with exact truth it raises Spearman slightly (TKSM SQ2
  0.950 to 0.953 for the default model, 0.954 to 0.956 with `--model-truncation`); on LongBench it
  improved 27 of 28 samples in a held-out sweep. On SIRV E2 it raises Spearman but slightly
  worsens MARD. `--score-prob-denom 5` restores the previous behaviour exactly.
- **Tail-neutral alignment scores, for all models.** The score credit a read's poly(A) tail earns
  by aligning to an annotated A-run at a transcript's end is removed, using the aligner's match
  score (from the minimap2 `@PG` line or the rammap preset). On TKSM SQ2 this improves the
  default model (Spearman 0.948 to 0.950, false positives 868 to 678); outputs on data without
  poly(A) tails are unchanged.
- **EM engines.** Without a per-read model term, the EM runs over equivalence classes (CSR
  layout); per-read models use an atomic-free CSR layout (`OARFISH_READ_EM=local|atomic` selects
  the alternatives; `--per-read-em` forces the per-read EM). Estimates match the previous
  release up to summation order.
- **Reads are ordered by likeliest transcript before inference**, which speeds up every EM;
  `--write-assignment-probs` still lists reads in input order.

### Performance

- Lower peak memory: the per-alignment record shrinks from 24 to 16 bytes (an unused field
  removed), coverage probabilities are only allocated by the coverage model, and the truncation
  EM no longer copies the alignments. Default model, versus 0.10.3: 2.36 to 1.62 GB on the NanoSim
  NA12878 cDNA simulation, 3.52 to 3.13 GB on TKSM SQ2.
- The truncation EM freezes its placement parameters once they stop improving the likelihood and
  then accelerates the abundance updates (SQUAREM).
- Genome-mode quantification of the NanoSim NA12878 cDNA simulation runs in 1:04 instead of
  2:26 (wall time), with identical estimates.

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
