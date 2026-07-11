# Changelog

All notable changes to `oarfish` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> Release tooling note: `dist` (cargo-dist) extracts the section whose header
> matches the release tag and uses it as the GitHub Release body. When tagging,
> rename the `## Unreleased` heading below to the chosen version, e.g.
> `## 0.10.0 - 2026-06-06`, so the notes are picked up automatically.

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
