---
title: Mapping backends
description: How oarfish maps raw reads in process, with rammap (the default) or spearfish.
---

## rammap

For raw-read input (`--reads`), `oarfish` no longer shells out to `minimap2`; mapping is performed in-process by [rammap](https://github.com/jwanglab/rammap), a pure-Rust long-read mapper. This removes the external alignment step and the `minimap2`/`htslib` build dependencies. Externally produced BAMs (`--alignments`, `--genome-alignments`) are of course still fully supported — those can come from `minimap2`, `pbmm2`, STAR, or any compatible aligner.

The `--seq-tech` value selects the mapping preset:

| `--seq-tech` | transcriptome preset (`--annotated`/`--index`) | genome preset (`--genome`) |
|--------------|-----------------------------------------------|----------------------------|
| `ont-cdna`   | `map-ont`  | `splice`     |
| `ont-drna`   | `map-ont`  | `splice`     |
| `pac-bio`    | `map-pb`   | `splice`     |
| `pac-bio-hifi` | `map-hifi` | `splice:hq` |

### Bounding peak memory

The mapper retains a per-thread DP alignment-scratch buffer between reads to avoid reallocation. Rare reads with large unaligned gaps can produce very large buffers, and at high thread counts the retained per-thread high-water marks can dominate peak RSS. `oarfish` caps that buffer at 128 MB by default (frees larger one-off buffers instead of pinning them), which on our genome benchmark cut peak RSS by ~20% with bit-identical output. Tune it with `--dp-cache-cap-mb <MB>` (or the `RAMMAP_DP_CACHE_CAP_MB` environment variable); pass `0` to disable the cap and restore the previous unbounded behavior.

## spearfish

With `--mapper spearfish`, read mode maps with [spearfish](https://github.com/COMBINE-lab/spearfish) instead of rammap, in process: a long-read transcriptome mapper that reports, for each read, every compatible transcript with its score difference from the best. It supports transcriptome read mode only (not the genome modes).

### The spearfish index

spearfish builds its index from `--annotated`/`--novel` in a few seconds; `--index-out transcripts.sfi` also writes it to a file (with `--only-index`, without quantifying), and `--index transcripts.sfi` opens a saved one instead of building it, memory-mapped. Reusing a saved index skips reading the FASTA and building, and lowers peak memory (by about 2 GB for a human transcriptome).

oarfish records the reference's digests in the index it writes, so quantifying from it records the same reference signature as building in process. An index written by `spearfish index` (spearfish installed separately, e.g. `conda install -c bioconda spearfish` or `cargo install spearfish-cli`) also works, but has no such digests, so oarfish records a name/length digest instead.

### Presets

Its preset follows `--seq-tech`:

| `--seq-tech` | spearfish preset |
|--------------|------------------|
| `ont-cdna`   | `ont`  |
| `ont-drna`   | `drna` |
| `pac-bio-hifi` | `hifi` |
| `pac-bio`    | `ont` (spearfish has no PacBio CLR preset; a warning says so) |

oarfish sets the mapping options its models need: with `--model-truncation`, for instance, spans are extended to the read ends. Hits are weighted by their score difference from the read's best on spearfish's scale (`--score-prob-denom` defaults accordingly).

### Keeping the mappings

To keep the mappings, map with spearfish itself and quantify the RAD file:

```bash
$ spearfish map transcripts.sfi sample1_reads.fq.gz --preset ont -t 16 --rad sample1.rad
$ oarfish -a sample1.rad -j 16 -o sample1
```
