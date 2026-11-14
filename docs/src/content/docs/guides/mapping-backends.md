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

spearfish builds its index from `--annotated`/`--novel` in a few seconds; `--index-out transcripts.sfi` also writes it to a file (with `--only-index`, without quantifying), and `--index transcripts.sfi` opens a saved one instead of building it, memory-mapped. Reusing a saved index skips reading the FASTA and building, and lowers peak memory (by about 2.6 GB for a human transcriptome).

oarfish records the reference's digests in the index it writes, so quantifying from it records the same reference signature as building in process. An index written by `spearfish index` (spearfish installed separately, e.g. `conda install -c bioconda spearfish` or `cargo install spearfish-cli`) also works, but has no such digests, so oarfish records a name/length digest instead.

### Presets

Its preset follows `--seq-tech`:

| `--seq-tech` | spearfish preset |
|--------------|------------------|
| `ont-cdna`   | `ont`  |
| `ont-drna`   | `drna` |
| `pac-bio-hifi` | `hifi` |
| `pac-bio`    | `ont` (spearfish has no PacBio CLR preset; a warning says so) |

oarfish sets the mapping options its models need: with `--model-truncation`, for instance, spans are extended to the read ends, also with the `ont-drna` and `pac-bio-hifi` presets, which otherwise leave them at the outermost seeds. That changes which reads the aligned-fraction filter (`-m`) drops; the log and `meta_info.json` (`spearfish_span_extension`) say whether spans were extended and why. A hit that scores `d` below the read's best counts `exp(-d / 2)` times as much: 2 is the scale spearfish suggests for its scores, and the one it records in its RAD files, so mapping in process and quantifying a RAD file weigh hits alike. (`--score-prob-denom` is in tenths of a spearfish score point with this backend, so its default is 20.) A read with only a chance match to the transcriptome, as reads from unannotated loci and contaminants have, is left unmapped by spearfish and counts as unmapped here.

### Oriented reads

spearfish maps both orientations of a cDNA read, and the read as given for direct RNA. When the reads were oriented to the transcript strand before quantification (Iso-Seq FLNC reads from `isoseq refine`, or cDNA reads turned by a tool that finds their primers), `--map-strand fw` maps only the orientation they are in, which takes less mapping time and changes little else. (`rc` maps only the reverse complement, and `both` maps both, also for direct RNA.)

```bash
$ oarfish --reads flnc.fastq.gz --index transcripts.sfi --mapper spearfish --seq-tech pac-bio-hifi \
    --map-strand fw -j 16 -o sample1
```

Use it only for reads that are oriented. It is not a filter: a read of the other strand is not simply left unmapped, since part of it may match some transcript in the mapped orientation (a gene on the opposite strand, a repeat), and that weak hit is then all the read has. To keep the alignments of one orientation of reads that were mapped in both, use `-d`/`--strand-filter`. `--map-strand` is for the spearfish backend only. oarfish warns when many reads stay unmapped under it; see [Recommended practices by sequencing technology](/oarfish/guides/sequencing-technologies/).

### Per-thread arenas

Mapping a read with spearfish makes a few hundred small allocations and frees nearly all of them before the next read. oarfish serves these (blocks of up to 8 KB) from memory that belongs to each mapping thread, so that a thread stays in the same few pages, instead of from the general allocator. Mapping takes a little less CPU time, for a few MB of memory a thread, and what is mapped does not change. This is on by default; `--bump-alloc=off` or `--bump-alloc=none` turns it off, and the option takes another limit than 8 KB as a value joined by `=` (`--bump-alloc=1KB`, or a number of bytes, at most 32 KB). The arenas are the [`nudgy`](https://crates.io/crates/nudgy) crate, the allocator `spearfish map` itself uses (there too on by default, from spearfish 0.2.1). They are for the spearfish backend only, and are there on 64-bit Linux and macOS; elsewhere oarfish maps without them (and says so when `--bump-alloc` was given).

### Keeping the mappings

To keep the mappings, map with spearfish itself and quantify the RAD file:

```bash
$ spearfish map transcripts.sfi sample1_reads.fq.gz --preset ont -t 16 --rad sample1.rad
$ oarfish -a sample1.rad -j 16 -o sample1
```
