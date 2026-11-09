---
title: Input
description: Read-based, alignment-based and genome-based input, and aligner options for BAM input.
---

## Read-based input (transcriptome)

The read-based transcriptome mode takes reference transcript sequences (`--annotated` and/or `--novel`, `FASTA` files) plus a set of reads (`--reads`) and a `--seq-tech`. `oarfish` maps the reads to the transcriptome internally with its built-in [rammap](/oarfish/guides/mapping-backends/) mapper and quantifies the result.

The `--annotated` and `--novel` flags distinguish the *provenance* of the underlying transcripts; their provenance is tracked separately in the output. It is recommended that you provide transcripts from known reference annotations (e.g. GENCODE) via `--annotated`, and transcripts assembled from your samples via `--novel`. **Importantly**, how transcripts are split between `--annotated` and `--novel` has no effect on the final quantification (they are joined and indexed together); only the per-source sequence signature in the output differs.

Alternatively you can provide a pre-built index via `--index` (an index built from a previous `oarfish` run is preferred, since it preserves the `--annotated`/`--novel` provenance tracking). The maximum multimapping rate is controlled by `--best-n` (default 100).

### Read-based input formats

`oarfish` accepts `FASTA`, `FASTQ`, or unaligned `BAM` (`uBAM`) reads. When you pass reads via `--reads`, the type is inferred from the file suffix: one of `.fa`, `.fasta`, `.fq`, `.fastq` (and their upper-case and `.gz` variants) is treated as (possibly compressed) `FASTA`/`FASTQ`; `.bam`/`.ubam` (any case) is treated as `uBAM`. If the format cannot be inferred from the suffix (e.g. with process substitution), it is parsed as (possibly compressed) `FASTA`/`FASTQ`.

## Alignment-based input (transcriptome BAM)

In alignment-based mode (`--alignments`), `oarfish` processes pre-computed alignments of the reads to the *transcriptome*. The input is a `BAM` file, name-sorted (the default order produced by `minimap2` is fine), with reads aligned against the transcriptome. `oarfish` relies on the `AS` tag in each record to obtain the alignment score used in probabilistic read assignment. See [Choosing aligner options](#choosing-aligner-options-for-bam-input) for recommended aligner flags.

## Genome-based input (projection onto transcripts)

The two genome-based modes both end by *projecting* genomic alignments onto the transcripts defined by `--annotation` (a GTF/GFF), using [bramble](https://github.com/COMBINE-lab/bramble), and then quantifying those projected transcriptome alignments. `--annotation` is **required** in genome mode.

- **Genome read-projection** (`--reads` + `--genome`): `oarfish` builds (or loads) a spliced genome index from `--genome` and spliced-aligns the reads to the genome, then projects. `--genome` may be a genome `FASTA` or a prebuilt genome index. The spliced-alignment preset is chosen from `--seq-tech` (`pac-bio-hifi` uses a high-quality splice preset; other technologies use the general splice preset).
- **Genome-alignment projection** (`--genome-alignments`): `oarfish` projects an existing spliced genome `BAM`. The BAM must be **collated by read name** (e.g. `samtools collate`). This lets you reuse alignments produced by `minimap2`, `pbmm2`, STAR, etc.

### Soft-clip rescue

Projection can lose discriminating sequence at read ends that were soft-clipped during genomic alignment. **Rescue is ON by default in genome mode**: bramble re-aligns soft-clipped read ends against the transcript's neighboring exons to recover that signal, which improves isoform-level accuracy at negligible cost. Rescue needs the genome reference sequence:

- In genome read-projection mode, when `--genome` is a FASTA the sequence is taken from it automatically.
- In genome-alignment (BAM) mode, or to override the source, supply `--genome-fasta`.
- Use `--no-rescue` to disable rescue entirely.

### Splice-junction hints

In genome read-projection mode, splice junctions are derived automatically from `--annotation` to guide spliced alignment. You may instead supply your own junctions/transcript models as a `BED12` file via `--junctions`, in which case that file is used.

### Expectations for the genome-based modes

These modes are newer than the transcriptome modes. On our human PacBio HiFi simulation benchmark:

- **Accuracy.** With rescue enabled, genome read-projection quantification closely tracks direct transcriptome quantification (essentially the same mapped-read counts and per-transcript estimates). Disabling rescue (`--no-rescue`) measurably lowers isoform-level accuracy, so leave it on unless you have a reason not to.
- **Resources.** Genome mode is heavier than transcriptome mode: it builds/holds a spliced genome index (tens of GB for a mammalian genome) and runs projection per read. On the 80k-read human benchmark at 48 threads, peak RSS was ~24–25 GB and wall time under a minute (excluding any one-time index build). Peak memory scales with thread count; see [`--dp-cache-cap-mb`](/oarfish/guides/mapping-backends/#bounding-peak-memory) and `-j` if you need to bound it.
- **Filters.** For genome projection we typically run with `--filter-group no-filters`; the NanoCount-style transcriptome filters are tuned for direct transcriptome alignments.
- **Probability weighting.** In genome mode, a read's projected alignments are weighted by bramble exonic-coverage similarity, so the transcriptome-only `--score-prob-denom` option is rejected (it is an error to pass it with `--genome`/`--genome-alignments`).

## Choosing aligner options (for BAM input)

When you provide your own BAM (`--alignments` or `--genome-alignments`), it is important that the alignments are generated in a way compatible with quantification. Primarily, the aligner should report as many optimal and near-optimal alignments as exist, so `oarfish` can see all of them and allocate reads accordingly. We recommend the following `minimap2` options:

- For transcriptome alignment (`--alignments`): ONT data (dRNA or cDNA) — `--eqx -N 100 -ax map-ont`; PacBio data — `--eqx -N 100 -ax map-pb` (or `-ax map-hifi` for very high-quality PacBio HiFi).
- For genome alignment (`--genome-alignments`): use a spliced preset, e.g. `--eqx -ax splice` (or `-ax splice:hq` for PacBio HiFi), and collate the BAM by read name.

:::note
A larger `-N` (e.g. the [TranSigner manuscript](https://www.biorxiv.org/content/10.1101/2024.04.13.589356v1.full) recommends `-N 181`) should not diminish accuracy, but will make alignment slower and the BAM larger.
:::
