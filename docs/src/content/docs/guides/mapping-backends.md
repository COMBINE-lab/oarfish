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

#### Gene labels

spearfish keeps a few candidate places for each read per gene, and collapses isoforms by gene in
its index. Build the index with gene labels: without them every transcript is its own gene, so
isoforms compete for those places one by one, ties among them go by their order in the FASTA, and
accuracy suffers. The labels can come from:

- **the FASTA headers**: GENCODE-style IDs, `transcript|gene|...`, give each transcript its gene;
- **a transcript-to-gene table**: `--t2g FILE` (whitespace-separated `transcript gene` per line;
  further columns are ignored) replaces the genes in the headers, for other FASTAs;
- **gene inference**: `--infer-genes` groups transcripts that share sequence, for a transcriptome
  without a gene annotation (spearfish's gene inference, with its default settings).

`--t2g` and `--infer-genes` apply when oarfish builds the index from `--annotated`/`--novel` with
`--mapper spearfish`; for `spearfish index`, use its own `--t2g` and `--infer-genes`. An index
without gene groups gets an `index_without_gene_groups` warning, in the log (also with `--quiet`)
and in `meta_info.json`, and every spearfish run records where the index's genes came from under
`spearfish_settings.index_genes` (`fasta`, `t2g`, `inferred` or `none`; see [the output
reference](/oarfish/reference/output/#spearfish-settings-spearfish_settings)). For a RAD file,
whose header does not record the index's genes in spearfish 0.3.0, it is `unknown` and there is no
warning.

```bash
$ oarfish --only-index --mapper spearfish --annotated transcripts.fa --t2g t2g.tsv \
    --index-out transcripts.sfi --seq-tech ont-cdna
```

oarfish records the reference's digests in the index it writes, so quantifying from it records the same reference signature as building in process. An index written by `spearfish index` (spearfish installed separately, e.g. `conda install -c bioconda spearfish` or `cargo install spearfish-cli`) also works, but has no such digests, so oarfish records a name/length digest instead.

### Presets

Its preset follows `--seq-tech`:

| `--seq-tech` | spearfish preset |
|--------------|------------------|
| `ont-cdna`   | `ont`  |
| `ont-drna`   | `drna` |
| `pac-bio-hifi` | `hifi` |
| `pac-bio`    | `ont` (spearfish has no PacBio CLR preset; a warning says so) |

oarfish sets the mapping options its models need: with `--model-truncation`, for instance, spans are extended to the read ends, also with the `ont-drna` and `pac-bio-hifi` presets, which otherwise leave them at the outermost seeds. That changes which reads the aligned-fraction filter (`-m`) drops; the log and `meta_info.json` (`spearfish_settings`: `extend_ends` and `extend_ends_reason`) say whether spans were extended and why. A hit that scores `d` below the read's best counts `exp(-d / 2)` times as much: 2 is the scale spearfish suggests for its scores, and the one it records in its RAD files, so mapping in process and quantifying a RAD file weigh hits alike. (`--score-prob-denom` is in tenths of a spearfish score point with this backend, so its default is 20.) A read with only a chance match to the transcriptome, as reads from unannotated loci and contaminants have, is left unmapped by spearfish and counts as unmapped here.

### Oriented reads

spearfish maps both orientations of a cDNA read, and the read as given for direct RNA. When the reads were oriented to the transcript strand before quantification (Iso-Seq FLNC reads from `isoseq refine`, or cDNA reads turned by a tool that finds their primers), `--map-strand fw` maps only the orientation they are in, which takes less mapping time and changes little else. (`rc` maps only the reverse complement, and `both` maps both, also for direct RNA.)

```bash
$ oarfish --reads flnc.fastq.gz --index transcripts.sfi --mapper spearfish --seq-tech pac-bio-hifi \
    --map-strand fw -j 16 -o sample1
```

Use it only for reads that are oriented. It is not a filter: a read of the other strand is not simply left unmapped, since part of it may match some transcript in the mapped orientation (a gene on the opposite strand, a repeat), and that weak hit is then all the read has. To keep the alignments of one orientation of reads that were mapped in both, use `-d`/`--strand-filter`. `--map-strand` is for the spearfish backend only. oarfish warns when many reads stay unmapped under it; see [Recommended practices by sequencing technology](/oarfish/guides/sequencing-technologies/).

### Per-thread arenas

Mapping a read with spearfish makes a few hundred small allocations and frees nearly all of them before the next read. oarfish serves these (blocks of up to 8 KB) from memory that belongs to each mapping thread, so that a thread stays in the same few pages, instead of from the general allocator. Mapping takes a little less CPU time, for a few MB of memory a thread, and what is mapped does not change. This is on by default; `--bump-alloc=off` or `--bump-alloc=none` turns it off, and the option takes another limit than 8 KB as a value joined by `=` (`--bump-alloc=1KB`, or a number of bytes, at most 32 KB). The arenas are the [`nudgy`](https://crates.io/crates/nudgy) crate, the allocator `spearfish map` itself uses (there too on by default, from spearfish 0.2.1). They are for the spearfish backend only, and are there on 64-bit Linux and macOS; elsewhere oarfish maps without them (and says so when `--bump-alloc` was given).

### Mapping and quantification scores

oarfish weighs each alignment of a read by `exp(-delta / D)`, where `delta` is how far the
alignment's score lies below the read's best and `D` is the score denominator
(`--score-prob-denom`).

With spearfish hits, oarfish uses the **quantification score**:

- **In process** (`--mapper spearfish`): the hit's `quant_score`.
- **From a RAD file** (`spearfish map --rad`): the `qscore` alignment tag when the file has it,
  otherwise the `score` tag. Files written by earlier spearfish versions have no `qscore` and
  are quantified as before.

The score filters (`--score-threshold`) see the same score, so they and the weights agree.

Which hits a read has, and spearfish's own primary hit and best set, come from the mapping
score; oarfish does not change them. The difference is only in how much weight each hit gets.
For a read with a poly(A) tail, or with a single candidate, the two scores are the same.

**When it matters.** With poly(A)-trimmed reads (or reads that lost their tail), the mapping
score has no 3'-end or strand evidence. Isoforms that differ only in their 3' end, and antisense
transcripts overlapping on opposite strands, then tie, and the EM tends to give such reads to
whichever transcript has more unique reads. The quantification score tilts these reads toward
the candidates whose 3' end the read reaches. In our tests on simulated ONT cDNA and real
spike-in data, this recovered much of the accuracy that trimming the poly(A) tails costs,
without changing spearfish's mapping output. It cannot separate antisense transcripts whose
reads fit the 3' ends on both strands (for example, a short antisense transcript lying wholly
inside the 3' end of the other); those still tie.

**Orientation evidence (spearfish 0.3.0).** The quantification score also carries evidence about
a read's orientation that spearfish 0.3.0 keeps out of its mapping score by default (its
`orient_quant` setting): a poly(A) tail followed by a primer and adapter, and, for best hits that
are one sequence on opposite strands, exact alignments of the whole read. It moves reads toward a
short isoform that lies antisense inside a longer one. This acts whether or not the end fit is on.
`--orient-evidence` chooses what the evidence does (see below); `--orient-evidence off` maps as
spearfish 0.2.1 did.

#### Options (in process)

With `--mapper spearfish`, oarfish maps with the preset `--seq-tech` names, end fit and orientation
evidence included. Three options change the end fit, and one the orientation evidence; each
defaults to the preset's value, so results do not change unless one is given:

| option | what it sets | default |
|---|---|---|
| `--end-fit SCORE` | points for the hits whose annotated 3' end the read reaches; `0` turns the end fit off (the orientation evidence still acts) | the preset's (5) |
| `--end-fit-tol BASES` | how far, in homopolymer-compressed bases, where a read ends may lie from a 3' end and still fit it | the preset's (4) |
| `--end-fit-from chain\|alignment` | how where a read ends is measured: from the chain (no alignment is run) or by aligning the read's right end | `chain` for `ont-drna`, `alignment` otherwise; `alignment` with `--model-truncation` |
| `--orient-evidence quant\|mapping\|off` | what the orientation evidence does: moves only the quantification score (`quant`), moves the mapping score and decides the orientation (`mapping`), or nothing (`off`, as spearfish 0.2.1 mapped) | the preset's (`quant`) |

The end-fit options are spearfish's `--end-fit`, `--end-fit-tol` and `--end-fit-chain` /
`--no-end-fit-chain`. `--orient-evidence mapping` is spearfish's `--no-orient-quant`, and
`--orient-evidence off` is spearfish's `--no-orient-same-exact --orient-tail-after=80`. All of
them need `--mapper spearfish` and `--reads`. For alignments read from a RAD file, the scores were
computed when the file was written, so set these on `spearfish map` instead.

oarfish has one switch for the orientation evidence rather than spearfish's separate settings:
`off` is the one combination that reproduces spearfish 0.2.1, and a separate `--no-orient-quant`
would be easy to mistake for it, since in spearfish it sends the evidence to the mapping score
instead of turning it off.

**With `--model-truncation`.** The truncation model needs where reads start and end, so oarfish
extends spans to the read ends in every preset, `ont-drna` included. The end fit then measures
every candidate by alignment: the alignments exist anyway, and the chain measurement is only
there to avoid running them for the end fit alone. `--end-fit-from chain` keeps the chain for the
candidates whose spans were not extended.

**Recorded settings.** `meta_info.json` records the spearfish settings under
`spearfish_settings`: `preset`, `extend_ends` (and why), `strands`, `quant_scores`, `end_fit`,
`end_fit_tol`, `end_fit_from`, `orient_quant` and `orient_evidence`, as resolved from the preset and the options, or,
for a RAD file, as its header records them (see [the output reference](/oarfish/reference/output/#spearfish-settings-spearfish_settings)).

**Mapping-score weights.** `--end-fit 0` alone turns only the end fit off: alignments are then
weighed by the mapping score, except for the reads whose orientation evidence moved their
quantification score (`quant_scores` stays `true`). With `--end-fit 0 --orient-evidence off`,
alignments are weighed by the mapping score alone (`quant_scores` is `false`), as oarfish weighed
spearfish 0.2.1's alignments. For RAD files, `spearfish map --end-fit 0 --no-orient-same-exact
--orient-tail-after=80` does the same. The `-fast` presets turn both off.

### Keeping the mappings

To keep the mappings, map with spearfish itself and quantify the RAD file:

```bash
$ spearfish map transcripts.sfi sample1_reads.fq.gz --preset ont -t 16 --rad sample1.rad
$ oarfish -a sample1.rad -j 16 -o sample1
```
