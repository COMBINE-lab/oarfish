# oarfish: transcript quantification from long-read RNA-seq data

`oarfish` is a program, written in Rust (https://www.rust-lang.org/), for quantifying transcript-level expression from long-read (i.e. Oxford nanopore cDNA and direct RNA and PacBio) sequencing technologies. It handles multi-mapping reads through the use of probabilistic allocation via an expectation-maximization (EM) algorithm.

`oarfish` can consume your data in several ways: pre-computed alignments to the *transcriptome* (a BAM file), raw reads that it maps to the transcriptome for you, or — added more recently — *genome*-based input, where reads are spliced-aligned to the genome (or a spliced genome BAM is provided) and the resulting genomic alignments are projected onto the annotated transcripts before quantification. See [Quantification modes](#quantification-modes) for an overview.

It optionally employs many filters to help discard alignments that may reduce quantification accuracy. Currently, the set of filters applied in `oarfish` are directly derived from the [`NanoCount`](https://github.com/a-slide/NanoCount)[^Gleeson] tool; both the filters that exist, and the way their values are set (with the exception of the `--three-prime-clip` filter, which is not set by default in `oarfish` but is in `NanoCount`).

Additionally, `oarfish` provides options to make use of coverage profiles derived from the aligned reads to improve quantification accuracy. The use of this coverage model is enabled with the `--model-coverage` flag. You can read more about `oarfish`[^preprint] in the [preprint](https://www.biorxiv.org/content/10.1101/2024.02.28.582591v1). Please cite the preprint if you use `oarfish` in your work or analysis.

Also, please note that `oarfish` is scientific software in active development. Therefore, please check the [GitHub Release](https://github.com/COMBINE-lab/oarfish/releases) page to make sure that you are using the latest version

## Installation

`oarfish` can be installed in a variety of ways.

### Precompiled binaries

Binaries are available via [GitHub Releases](https://github.com/COMBINE-lab/oarfish/releases).

You can quickly install the latest release using the following helper script:

```sh
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/COMBINE-lab/oarfish/releases/latest/download/oarfish-installer.sh | sh
```

### Using `cargo`

If you have `cargo` installed, you can install `oarfish` directly from the source code:

```sh
cargo install oarfish
```

You can find the crate on [crates.io](https://crates.io/crates/oarfish).

### Bioconda

`oarfish` is available via [Bioconda](https://anaconda.org/bioconda/oarfish):

```sh
conda install -c bioconda oarfish
```

## Quantification modes

`oarfish` supports four ways of providing input for **bulk** quantification. All of them feed the same EM-based quantification core; they differ only in how alignments to the transcriptome are obtained. The mode is selected implicitly by which input flag you pass. (Single-cell quantification, `--single-cell`, is separate and accepts **only** a transcriptome BAM — see [Notes about single-cell mode](#notes-about-single-cell-mode).)

| Mode | Selected by | What it does |
|------|-------------|--------------|
| **Alignment** (transcriptome BAM) | `--alignments` | Quantify pre-computed alignments of reads to the **transcriptome** (a name-sorted BAM, e.g. from `minimap2`/`pbmm2`). |
| **Read** (transcriptome mapping) | `--reads` + `--annotated`/`--index` | Map raw reads to the **transcriptome** internally (with the built-in [rammap](#the-mapping-backend) mapper), then quantify. |
| **Genome read-projection** | `--reads` + `--genome` + `--annotation` | Spliced-align raw reads to the **genome**, project each read's genomic alignments onto the annotated transcripts (via [bramble](https://github.com/COMBINE-lab/bramble)), then quantify. |
| **Genome-alignment projection** | `--genome-alignments` + `--annotation` | Project an existing **spliced genome BAM** (e.g. from `minimap2 -ax splice`) onto the annotated transcripts, then quantify. |

> **Which should I use?** The transcriptome read/alignment modes are the most mature and are the right default when you have (or can build) a transcriptome reference. The genome-based modes are useful when you want to align against the genome — for example to keep one genomic alignment pass for multiple downstream analyses, to start from genome BAMs you already have, or to capture reads that fall outside annotated transcript boundaries (recovered via soft-clip [rescue](#soft-clip-rescue)). On our benchmarks, genome read-projection quantification closely matches direct transcriptome quantification when rescue is enabled (see [expectations](#expectations-for-the-genome-based-modes)).

## Basic usage

The usage can be provided by passing `-h` at the command line.

```
A fast, accurate and versatile tool for long-read transcript quantification.

Usage: oarfish [OPTIONS] <--alignments <ALIGNMENTS>|--reads <READS>|--only-index|--genome-alignments <GENOME_ALIGNMENTS>>

Options:
      --quiet              be quiet (only warnings and errors)
      --verbose            be verbose (all non-developer logging)
  -o, --output <OUTPUT>    location where output quantification file should be written
      --single-cell        input is a single-cell BAM with the `CB:z` tag on all records
  -j, --threads <THREADS>  number of cores oarfish will use [default: 3]
      --num-bootstraps <NUM_BOOTSTRAPS>
                           number of bootstrap replicates for uncertainty [default: 0]
  -h, --help               Print help
  -V, --version            Print version

alignment mode:
  -a, --alignments <ALIGNMENTS>  path to a transcriptome-aligned BAM to quantify

raw read mode:
      --reads <READS>          input reads: FASTA/Q (optionally gzipped) or uBAM; format inferred from suffix
      --annotated <ANNOTATED>  transcriptome FASTA (e.g. GENCODE) to map against
      --novel <NOVEL>          additional novel/assembled transcripts (indexed together, tracked separately)
      --index <INDEX>          an existing index (oarfish-created preferred, or a compatible prebuilt index)
      --seq-tech <SEQ_TECH>    [possible values: ont-cdna, ont-drna, pac-bio, pac-bio-hifi]
      --best-n <BEST_N>        max secondary mappings to consider [default: 100]
      --dp-cache-cap-mb <DP_CACHE_CAP_MB>
                               cap (MB) on the per-thread DP scratch buffer the mapper retains; bounds peak
                               RSS. Unset = default 128 MB / `RAMMAP_DP_CACHE_CAP_MB` env var; `0` = unbounded

genome mode:
      --genome-alignments <GENOME_ALIGNMENTS>
                               a spliced, genome-aligned BAM (name-collated) to project onto `--annotation`
      --genome <GENOME>        a genome FASTA or prebuilt genome index; with `--reads`, spliced-align then project
      --annotation <ANNOTATION>
                               transcript annotation (GTF/GFF); required in genome mode (projection target)
      --genome-fasta <GENOME_FASTA>
                               genome FASTA for soft-clip rescue (auto-taken from `--genome` when it is a FASTA)
      --no-rescue              disable bramble soft-clip rescue during projection (on by default)
      --junctions <JUNCTIONS>  optional BED12 of splice junctions to hint alignment (else derived from `--annotation`)

filters:
      --score-prob-denom <SCORE_PROB_DENOM>
                               denominator D in exp((score-best)/D) weighting (default 5). Transcriptome mode only
      --filter-group <FILTER_GROUP>     [possible values: no-filters, nanocount-filters]
  -t, --three-prime-clip <THREE_PRIME_CLIP>     max 3' end distance [default: *4294967295]
  -f, --five-prime-clip <FIVE_PRIME_CLIP>       max 5' end distance [default: *4294967295]
  -s, --score-threshold <SCORE_THRESHOLD>       min fraction of best score for a secondary aln [default: *0.95]
  -m, --min-aligned-fraction <MIN_ALIGNED_FRACTION>  min mapped fraction of a query [default: *0.5]
  -l, --min-aligned-len <MIN_ALIGNED_LEN>       min aligned nucleotides [default: *50]
  -d, --strand-filter <STRAND_FILTER>           allowed strand (fw/+, rc/-, both/.) [default: .]

indexing:
      --only-index             only build the index, do not quantify
      --index-out <INDEX_OUT>  path where the index will be written (if provided)

coverage model:
      --model-coverage             apply the coverage model
  -k, --growth-rate <GROWTH_RATE>  logistic `k` [default: 2]
  -b, --bin-width <BIN_WIDTH>      coverage bin width [default: 100]

output read-txps probabilities:
      --write-assignment-probs[=<WRITE_ASSIGNMENT_PROBS>]  write per-read assignment probs (uncompressed|compressed)
      --display-thresh <DISPLAY_THRESH>  min posterior prob to write to .prob, or 'none' [default: 0.000001]

EM:
      --max-em-iter <MAX_EM_ITER>            max EM iterations [default: 1000]
      --convergence-thresh <CONVERGENCE_THRESH>  EM convergence threshold [default: 0.001]
  -q, --short-quant <SHORT_QUANT>            short-read quantification (if provided)
```

> The block above is abridged for readability; run `oarfish --help` for the exact, complete option text for your installed version.

## Usage examples

### Alignment-mode example (transcriptome BAM)

Assume you have ONT cDNA reads in `sample1_reads.fq.gz` and a *transcriptome* reference `transcripts.fa`, and you have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/) installed:

```bash
$ minimap2 -t 16 --eqx -N 100 -ax map-ont transcripts.fa sample1_reads.fq.gz | samtools view -@4 -b -o alignments.bam
$ oarfish -j 16 -a alignments.bam -o sample1 --filter-group no-filters --model-coverage
```

### Read-mode example (oarfish maps to the transcriptome)

Here `oarfish` maps the reads to the transcriptome for you (no external aligner needed):

```bash
$ oarfish -j 16 --reads sample1_reads.fq.gz --annotated transcripts.fa --seq-tech ont-cdna -o sample1 --filter-group no-filters --model-coverage
```

If you will quantify more than one sample against the same reference, save the index that the above command builds and reuse it:

```bash
# build + quantify, writing the index out
$ oarfish -j 16 --reads sample1_reads.fq.gz --annotated transcripts.fa --index-out transcripts.oar --seq-tech ont-cdna -o sample1 --filter-group no-filters --model-coverage
# subsequent samples reuse the index
$ oarfish -j 16 --reads sample2_reads.fq.gz --index transcripts.oar --seq-tech ont-cdna -o sample2 --filter-group no-filters --model-coverage
```

### Genome read-projection example (reads → genome → transcripts)

Spliced-align the reads to the genome and project onto the annotated transcripts. You need a genome FASTA (or a prebuilt genome index) and a transcript annotation (GTF/GFF):

```bash
$ oarfish -j 16 --reads sample1_reads.fq.gz --genome genome.fa --annotation annotation.gtf \
    --seq-tech pac-bio-hifi -o sample1 --filter-group no-filters
```

Soft-clip rescue is on by default and, when `--genome` is a FASTA, the rescue reference sequence is taken from it automatically. Add `--no-rescue` to disable rescue, or `--junctions models.bed12` to supply your own splice-junction hints instead of deriving them from the annotation.

### Genome-alignment projection example (project an existing genome BAM)

If you already have a spliced genome BAM (e.g. produced with `minimap2 -ax splice`), project it onto the transcripts directly — no reads or genome FASTA required for the projection itself, though a genome FASTA is needed if you want soft-clip rescue:

```bash
# BAM must be collated by read name
$ minimap2 -t 16 --eqx -ax splice genome.fa sample1_reads.fq.gz | samtools collate -@4 -O -u - | samtools view -b -o genome.bam
$ oarfish -j 16 --genome-alignments genome.bam --annotation annotation.gtf --genome-fasta genome.fa \
    -o sample1 --filter-group no-filters
```

All of these commands produce the output files described [below](#output).

## Input to `oarfish`

### Read-based input (transcriptome)

The read-based transcriptome mode takes reference transcript sequences (`--annotated` and/or `--novel`, `FASTA` files) plus a set of reads (`--reads`) and a `--seq-tech`. `oarfish` maps the reads to the transcriptome internally with its built-in [rammap](#the-mapping-backend) mapper and quantifies the result.

The `--annotated` and `--novel` flags distinguish the *provenance* of the underlying transcripts; their provenance is tracked separately in the output. It is recommended that you provide transcripts from known reference annotations (e.g. GENCODE) via `--annotated`, and transcripts assembled from your samples via `--novel`. **Importantly**, how transcripts are split between `--annotated` and `--novel` has no effect on the final quantification (they are joined and indexed together); only the per-source sequence signature in the output differs.

Alternatively you can provide a pre-built index via `--index` (an index built from a previous `oarfish` run is preferred, since it preserves the `--annotated`/`--novel` provenance tracking). The maximum multimapping rate is controlled by `--best-n` (default 100).

#### Read-based input formats

`oarfish` accepts `FASTA`, `FASTQ`, or unaligned `BAM` (`uBAM`) reads. When you pass reads via `--reads`, the type is inferred from the file suffix: one of `.fa`, `.fasta`, `.fq`, `.fastq` (and their upper-case and `.gz` variants) is treated as (possibly compressed) `FASTA`/`FASTQ`; `.bam`/`.ubam` (any case) is treated as `uBAM`. If the format cannot be inferred from the suffix (e.g. with process substitution), it is parsed as (possibly compressed) `FASTA`/`FASTQ`.

### Alignment-based input (transcriptome BAM)

In alignment-based mode (`--alignments`), `oarfish` processes pre-computed alignments of the reads to the *transcriptome*. The input is a `BAM` file, name-sorted (the default order produced by `minimap2` is fine), with reads aligned against the transcriptome. `oarfish` relies on the `AS` tag in each record to obtain the alignment score used in probabilistic read assignment. See [Choosing aligner options](#choosing-aligner-options-for-bam-input) for recommended aligner flags.

### Genome-based input (projection onto transcripts)

The two genome-based modes both end by *projecting* genomic alignments onto the transcripts defined by `--annotation` (a GTF/GFF), using [bramble](https://github.com/COMBINE-lab/bramble), and then quantifying those projected transcriptome alignments. `--annotation` is **required** in genome mode.

* **Genome read-projection** (`--reads` + `--genome`): `oarfish` builds (or loads) a spliced genome index from `--genome` and spliced-aligns the reads to the genome, then projects. `--genome` may be a genome `FASTA` or a prebuilt genome index. The spliced-alignment preset is chosen from `--seq-tech` (`pac-bio-hifi` uses a high-quality splice preset; other technologies use the general splice preset).

* **Genome-alignment projection** (`--genome-alignments`): `oarfish` projects an existing spliced genome `BAM`. The BAM must be **collated by read name** (e.g. `samtools collate`). This lets you reuse alignments produced by `minimap2`, `pbmm2`, STAR, etc.

#### Soft-clip rescue

Projection can lose discriminating sequence at read ends that were soft-clipped during genomic alignment. **Rescue is ON by default in genome mode**: bramble re-aligns soft-clipped read ends against the transcript's neighboring exons to recover that signal, which improves isoform-level accuracy at negligible cost. Rescue needs the genome reference sequence:

* In genome read-projection mode, when `--genome` is a FASTA the sequence is taken from it automatically.
* In genome-alignment (BAM) mode, or to override the source, supply `--genome-fasta`.
* Use `--no-rescue` to disable rescue entirely.

#### Splice-junction hints

In genome read-projection mode, splice junctions are derived automatically from `--annotation` to guide spliced alignment. You may instead supply your own junctions/transcript models as a `BED12` file via `--junctions`, in which case that file is used.

#### Expectations for the genome-based modes

These modes are newer than the transcriptome modes. On our human PacBio HiFi simulation benchmark:

* **Accuracy.** With rescue enabled, genome read-projection quantification closely tracks direct transcriptome quantification (essentially the same mapped-read counts and per-transcript estimates). Disabling rescue (`--no-rescue`) measurably lowers isoform-level accuracy, so leave it on unless you have a reason not to.
* **Resources.** Genome mode is heavier than transcriptome mode: it builds/holds a spliced genome index (tens of GB for a mammalian genome) and runs projection per read. On the 80k-read human benchmark at 48 threads, peak RSS was ~24–25 GB and wall time under a minute (excluding any one-time index build). Peak memory scales with thread count; see [`--dp-cache-cap-mb`](#the-mapping-backend) and `-j` if you need to bound it.
* **Filters.** For genome projection we typically run with `--filter-group no-filters`; the NanoCount-style transcriptome filters are tuned for direct transcriptome alignments.
* **Probability weighting.** In genome mode, a read's projected alignments are weighted by bramble exonic-coverage similarity, so the transcriptome-only `--score-prob-denom` option is rejected (it is an error to pass it with `--genome`/`--genome-alignments`).

### Choosing aligner options (for BAM input)

When you provide your own BAM (`--alignments` or `--genome-alignments`), it is important that the alignments are generated in a way compatible with quantification. Primarily, the aligner should report as many optimal and near-optimal alignments as exist, so `oarfish` can see all of them and allocate reads accordingly. We recommend the following `minimap2` options:

* For transcriptome alignment (`--alignments`): ONT data (dRNA or cDNA) — `--eqx -N 100 -ax map-ont`; PacBio data — `--eqx -N 100 -ax map-pb` (or `-ax map-hifi` for very high-quality PacBio HiFi).
* For genome alignment (`--genome-alignments`): use a spliced preset, e.g. `--eqx -ax splice` (or `-ax splice:hq` for PacBio HiFi), and collate the BAM by read name.

**Note**: a larger `-N` (e.g. the [TranSigner manuscript](https://www.biorxiv.org/content/10.1101/2024.04.13.589356v1.full) recommends `-N 181`) should not diminish accuracy, but will make alignment slower and the BAM larger.

## The mapping backend

For raw-read input (`--reads`), `oarfish` no longer shells out to `minimap2`; mapping is performed in-process by [rammap](https://github.com/jwanglab/rammap), a pure-Rust long-read mapper. This removes the external alignment step and the `minimap2`/`htslib` build dependencies. Externally produced BAMs (`--alignments`, `--genome-alignments`) are of course still fully supported — those can come from `minimap2`, `pbmm2`, STAR, or any compatible aligner.

The `--seq-tech` value selects the mapping preset:

| `--seq-tech` | transcriptome preset (`--annotated`/`--index`) | genome preset (`--genome`) |
|--------------|-----------------------------------------------|----------------------------|
| `ont-cdna`   | `map-ont`  | `splice`     |
| `ont-drna`   | `map-ont`  | `splice`     |
| `pac-bio`    | `map-pb`   | `splice`     |
| `pac-bio-hifi` | `map-hifi` | `splice:hq` |

**Bounding peak memory.** The mapper retains a per-thread DP alignment-scratch buffer between reads to avoid reallocation. Rare reads with large unaligned gaps can produce very large buffers, and at high thread counts the retained per-thread high-water marks can dominate peak RSS. `oarfish` caps that buffer at 128 MB by default (frees larger one-off buffers instead of pinning them), which on our genome benchmark cut peak RSS by ~20% with bit-identical output. Tune it with `--dp-cache-cap-mb <MB>` (or the `RAMMAP_DP_CACHE_CAP_MB` environment variable); pass `0` to disable the cap and restore the previous unbounded behavior.

## Other notes on `oarfish` parameters

The parameters above should be explained by their relevant help option, but the `-d`/`--strand-filter` is worth noting explicitly. By default, alignments to both strands of a transcript will be considered valid.  You can use this option to allow only alignments in the specified orientation; for example `-d fw` will allow only alignments in the forward orientation and `-d rc` will allow only alignments in the reverse-complement orientation and `-d both` (the default) will allow both.  The `-d` filter, if explicitly provided, overrides the orientation filter in any provided "filter group" so e.g. passing `--filter-group no-filters -d fw` will disable other filters, but will still only admit alignments in the forward orientation.

**In general**, if you apply a `filter-group`, the group options will be applied first and then any explicitly provided options given will override the corresponding option in the `filter-group`.

### Read-level assignment probabilities

`oarfish` has the ability to output read-level assignment probabilities.  That is, for each input read, what is the probability, conditioned on the final estimate of transcript abundances, that the read was sequenced from each transcript to which it aligned. By default, this information is not recorded (as it's not required, or commonly used, for most standard analyses). To enable this output, you should pass the `--write-assignment-probs` option to `oarfish`.  Optionally, you may also pass `--write-assignment-probs=compressed` to write the output to a compressed ([lz4](https://github.com/lz4/lz4)) stream --- the default
output is to an uncompressed text file.

The format of the read assignment probabilities is as follows --- where all fields below on a single line are `\t` delimited:

```
<m> <n>
<tname_1>
<tname_2>
...
<tname_m>
<rname_1> <k_1> <tid_11> <tid_21> ... <tid_{k_1}1> <p_11> <p_21> ... <p_{k_1}1>
...
<rname_n> <k_1> <tid_1n> <tid_2n> ... <tid_{k_n}n> <p_1n> <p_2n> ... <p_{k_n}n>
```

Here, `<m>` is the number of transcripts in the reference, `<n>` is the number of *mapped* reads. The following `m` lines consist of the names of the transcripts in the order they will be referred to in the file. The following `n` lines after that provide the actual alignment information and assignment probabilities for each read.

The format of each of these lines is as follows; `<rname>` is the name of the read, `<k>` is the number of alignments for which probabilities are reported (if it is determined an alignment is invalid under the model, it may not be reported). Subsequently, there is a list of `k` integers, and `k` floating point numbers. Each of the `k` integers is the index of some transcript in the list of `m` transcripts given at the start of the file, and the subsequent list of `k` floating point numbers are the assignment probabilities of this read to each of the transcripts.

For example:

```
5 3
txpA
txpB
txpC
txpD
txpE
readX 2 0 2 0.4 0.6
readY 3 1 3 4 0.1 0.2 0.7
readZ 1 4 1.0
```
 
This file provides an example with 5 reference transcripts and 3 mapped reads. The first read (readX) maps to transcripts 0 and 2 (so txpA and txpC) with probabilities 0.4 and 0.6 respectively. The next read (readY) maps to 3 transcripts, txpB, txpD and txpE with probabilities 0.1, 0.2, and 0.7 respectively. Finally, the last read (readZ) maps uniquely to transcript txpE with probability 1.

The compressed output (i.e. what is generated if one passes `--write-assignment-probs=compressed`) is exactly the same format, except instead of residing in a plain text file, it is written to an lz4 compressed text file.  You can either decompress this file first with an lz4 decompressor, or decompress it on-the-fly as you are parsing the file using the lz4 library in your favorite language.

## Notes about single-cell mode

Starting with version 0.6.1 `oarfish` incorporates the first single-cell quantification capabilities. Given a `bam` file, **collated by cell barcode and with already (UMI) deduplicated reads**, this mode, enabled with the `--single-cell` flag, will allow `oarfish` to produce a single-cell quantification matrix.

> **Single-cell mode requires a *transcriptome*-aligned BAM, supplied with `--alignments`.** Unlike bulk quantification — which can take a transcriptome BAM (`--alignments`), raw reads (`--reads`), a genome FASTA + reads (`--genome`), or a genome BAM (`--genome-alignments`) — single-cell mode supports **only** the transcriptome-BAM path. It does not support raw-read mapping or either of the genome (projection) modes. The CLI enforces this: `--single-cell` requires `--alignments` and conflicts with `--reads`, `--genome`, and `--genome-alignments`. The input BAM must therefore be reads aligned to the *transcriptome* (e.g. with minimap2/pbmm2), not to the genome.

**Formatting requirements of BAM input in single-cell mode**: All alignment records for the same cell barcode should be adjacent in the `bam` file, and a count will be obtained for each read record, so UMI de-duplication should have been performed if those are the counts you want. In the future, counting UMIs directly may be supported, and some of these other restrictions may be lifted.

## Inferential Replicates

`oarfish` has the ability to compute [_inferential replicates_](https://academic.oup.com/nar/article/47/18/e105/5542870) of its quantification estimates. This is performed by bootstrap sampling of the original read mappings, and subsequently performing inference under each resampling.  These inferential replicates allow assessing the variance of the point estimate of transcript abundance, and can lead to improved differential analysis at the transcript level, if using a differential testing tool that takes advantage of this information. The generation of inferential replicates is controlled by the `--num-bootstraps` argument to `oarfish`.  The default value is `0`, meaning that no inferential replicates are generated.  If you set this to some value greater than `0`, the the requested number of inferential replicates will be generated. It is recommended, if generating inferential replicates, to run `oarfish` with multiple threads, since replicate generation is highly-parallelized. Finally, if replicates are generated, they are written to a [`Parquet`](https://parquet.apache.org/), starting with the specified output stem and ending with `infreps.pq`.

## Output

The `--output` option passed to `oarfish` corresponds to a path prefix (this prefix can contain the path separator character and if it refers to a directory that does not yeat exist, that directory will be created). Based on this path prefix, say `P`, `oarfish` will create 2 files:

  * `P.meta_info.json` - a JSON format file containing information about relevant parameters with which `oarfish` was run, and other relevant inforamtion from the processed sample apart from the actual transcript quantifications.
  * `P.quant` - a tab separated file listing the quantified targets, as well as information about their length and other metadata. The `num_reads` column provides the estimate of the number of reads originating from each target.
  * `P.infreps.pq` - a [`Parquet`](https://parquet.apache.org/) table where each row is a transcript and each column is an inferential replicate, containing the estimated counts for each transcript under each computed inferential replicate.
  * `P.ambig_info.tsv` - a tab separated file listing, for each transcript (in the same order in which they appear in `P.quant`) the number of uniquely mapped, ambiguously mapped, and total reads.  The quantification estimate for each transcript, in general, should reside between the number of uniquely aligned reads and the total number of reads (i.e. these provide, respectively lower and upper bounds for the number of reads assigned to each transcript).  Note that the total in this file is the total number of reads that align to this transcript with a sufficiently high alignment score --- it is _not_, in general, an estimate of the number of reads originating from this transcript as many of those reads can be multimapping and, in fact, potentially better described by other transcripts.
  * `P.prob[.lz4]` - a file encoding the assignment probability of each read to each transcript to which it had a valid alignment (optionally compressed using [`lz4`](https://github.com/lz4/lz4)). This file is optional and is generated only if `--write-assignment-probs` is passed to `oarfish`.

## References

[^Gleeson]: Josie Gleeson, Adrien Leger, Yair D J Prawer, Tracy A Lane, Paul J Harrison, Wilfried Haerty, Michael B Clark, Accurate expression quantification from nanopore direct RNA sequencing with NanoCount, Nucleic Acids Research, Volume 50, Issue 4, 28 February 2022, Page e19, [https://doi.org/10.1093/nar/gkab1129](https://doi.org/10.1093/nar/gkab1129)

[^preprint]: Zahra Zare Jousheghani, Rob Patro. Oarfish: Enhanced probabilistic modeling leads to improved accuracy in long read transcriptome quantification, bioRxiv 2024.02.28.582591; doi: [https://doi.org/10.1101/2024.02.28.582591](https://doi.org/10.1101/2024.02.28.582591)
