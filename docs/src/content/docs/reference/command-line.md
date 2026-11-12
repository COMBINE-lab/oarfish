---
title: Command-line options
description: oarfish's options, filter groups and the strand filter.
---

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
  -a, --alignments <ALIGNMENTS>  path to a transcriptome-aligned BAM, or a RAD file from `spearfish map --rad`, to quantify

raw read mode:
      --reads <READS>          input reads: FASTA/Q (optionally gzipped) or uBAM; format inferred from suffix
      --annotated <ANNOTATED>  transcriptome FASTA (e.g. GENCODE) to map against
      --novel <NOVEL>          additional novel/assembled transcripts (indexed together, tracked separately)
      --index <INDEX>          an existing index (oarfish-created preferred, or a compatible prebuilt index)
      --mapper <MAPPER>        read-mode mapping backend: rammap (default) or spearfish
      --map-strand <MAP_STRAND>
                               with --mapper spearfish, the read orientation to map (fw, rc or both);
                               fw for reads already oriented to the transcript strand
      --bump-alloc[=<SIZE>]    with --mapper spearfish, serve the mapper's small allocations from
                               per-thread arenas (a little less mapping time); SIZE is the
                               largest one served (8KB without it, at most 32KB), off for none
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
                               denominator D in exp((score-best)/D) weighting (default 3; with
                               --mapper spearfish 20, in tenths of a spearfish score point).
                               Transcriptome mode only
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
      --coverage-model <MODEL>     none | logistic | truncation (where reads land on transcripts) [default: none]
      --model-coverage             same as --coverage-model logistic
      --model-truncation           same as --coverage-model truncation
  -k, --growth-rate <GROWTH_RATE>  logistic `k` [default: 2]
  -b, --bin-width <BIN_WIDTH>      coverage bin width [default: 100]

output read-txps probabilities:
      --write-assignment-probs[=<WRITE_ASSIGNMENT_PROBS>]  write per-read assignment probs (uncompressed|compressed)
      --display-thresh <DISPLAY_THRESH>  min posterior prob to write to .prob, or 'none' [default: 0.000001]

EM:
      --max-em-iter <MAX_EM_ITER>            max EM iterations [default: 1000]
      --convergence-thresh <CONVERGENCE_THRESH>  EM convergence threshold [default: 0.001]
      --per-read-em                          per-read EM even without a per-read model term
      --presence-model <PRESENCE_MODEL>      none | spike-slab | bayes (experimental) [default: none]
  -q, --short-quant <SHORT_QUANT>            short-read quantification (if provided)
```

:::note
The block above is abridged for readability; run `oarfish --help` for the exact, complete option text for your installed version.
:::

## Filters

Currently, the set of filters applied in `oarfish` are directly derived from the [`NanoCount`](https://github.com/a-slide/NanoCount) tool; both the filters that exist, and the way their values are set (with the exception of the `--three-prime-clip` filter, which is not set by default in `oarfish` but is in `NanoCount`).

**In general**, if you apply a `filter-group`, the group options will be applied first and then any explicitly provided options given will override the corresponding option in the `filter-group`.

### Strand filter

The parameters above should be explained by their relevant help option, but the `-d`/`--strand-filter` is worth noting explicitly. By default, alignments to both strands of a transcript will be considered valid. You can use this option to allow only alignments in the specified orientation; for example `-d fw` will allow only alignments in the forward orientation and `-d rc` will allow only alignments in the reverse-complement orientation and `-d both` (the default) will allow both. The `-d` filter, if explicitly provided, overrides the orientation filter in any provided "filter group" so e.g. passing `--filter-group no-filters -d fw` will disable other filters, but will still only admit alignments in the forward orientation.
