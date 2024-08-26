# oarfish: transcript quantification from long-read RNA-seq data

### Basic usage

`oarfish` is a program, written in [`rust`](https://www.rust-lang.org/), for quantifying transcript-level expression from long-read (i.e. Oxford nanopore cDNA and direct RNA and PacBio) sequencing technologies. `oarfish` requires a sample of sequencing reads aligned to the *transcriptome* (currntly not to the genome). It handles multi-mapping reads through the use of probabilistic allocation via an expectation-maximization (EM) algorithm.  

It optionally employs many filters to help discard alignments that may reduce quantification accuracy.  Currently, the set of filters applied in `oarfish` are directly derived from the [`NanoCount`](https://github.com/a-slide/NanoCount)[^Gleeson] tool; both the filters that exist, and the way their values are set (with the exception of the `--three-prime-clip` filter, which is not set by default in `oarfish` but is in `NanoCount`).

Additionally, `oarfish` provides options to make use of coverage profiles derived from the aligned reads to improve quantification accuracy.  The use of this coverage model is enabled with the `--model-coverage` flag. You can read more about `oarfish`[^preprint] in the [preprint](https://www.biorxiv.org/content/10.1101/2024.02.28.582591v1). Please cite the preprint if you use `oarfish` in your work or analysis.

Also, please note that `oarfish` is scientific software in active development.  Therefore, please check the [GitHub Release](https://github.com/COMBINE-lab/oarfish/releases) page to make sure that you are using the latest version 
(also, the `dev` branch should compile from source at all times so feel free to use it, but let us know if you run into any issues).

The usage can be provided by passing `-h` at the command line.
```
A fast, accurate and versatile tool for long-read transcript quantification.

Usage: oarfish [OPTIONS] --alignments <ALIGNMENTS> --output <OUTPUT>

Options:
      --quiet
          be quiet (i.e. don't output log messages that aren't at least warnings)
      --verbose
          be verbose (i.e. output all non-developer logging messages)
  -a, --alignments <ALIGNMENTS>
          path to the file containing the input alignments
  -o, --output <OUTPUT>
          location where output quantification file should be written
  -j, --threads <THREADS>
          maximum number of cores that the oarfish can use to obtain binomial probability [default: 1]
      --num-bootstraps <NUM_BOOTSTRAPS>
          number of bootstrap replicates to produce to assess quantification uncertainty [default: 0]
  -h, --help
          Print help
  -V, --version
          Print version

filters:
      --filter-group <FILTER_GROUP>
          [possible values: no-filters, nanocount-filters]
  -t, --three-prime-clip <THREE_PRIME_CLIP>
          maximum allowable distance of the right-most end of an alignment from the 3' transcript end [default: 4294967295]
  -f, --five-prime-clip <FIVE_PRIME_CLIP>
          maximum allowable distance of the left-most end of an alignment from the 5' transcript end [default: 4294967295]
  -s, --score-threshold <SCORE_THRESHOLD>
          fraction of the best possible alignment score that a secondary alignment must have for consideration [default: 0.95]
  -m, --min-aligned-fraction <MIN_ALIGNED_FRACTION>
          fraction of a query that must be mapped within an alignemnt to consider the alignemnt valid [default: 0.5]
  -l, --min-aligned-len <MIN_ALIGNED_LEN>
          minimum number of nucleotides in the aligned portion of a read [default: 50]
  -d, --strand-filter <STRAND_FILTER>
          only alignments to this strand will be allowed; options are (fw /+, rc/-, or both/.) [default: .]
```


## Input to `oarfish`

`Oarfish` can accept as input either a `bam` file containing reads aligned to the transcriptome as specified [below](index.md#alignment-based-input), or
raw sequencing reads themselves (along with a reference transcriptome), which are then mapped to the reference using [minimap2-rs](https://github.com/jguhlin/minimap2-rs)
and subsequently processed with `oarfish`.  With equivalent alignment options, the results of these input modes should be equivalent, so which to use is therefore
based on the preference of the user.


### Read-based input

The read-based input mode takes as input a reference (specified with the `--reference` argument), which can be either a `FASTA` file containing a transcriptome reference
or an pre-build `minimap2` index, as well as a set of reads (specified with the `--reads` argument), and a `--seq-tech` argument specifying the sequencing technology 
type of the reads to be mapped.

The mapping between the potential values that can be passed to `oarfish`'s `--seq-tech` argument and the `minimap2` presets is as follows:

  - `oarfish` seq-tech `ont-cdna` corresponds to `minimap2` preset `map-ont`
  - `oarfish` seq-tech `ont-drna` corresponds to `minimap2` preset `map-ont`
  - `oarfish` seq-tech `pac-bio` corresponds to `minimap2` preset `map-pb`
  - `oarfish` seq-tech `pac-bio-hifi` corresponds to `minimap2` preset `map-hifi`

Given these inputs, `oarfish` will either load the pre-built `minimap2` index, or build one according to the parameter specified by `--seq-tech`, and will then align
the reads to this index using [`minimap2-rs`](https://github.com/jguhlin/minimap2-rs).  Optionally, the maximum multimapping rate (i.e. the number of secondary alignments 
corresponding to the `minimap2` parameter `-N`) can be specified with the command line parameter `--best-n`. The default value of this parameter is 100.

### Alignmment-based input

In alignment-based mode, `oarfish` processes pre-computed alignments of hte read to the transcriptome. The input should be a `bam` format file, with reads aligned using [`minimap2`](https://github.com/lh3/minimap2) against the _transcriptome_. That is, `oarfish` does not currently handle spliced alignment to the genome. Further, the output alignments should be name sorted (the default order produced by `minimap2` should be fine). _Specifically_, `oarfish` relies on the existence of the `AS` tag in the `bam` records that encodes the alignment score in order to obtain the score for each alignment (which is used in probabilistic read assignment), and the score of the best alignment, overall, for each read. ### Choosing `minimap2` alignment options Since the purpose of `oarfish` is to estimate transcript abundance from a collection of alignments to the target transcriptome, it is important that the alignments are generated in a fashion that is compatible with this goal.  Primarily, this means that the aligner should be configured to report as many optimal (and near-optimal) alignments as exist, so that `oarfish` can observe all of this information and determine how to allocate reads to transcripts.  We recommend using the following options with `minimap2` when aligning data for later processing by `oarfish` * For ONT data (either dRNA or cDNA): please use the flags `--eqx -N 100 -ax map-ont` For PacBio data: please use the flags `--eqx -N 100 -ax pacbio` **Note (1)**: It may be worthwile using an even larger `N` value (e.g. the [TranSigner manuscript](https://www.biorxiv.org/content/10.1101/2024.04.13.589356v1.full) recommends `-N 181`). A larger value should not diminish the accuracy of `oarfish`, but it may make alignment take longer and produce a larger `bam` file.

**Note (2)**: For very high quality PacBio data, it may be most appropriate to use the `-ax map-hifi` flag in place of `-ax pacbio`.  We are currently evaluating the effect of this option, and also welcome feedback if you have experiences to share on the use of data aligned with these different flags with `oarfish`.

### Other notes on `oarfish` parameters

The parameters above should be explained by their relevant help option, but the `-d`/`--strand-filter` is worth noting explicitly. By default, alignments to both strands of a transcript will be considered valid.  You can use this option to allow only alignments in the specified orientation; for example `-d fw` will allow only alignments in the forward orientation and `-d rc` will allow only alignments in the reverse-complement orientation and `-d both` (the default) will allow both.  The `-d` filter, if explicitly provided, overrides the orientation filter in any provided "filter group" so e.g. passing `--filter-group no-filters -d fw` will disable other filters, but will still only admit alignments in the forward orientation.

**In general**, if you apply a `filter-group`, the group options will be applied first and then any explicitly provided options given will override the corresponding option in the `filter-group`.

### Inferential Replicates

`oarfish` has the ability to compute [_inferential replicates_](https://academic.oup.com/nar/article/47/18/e105/5542870) of its quantification estimates. This is performed by bootstrap sampling of the original read mappings, and subsequently performing inference under each resampling.  These inferential replicates allow assessing the variance of the point estimate of transcript abundance, and can lead to improved differential analysis at the transcript level, if using a differential testing tool that takes advantage of this information. The generation of inferential replicates is controlled by the `--num-bootstraps` argument to `oarfish`.  The default value is `0`, meaning that no inferential replicates are generated.  If you set this to some value greater than `0`, the the requested number of inferential replicates will be generated. It is recommended, if generating inferential replicates, to run `oarfish` with multiple threads, since replicate generation is highly-parallelized. Finally, if replicates are generated, they are written to a [`Parquet`](https://parquet.apache.org/), starting with the specified output stem and ending with `infreps.pq`.

### Output

The `--output` option passed to `oarfish` corresponds to a path prefix (this prefix can contain the path separator character and if it refers to a directory that does not yeat exist, that directory will be created). Based on this path prefix, say `P`, `oarfish` will create 2 files:

  * `P.meta_info.json` - a JSON format file containing information about relevant parameters with which `oarfish` was run, and other relevant inforamtion from the processed sample apart from the actual transcript quantifications.
  * `P.quant` - a tab separated file listing the quantified targets, as well as information about their length and other metadata. The `num_reads` column provides the estimate of the number of reads originating from each target.
  * `P.infreps.pq` - a [`Parquet`](https://parquet.apache.org/) table where each row is a transcript and each column is an inferential replicate, containing the estimated counts for each transcript under each computed inferential replicate.

### References

[^Gleeson]: Josie Gleeson, Adrien Leger, Yair D J Prawer, Tracy A Lane, Paul J Harrison, Wilfried Haerty, Michael B Clark, Accurate expression quantification from nanopore direct RNA sequencing with NanoCount, Nucleic Acids Research, Volume 50, Issue 4, 28 February 2022, Page e19, [https://doi.org/10.1093/nar/gkab1129](https://doi.org/10.1093/nar/gkab1129)

[^preprint]: Zahra Zare Jousheghani, Rob Patro. Oarfish: Enhanced probabilistic modeling leads to improved accuracy in long read transcriptome quantification, bioRxiv 2024.02.28.582591; doi: [https://doi.org/10.1101/2024.02.28.582591](https://doi.org/10.1101/2024.02.28.582591)
