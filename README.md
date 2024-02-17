# oarfish: Accurate, fast, and versatile transcript quantification from long-read RNA-seq data

`oarfish` is a program, written in `rust`, for quantifying transcript-level expression from long-read (i.e. Oxford nanopore cDNA and direct RNA and PacBio) sequencing technologies. `oarfish` currently requires a sample of sequencing reads aligned to the *transcriptome* (not the genome). It handles multi-mapping reads through the use of probabilistic allocation via an expectation-maximization (EM) algorithm.  Further, it employs many filters to help discard alignments that may reduce quantification accuracy.  Currently, the set of filters applied in `oarfish` are directly derived from the [`NanoCount`](https://github.com/a-slide/NanoCount)[^Gleeson] tool; both the filters that exist, and the way their values are set (with the exception of the `--three-prime-clip` filter, which is not set by default in `oarfish` but is in `NanoCount`).  

Additionally, as a core component of `oarfish` we are developing novel coverage models that further improve quantification estimates by taking into account the coverage profile of the long mapped reads over the transcripts.  A detailed writeup of these coverage models is currently being drafted. In the meantime, if you make use of `oarfish`, please cite this repository.

The usage can be provided by passing `-h` at the command line.
```
accurate transcript quantification from long-read RNA-seq data

Usage: oarfish [OPTIONS] --alignments <ALIGNMENTS> --output <OUTPUT>

Options:
      --quiet                    be quiet (i.e. don't output log messages that aren't at least warnings)
      --verbose                  be verbose (i.e. output all non-developer logging messages)
  -a, --alignments <ALIGNMENTS>  path to the file containing the input alignments
  -o, --output <OUTPUT>          location where output quantification file should be written
  -t, --threads <THREADS>        maximum number of cores that the oarfish can use to obtain binomial probability [default: 1]
  -h, --help                     Print help
  -V, --version                  Print version

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
  -n, --allow-negative-strand
          allow both forward-strand and reverse-complement alignments

coverage model:
      --model-coverage  apply the coverage model
  -b, --bins <BINS>     number of bins to use in coverage model [default: 10]

EM:
      --max-em-iter <MAX_EM_ITER>
          maximum number of iterations for which to run the EM algorithm [default: 1000]
      --convergence-thresh <CONVERGENCE_THRESH>
          maximum number of iterations for which to run the EM algorithm [default: 0.001]
  -q, --short-quant <SHORT_QUANT>
          location of short read quantification (if provided)
```

The input should be a `bam` format file, with reads aligned using `minimap2` against the _transcriptome_. That is, `oarfish` does not currently handle spliced alignment to the genome.  Further, the output alignments should be name sorted (the default order produced by `minimap2` should be fine).

The output is a tab separated file listing the quantified targets, as well as information about their length and other metadata. The `num_reads` column provides the estimate of the number of reads originating from each target.


### References

[^Gleeson]: Josie Gleeson, Adrien Leger, Yair D J Prawer, Tracy A Lane, Paul J Harrison, Wilfried Haerty, Michael B Clark, Accurate expression quantification from nanopore direct RNA sequencing with NanoCount, Nucleic Acids Research, Volume 50, Issue 4, 28 February 2022, Page e19, [https://doi.org/10.1093/nar/gkab1129](https://doi.org/10.1093/nar/gkab1129)
