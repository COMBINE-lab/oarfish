# oarfish: transcript quantification from long-read RNA-seq data

The usage can be provided by passing `-h` at the command line.
```
Usage: oarfish [OPTIONS] --alignments <ALIGNMENTS> --output <OUTPUT>

Options:
  -a, --alignments <ALIGNMENTS>
          Name of the person to greet
  -o, --output <OUTPUT>
          Location where output quantification file should be written
  -t, --three-prime-clip <THREE_PRIME_CLIP>
          [default: 4294967295]
  -f, --five-prime-clip <FIVE_PRIME_CLIP>
          [default: 4294967295]
  -s, --score-threshold <SCORE_THRESHOLD>
          [default: 0.95]
  -m, --min-aligned-fraction <MIN_ALIGNED_FRACTION>
          [default: 0.5]
  -l, --min-aligned-len <MIN_ALIGNED_LEN>
          [default: 50]
  -n, --allow-negative-strand

      --model-coverage

  -h, --help
          Print help
  -V, --version
          Print version
```

The input should be a `bam` format file, with reads aligned using `minimap2` against the _transcriptome_. That is, `oarfish` does not currently handle spliced alignment to the genome.  Further, the output alignments should be name sorted (the default order produced by `minimap2` should be fine).

The output is a tab separated file listing the quantified targets, as well as information about their length and other metadata. The `num_reads` column provides the estimate of the number of reads originating from each target.
