---
title: Output files
description: The files oarfish writes, and the format of read-level assignment probabilities.
---

The `--output` option passed to `oarfish` corresponds to a path prefix (this prefix can contain the path separator character and if it refers to a directory that does not yet exist, that directory will be created). Based on this path prefix, say `P`, `oarfish` will create the following files:

- `P.meta_info.json` - a JSON format file containing information about relevant parameters with which `oarfish` was run, and other relevant information from the processed sample apart from the actual transcript quantifications. As in `salmon`, it records the version that produced the output (`oarfish_version`) and when the run started and ended (`start_time`, `end_time`).
- `P.quant` - a tab separated file listing the quantified targets, as well as information about their length and other metadata. The `num_reads` column provides the estimate of the number of reads originating from each target.
- `P.infreps.pq` - a [`Parquet`](https://parquet.apache.org/) table where each row is a transcript and each column is an inferential replicate, containing the estimated counts for each transcript under each computed [inferential replicate](/oarfish/guides/inferential-replicates/).
- `P.ambig_info.tsv` - a tab separated file listing, for each transcript (in the same order in which they appear in `P.quant`) the number of uniquely mapped, ambiguously mapped, and total reads. The quantification estimate for each transcript, in general, should reside between the number of uniquely aligned reads and the total number of reads (i.e. these provide, respectively lower and upper bounds for the number of reads assigned to each transcript). Note that the total in this file is the total number of reads that align to this transcript with a sufficiently high alignment score --- it is _not_, in general, an estimate of the number of reads originating from this transcript as many of those reads can be multimapping and, in fact, potentially better described by other transcripts.
- `P.prob[.lz4]` - a file encoding the assignment probability of each read to each transcript to which it had a valid alignment (optionally compressed using [`lz4`](https://github.com/lz4/lz4)). This file is optional and is generated only if `--write-assignment-probs` is passed to `oarfish`.

## Read-level assignment probabilities

`oarfish` has the ability to output read-level assignment probabilities. That is, for each input read, what is the probability, conditioned on the final estimate of transcript abundances, that the read was sequenced from each transcript to which it aligned. By default, this information is not recorded (as it's not required, or commonly used, for most standard analyses). To enable this output, you should pass the `--write-assignment-probs` option to `oarfish`. Optionally, you may also pass `--write-assignment-probs=compressed` to write the output to a compressed ([lz4](https://github.com/lz4/lz4)) stream --- the default output is to an uncompressed text file.

### Format

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

The reported probabilities are renormalized over the alignments that pass `--display-thresh`, so the `k` values on a line sum to 1. They are printed with enough decimal places that an assignment passing `--display-thresh` is never rounded to `0`: the precision tracks the threshold (the default `--display-thresh 1e-6` prints 6 decimal places, and thresholds of `1e-3` or larger keep 3), so lowering the threshold both reports more alignments and prints them at a resolution that can represent them.

### Example

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

### Compressed output

The compressed output (i.e. what is generated if one passes `--write-assignment-probs=compressed`) is exactly the same format, except instead of residing in a plain text file, it is written to an lz4 compressed text file. You can either decompress this file first with an lz4 decompressor, or decompress it on-the-fly as you are parsing the file using the lz4 library in your favorite language.
