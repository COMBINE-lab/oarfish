---
title: Usage examples
description: Quantify from a transcriptome BAM, from raw reads, or from the genome.
---

All of these commands produce the output files described in [Output files](/oarfish/reference/output/). For how the modes differ, see [Quantification modes](/oarfish/guides/quantification-modes/).

## Alignment mode (transcriptome BAM)

Assume you have ONT cDNA reads in `sample1_reads.fq.gz` and a *transcriptome* reference `transcripts.fa`, and you have [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/) installed:

```bash
$ minimap2 -t 16 --eqx -N 100 -ax map-ont transcripts.fa sample1_reads.fq.gz | samtools view -@4 -b -o alignments.bam
$ oarfish -j 16 -a alignments.bam -o sample1 --filter-group no-filters --model-coverage
```

## Read mode (oarfish maps to the transcriptome)

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

To map with [spearfish](/oarfish/guides/mapping-backends/#spearfish) instead of rammap (the index is saved and reused the same way, as a spearfish `.sfi`):

```bash
$ oarfish -j 16 --reads sample1_reads.fq.gz --annotated transcripts.fa --seq-tech ont-cdna --mapper spearfish -o sample1
# or: write the index once, then reuse it
$ oarfish --only-index --mapper spearfish --annotated transcripts.fa --index-out transcripts.sfi --seq-tech ont-cdna
$ oarfish -j 16 --reads sample2_reads.fq.gz --index transcripts.sfi --seq-tech ont-cdna --mapper spearfish -o sample2
```

## Genome read-projection (reads → genome → transcripts)

Spliced-align the reads to the genome and project onto the annotated transcripts. You need a genome FASTA (or a prebuilt genome index) and a transcript annotation (GTF/GFF):

```bash
$ oarfish -j 16 --reads sample1_reads.fq.gz --genome genome.fa --annotation annotation.gtf \
    --seq-tech pac-bio-hifi -o sample1 --filter-group no-filters
```

Soft-clip rescue is on by default and, when `--genome` is a FASTA, the rescue reference sequence is taken from it automatically. Add `--no-rescue` to disable rescue, or `--junctions models.bed12` to supply your own splice-junction hints instead of deriving them from the annotation. See [Genome-based input](/oarfish/guides/input/#genome-based-input-projection-onto-transcripts).

## Genome-alignment projection (an existing genome BAM)

If you already have a spliced genome BAM (e.g. produced with `minimap2 -ax splice`), project it onto the transcripts directly — no reads or genome FASTA required for the projection itself, though a genome FASTA is needed if you want soft-clip rescue:

```bash
# BAM must be collated by read name
$ minimap2 -t 16 --eqx -ax splice genome.fa sample1_reads.fq.gz | samtools collate -@4 -O -u - | samtools view -b -o genome.bam
$ oarfish -j 16 --genome-alignments genome.bam --annotation annotation.gtf --genome-fasta genome.fa \
    -o sample1 --filter-group no-filters
```
