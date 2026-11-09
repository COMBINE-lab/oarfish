---
title: Single-cell mode
description: Single-cell quantification from a transcriptome BAM collated by cell barcode.
---

Starting with version 0.6.1 `oarfish` incorporates the first single-cell quantification capabilities. Given a `bam` file, **collated by cell barcode and with already (UMI) deduplicated reads**, this mode, enabled with the `--single-cell` flag, will allow `oarfish` to produce a single-cell quantification matrix.

:::caution[Single-cell mode requires a transcriptome-aligned BAM, supplied with `--alignments`]
Unlike bulk quantification — which can take a transcriptome BAM (`--alignments`), raw reads (`--reads`), a genome FASTA + reads (`--genome`), or a genome BAM (`--genome-alignments`) — single-cell mode supports **only** the transcriptome-BAM path. It does not support raw-read mapping or either of the genome (projection) modes. The CLI enforces this: `--single-cell` requires `--alignments` and conflicts with `--reads`, `--genome`, and `--genome-alignments`. The input BAM must therefore be reads aligned to the *transcriptome* (e.g. with minimap2/pbmm2), not to the genome.
:::

**Formatting requirements of BAM input in single-cell mode**: All alignment records for the same cell barcode should be adjacent in the `bam` file, and a count will be obtained for each read record, so UMI de-duplication should have been performed if those are the counts you want. In the future, counting UMIs directly may be supported, and some of these other restrictions may be lifted.
