---
title: Quantification modes
description: The five ways of giving oarfish input for bulk quantification, and which to use.
---

`oarfish` supports five ways of providing input for **bulk** quantification. All of them feed the same EM-based quantification core; they differ only in how alignments to the transcriptome are obtained. The mode is selected implicitly by which input flag you pass. (Single-cell quantification, `--single-cell`, is separate and accepts **only** a transcriptome BAM — see [Single-cell mode](/oarfish/guides/single-cell/).)

| Mode | Selected by | What it does |
|------|-------------|--------------|
| **Alignment** (transcriptome BAM) | `--alignments` | Quantify pre-computed alignments of reads to the **transcriptome** (a name-sorted BAM, e.g. from `minimap2`/`pbmm2`). |
| **Read** (transcriptome mapping) | `--reads` + `--annotated`/`--index` | Map raw reads to the **transcriptome** internally (with the built-in [rammap](/oarfish/guides/mapping-backends/) mapper, or [spearfish](/oarfish/guides/mapping-backends/#spearfish) with `--mapper spearfish`), then quantify. |
| **RAD** (spearfish mappings) | `--alignments` with a RAD file | Quantify the mappings in a RAD file written by [spearfish](/oarfish/guides/mapping-backends/#spearfish) (`spearfish map --rad`). The file is recognized by its contents, and carries the transcripts, the score scale, the preset and whether spans reach the read ends, so no other options are needed. |
| **Genome read-projection** | `--reads` + `--genome` + `--annotation` | Spliced-align raw reads to the **genome**, project each read's genomic alignments onto the annotated transcripts (via [bramble](https://github.com/COMBINE-lab/bramble)), then quantify. |
| **Genome-alignment projection** | `--genome-alignments` + `--annotation` | Project an existing **spliced genome BAM** (e.g. from `minimap2 -ax splice`) onto the annotated transcripts, then quantify. |

:::tip[Which should I use?]
The transcriptome read/alignment modes are the most mature and are the right default when you have (or can build) a transcriptome reference. The genome-based modes are useful when you want to align against the genome — for example to keep one genomic alignment pass for multiple downstream analyses, to start from genome BAMs you already have, or to capture reads that fall outside annotated transcript boundaries (recovered via soft-clip [rescue](/oarfish/guides/input/#soft-clip-rescue)). On our benchmarks, genome read-projection quantification closely matches direct transcriptome quantification when rescue is enabled (see [expectations](/oarfish/guides/input/#expectations-for-the-genome-based-modes)).
:::

For details of each input, see [Input](/oarfish/guides/input/); for worked commands, see [Usage examples](/oarfish/getting-started/usage-examples/).
