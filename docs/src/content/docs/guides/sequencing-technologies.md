---
title: Recommended practices by sequencing technology
description: How to prepare ONT and PacBio reads for oarfish, which options to set for each, and what the run warnings mean.
---

What `oarfish` needs from the reads differs a little between technologies. This page collects our recommendations, and explains the [warnings](#run-warnings) `oarfish` gives when a run looks as if one of them was missed.

In short:

| Reads | Before quantification | Options |
|-------|----------------------|---------|
| ONT cDNA (PCR or direct) | trim adapters and primers; **keep the poly(A) tails** | `--seq-tech ont-cdna` |
| ONT direct RNA | trim the adapter (the basecaller does) | `--seq-tech ont-drna` |
| PacBio HiFi, oriented (Iso-Seq / Kinnex FLNC) | nothing more | `--seq-tech pac-bio-hifi`, and `--map-strand fw` with `--mapper spearfish` |
| PacBio HiFi, not oriented | nothing more | `--seq-tech pac-bio-hifi`; **no** `--map-strand` or `-d` |

## Oxford Nanopore cDNA

### Trim adapters and primers

A cDNA read starts and ends with sequence that is not part of the transcript: the sequencing adapters, and the primers of the library kit (the strand-switching primer at one end, the anchored poly(T) primer at the other). If it is left on the reads, it cannot align, and a read with enough of it fails the aligned-fraction filter (`-m`/`--min-aligned-fraction`, 0.5 by default): the read is dropped although it maps well. Short reads are hit hardest, since the same primers are a larger part of them.

- Reads basecalled with a current `dorado` are trimmed of adapters and primers by default (unless it ran with `--no-trim`). Nothing more is needed.
- Reads from older basecallers, reads basecalled with trimming off, and many reads downloaded from public archives still carry them. Trim these before quantification with a tool that removes adapters and primers, such as `dorado trim`, [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI), [cutadapt](https://cutadapt.readthedocs.io/) with the kit's primer sequences, or [Pychopper](https://github.com/epi2me-labs/pychopper).

If you cannot tell which case you are in, run `oarfish` and look at the discard table: when the aligned-fraction filter drops more than a tenth of the mapped reads, `oarfish` [warns](#run-warnings), and untrimmed reads are the usual cause.

Trimming is better than turning the filter off. `--min-aligned-fraction 0` (or `--filter-group no-filters`) keeps the reads, and is a reasonable choice when the reads cannot be trimmed again, but the filter also removes reads that are worth removing, such as chimeric reads.

### Do not trim the poly(A) tails

Remove adapters and primers, but leave the poly(A) (or, in reads of the other strand, poly(T)) tails on the reads. Trimming them does not help quantification, and with the spearfish backend it makes the estimates for isoforms worse: the tail marks where the transcript ends, and that helps to tell apart isoforms that differ at their 3′ ends.

Some trimmers remove the tails by default or as an option; turn that off.

### Keeping only full-length reads

Tools such as Pychopper can also classify reads, and keep only those with both primers ("full-length" reads). For PCR cDNA this discards few reads and is harmless. For direct cDNA, and for older libraries in general, it can discard most of the reads, and those reads carry usable signal: prefer to trim and keep all reads, rather than to keep only the full-length ones.

### Orientation

cDNA reads come from either strand, so `oarfish` maps both orientations and keeps the alignments of both, which is the default. Leave `--map-strand` and `-d`/`--strand-filter` alone, unless *every* read was oriented to the transcript strand beforehand (by a tool that finds the primers and turns the read, and that was told to discard the reads it could not orient). Then the advice [for oriented PacBio reads](#oriented-reads-iso-seq-and-kinnex-flnc) applies.

## Oxford Nanopore direct RNA

Direct RNA reads are the RNA molecule itself, so all of them are in the orientation of the transcript. `--seq-tech ont-drna` selects the presets for them; with `--mapper spearfish`, only the orientation the read is given in is mapped. The basecaller removes the adapter; no other trimming is needed, and the poly(A) tail stays on the read.

## PacBio

HiFi reads (`--seq-tech pac-bio-hifi`) that went through the Iso-Seq workflow have had their primers removed, so there is nothing to trim. What matters is whether they are oriented.

### Oriented reads (Iso-Seq and Kinnex FLNC)

Full-length non-concatemer (FLNC) reads, as `isoseq refine` writes them, are all oriented to the transcript strand. For such reads it is useful to say so:

- With `--mapper spearfish`, pass `--map-strand fw`. Only the orientation the reads are in is mapped, which takes less mapping time and changes little else.

  ```bash
  $ oarfish --reads flnc.fastq.gz --index transcripts.sfi --mapper spearfish --seq-tech pac-bio-hifi \
      --map-strand fw -j 16 -o sample1
  ```

- With the default backend (rammap) or with BAM input, there is no such mapping option. `-d fw` (`--strand-filter fw`) removes the alignments in the other orientation after mapping; that saves no time, and removes the few spurious antisense alignments.

### Reads that are not oriented

Do **not** pass `--map-strand fw` (or `-d fw`) for reads that were not oriented: HiFi/CCS reads that did not go through `lima` and `isoseq refine`, or any set in which reads of both strands remain. About half of such reads are the reverse complement of their transcript. With `--map-strand fw` those reads either do not map, or map only where some part of them happens to match a transcript in the mapped orientation (a gene on the opposite strand, a repeat), which is a wrong assignment. With `-d fw` their alignments are removed. Either way about half of the data is lost; `oarfish` [warns](#run-warnings) when it sees this.

The defaults (both orientations mapped, both kept) are right for these reads.

## Run warnings

When a run finishes, `oarfish` checks the counts it kept of what happened to the reads, and reports conditions that deserve a look. Each warning is written to the log, and is recorded in the `warnings` array of [`meta_info.json`](/oarfish/reference/output/#run-warnings), so that it stays with the results and a pipeline can check for it. A warning does not stop the run or change the estimates.

| Code | Issued when | What to check |
|------|-------------|---------------|
| `aligned_fraction_filter` | more than 10% of the mapped reads were dropped because too little of the read is in its best alignment | Are adapters and primers still on the reads? [Trim them](#trim-adapters-and-primers), or keep the reads with `--min-aligned-fraction 0`. |
| `strand_filter` | `-d fw` or `-d rc` was given, and more than 10% of the mapped reads align only to the strand it removes | Are the reads really all in one orientation? If not, use the default (`-d both`). |
| `alignment_filters` | more than 10% of the mapped reads have no alignment that passes the other alignment filters (aligned length, distance from the transcript ends, score) | The filter settings, recorded under `filter_options`. Filters meant for one kind of data (say, a 3′-end distance limit meant for direct RNA) can remove much of another. |
| `map_strand_unmapped` | `--map-strand fw` or `rc` was given, and more than 30% of the reads did not map | Are the reads really oriented? If not, drop `--map-strand`. |
| `low_mapping_rate` | more than half of the reads did not map | Does the reference match the sample, and `--seq-tech` the reads? |

The thresholds are deliberately loose: a run without warnings is not thereby a good one, and the discard table, in the log and in `meta_info.json`, has the full counts.
