---
title: Coverage and truncation models
description: Optional models that use where reads land on each transcript.
---

By default, `oarfish` weighs each of a read's alignments by its alignment score. Two optional models also use where reads land on each transcript; they are chosen with `--coverage-model` and cannot be combined.

## Coverage model

`--coverage-model logistic` (or `--model-coverage`) makes use of coverage profiles derived from the aligned reads to improve quantification accuracy. Its options are the logistic growth rate (`-k`/`--growth-rate`, default 2) and the coverage bin width (`-b`/`--bin-width`, default 100).

## Truncation model

Long-read protocols produce about one read per captured molecule, so `oarfish` does not divide abundances by an effective length. What does depend on the transcript is where a read can land on it, and that differs between protocols: direct RNA and poly(A)-primed cDNA reads end at the 3' end and are truncated toward the 5' end; some libraries are dominated by full-length reads; others are truncated at both ends. `--coverage-model truncation` (or `--model-truncation`) multiplies each alignment's score probability by the probability of its placement on the transcript under a mixture learned from the data, so it needs no setting for the protocol:

- **3'-anchored** reads end at the transcript's 3' end and extend toward the 5' end by a length drawn from a learned (nonparametric) length law;
- **5'-anchored** reads start at the 5' end, likewise;
- **full-length** reads span the transcript;
- **uniformly placed** reads start anywhere, with their own learned length law.

Where a read is "at" an end, the gap to that end follows a learned end-gap law (nonparametric on log-spaced bins, with a part that scales with transcript length), which absorbs aligner clipping and annotation imprecision; a read longer than its transcript is treated as censored. The mixture weights and laws are estimated in the EM together with the abundances, and the fitted mixture is written to the log. With more than a million reads they are estimated on a sample of about a million (every k-th read) and then held fixed while the abundances are estimated from all reads, which gives the same estimates as fitting them on all reads in about half the time (`OARFISH_TRUNC_SAMPLE=0` fits them on all reads). Because each transcript's placement distribution sums to one over the placements it can produce, a read that fits a short isoform end to end but would be truncated on a longer one favors the short isoform, by as much as the fitted mixture says truncation is unlikely.

The model cannot be combined with the coverage model (`--coverage-model logistic`), and is not available in single-cell mode. With `--num-bootstraps`, each replicate resamples the reads and re-fits the placement parameters along with the abundances. It works in transcriptome-alignment, read, and genome (projection) modes; it is per read, and its EM is slower than the default equivalence-class EM (tens of seconds to minutes, depending on the number of alignments).

## Tail-neutral scores

Alignment scores are also made *tail-neutral* (for all models): when a read's poly(A) tail aligns to an A-run at the end of a transcript, the score credit those bases earn is removed, so that transcripts annotated with a poly(A) stretch do not attract reads from siblings without one. The aligner's match score is taken from the BAM header (minimap2) or the mapping preset.
