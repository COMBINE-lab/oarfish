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

**Convergence.** The truncation EM stops when no transcript with more than 0.01 reads changes by more than `--convergence-thresh` (relative, in either direction) between iterations. Transcripts below 0.01 reads are left out of the test (as salmon does): their abundances decay toward zero for thousands of iterations while holding a negligible number of reads. The test counts decreases as well as increases, unlike the default EM's, because the truncation EM is accelerated (SQUAREM), and after an accelerated step most of the motion is in decreases. It typically converges after 1,000–1,300 iterations, so with this model `--max-em-iter` defaults to 5,000 instead of 1,000. A run that reaches the limit before converging gets a [`truncation_em_not_converged` warning](/oarfish/guides/sequencing-technologies/#run-warnings), and `meta_info.json` records the iteration count under `em`.

**Read ends with the spearfish backend.** The model needs to know where each read starts and ends on the transcript, so with `--mapper spearfish` it turns on the extension of alignment spans to the read ends for the presets that leave spans at the outermost seeds (`--seq-tech ont-drna` and `pac-bio-hifi`; the cDNA preset extends them anyway). Longer spans give reads a larger aligned fraction, so the aligned-fraction filter (`-m`) drops a different, usually smaller, set of reads than in a run of the same preset without the model, and a comparison of the two runs is not purely a comparison of the models. The log says at the start of the run whether spans were extended and why, and `meta_info.json` records it as `spearfish_span_extension`. For a RAD file, the mapper decides: map with `spearfish map --extend-ends true` for this model (`oarfish` warns when the file's spans were not extended).

The model cannot be combined with the coverage model (`--coverage-model logistic`), and is not available in single-cell mode. With `--num-bootstraps`, each replicate resamples the reads and re-fits the placement parameters along with the abundances. It works in transcriptome-alignment, read, and genome (projection) modes; it is per read, and its EM is slower than the default equivalence-class EM (tens of seconds to minutes, depending on the number of alignments).

## Tail-neutral scores

Alignment scores are also made *tail-neutral* (for all models): when a read's poly(A) tail aligns to an A-run at the end of a transcript, the score credit those bases earn is removed, so that transcripts annotated with a poly(A) stretch do not attract reads from siblings without one. The aligner's match score is taken from the BAM header (minimap2) or the mapping preset.
