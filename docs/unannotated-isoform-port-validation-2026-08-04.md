# Validating the ported annotation-omission model (2026-08-04)

## What this validates

`--model-unannotated-isoforms` was ported onto the retired-coverage base
(`2d64016`) out of `wip/coverage-and-annotation-2026-07-25`, separating it from
three other threads that share the same files. This document checks that the
port reproduces the behaviour the original evaluation
(`docs/annotation-omission-evaluation-2026-07-25.md`, on that branch) measured.

**Result: both the detection and the accuracy claims reproduce.** An initial
report that accuracy regressed was a scoring error, documented in §4.

## Setup

Genome-alignment mode from the existing spliced-genome BAM
(`eval/parity/genome.bam`, 1,356,492 reads), GRCh38.p14 with
`GCF.pc_lncrna.matched.gtf`, `--filter-group no-filters`, 8 threads.

The read names are NanoSim-style (`NM_002738_..._aligned_...`), so the
generating transcript is known per read.

Holdout follows the original design: transcripts are deleted from the set the
reads were actually simulated *from*. `origin50` deletes a seeded random 50% of
the 30,245 expressed-and-annotated transcripts — 15,122 transcripts across
9,720 genes.

## 1. Projection accounting reproduces to the read

Under the correct annotation, the port reports:

| quantity | this port | original evaluation (§5) |
|---|---|---|
| reads projected | 1,356,492 | 1,356,492 |
| junction-informative | 900,174 | 900,174 |
| disagree with **every** candidate | 7,876 (0.87%) | 7,876 (0.87%) |

Exact agreement, so the junction-evidence path survived the merge intact.

## 2. Detection reproduces — high precision, very low recall

Flagged loci under the holdout, scored against the manifest (a locus is a true
positive if any member transcript's gene lost an expressed isoform):

| min flagged reads | loci | TP | FP | precision |
|---|---|---|---|---|
| 1 | 83 | 83 | 0 | **1.000** |
| 3 | 83 | 83 | 0 | **1.000** |
| 10 | 63 | 63 | 0 | **1.000** |
| 50 | 21 | 21 | 0 | **1.000** |

Recall is 89 of 9,720 genes = **0.0092**.

The original evaluation reported precision 0.913 → 1.000 over the same sweep and
recall 0.0123–0.0195. Precision here is 1.000 even at the loosest threshold
because the model's own `--novel-min-locus-reads` default of 5 already filters
the low end; loci below that are never created.

**False-positive floor.** Under the *correct* annotation the model still flags 7
loci — these are false positives by construction, and correspond to the 0.87%
all-mismatch rate attributable to alignment/projection error. Against 83 under
the holdout that is an **11.9x enrichment**, and 77 of the 83 are loci not
flagged under the correct annotation.

The "high precision, very low recall" conclusion is confirmed: the model
identifies *which* loci are incomplete with near-perfect precision, and finds
about 1% of them.

## 3. Localisation

Of 14,908 transcripts not belonging to any flagged locus, **91 (0.61%)** change
by more than 1e-6 between model-off and model-on. The model is well localised
and does not perturb the rest of the quantification.

## 4. The accuracy claim reproduces (after a scoring correction)

Mean |log2(est/truth)| over transcripts belonging to acted-on loci, with
estimates renormalised to the truth library size:

| scope | measure | model OFF → ON | change |
|---|---|---|---|
| acted-on loci (n=215) | mean \|log2 err\| | 1.5854 → 1.4158 | **−10.7%** |
| acted-on loci | median signed error | +0.5218 → +0.4190 | over-attribution reduced |
| acted-on loci | % over-attributed | 61.9% → 58.6% | |
| global (56,098 survivors) | CCC on log1p | 0.83578 → 0.83616 | **+0.00038** |
| global | MARD | 0.23330 → 0.23306 | −0.00024 |
| global | RMSE | 674.638 → 668.759 | −5.879 |
| survivors at affected genes (n=24,935) | CCC | 0.76668 → 0.76710 | +0.00042 |
| *original evaluation, origin50* | *acted-on* | *1.2255 → 1.0512* | *−14.2%* |
| *original evaluation, origin50* | *global CCC* | | *+0.00118* |

The direction and mechanism match the original throughout: the baseline
over-attributes (median +0.52 here, +0.43 there) and the model reduces it. The
magnitude is smaller, consistent with a differently-constructed holdout.

### The scoring error that initially inverted this

A first pass reported **+4.5% (worse)** and concluded the accuracy claim failed
to reproduce. That was a defect in the scoring, not in the feature.

The holdout drops assigned reads from 1,356,492 to 881,436. Comparing raw
estimated counts against a truth table for the *full* library therefore imposes
a uniform deflation on every surviving transcript — the acted-on population
measured a median signed error of **−1.70**, i.e. estimates roughly 3x below
truth, where the original evaluation's comparable figure was **+0.43**. On a
population that starved, diverting further mass to a novel state can only
increase the error, so the measured sign was an artifact of the missing
renormalisation.

The mismatch in baseline sign (−1.70 versus +0.43) was the signal that the two
analyses were not measuring the same thing; it should have been checked before
concluding non-reproduction. Renormalising both arms by the same factor,
computed over all surviving expressed transcripts so the normaliser is not
itself perturbed by the intervention, restores agreement.

A second hypothesis was also tested and rejected along the way: that the
regression came from genes losing *every* expressed isoform, which the original
calls harmful by construction. 211 of 215 acted-on transcripts already lie in
survivor-retaining genes, so that mixture was not the cause either.

## Bearing on the feature

The original document's own headline is *"not 'the novel state improves
accuracy' — that is a 4th-decimal global effect on top of an already-robust
baseline"*, but rather that the model **contributes reporting**: near-perfect
precision on which loci are incomplete. That claim is validated here at
precision 1.000 with an 11.9x enrichment over the false-positive floor. The
secondary accuracy claim also reproduces, at -10.7% error on acted-on loci
against the original's -14.2%.

## Operational note

The parallel M-step does not implement the novel latent state, so enabling
`--model-unannotated-isoforms` forces serial EM (`novel_loci == 0 &&
args.threads > 4`). On large inputs this costs parallelism silently. Inherited
behaviour, worth addressing before the flag is considered for any default.

## Reproducing

```
scripts/omission-validation-2026-08-04/make_holdout.py    # build origin50
scripts/omission-validation-2026-08-04/run_omission.sh    # three runs
scripts/omission-validation-2026-08-04/score_omission.py  # precision/recall
```
