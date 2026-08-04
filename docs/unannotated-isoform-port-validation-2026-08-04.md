# Validating the ported annotation-omission model (2026-08-04)

## What this validates

`--model-unannotated-isoforms` was ported onto the retired-coverage base
(`2d64016`) out of `wip/coverage-and-annotation-2026-07-25`, separating it from
three other threads that share the same files. This document checks that the
port reproduces the behaviour the original evaluation
(`docs/annotation-omission-evaluation-2026-07-25.md`, on that branch) measured.

**Result: the detection claim reproduces exactly. The accuracy claim does not.**

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

## 4. The accuracy claim does NOT reproduce

Mean |log2(est/truth)| over transcripts belonging to acted-on loci:

| population | n | model OFF → ON | change |
|---|---|---|---|
| all acted-on | 215 | 1.8413 → 1.9207 | **+4.3% (worse)** |
| genes retaining a survivor | 211 | 1.8263 → 1.9082 | **+4.5% (worse)** |
| genes losing every isoform | 0 | — | — |
| *original evaluation, origin50* | *179* | *1.2255 → 1.0512* | *−14.2%* |

**A hypothesis was tested and rejected.** The original document notes that
designs deleting *every* isoform of a gene (`locus*`) are "harmful by
construction" (+42.0% on locus25) because no survivor exists to absorb the
correction. 52.3% of the genes affected by this holdout lose all their expressed
isoforms, so the mixture seemed a likely explanation. It is not: splitting the
two populations shows 211 of 215 acted-on transcripts already lie in
survivor-retaining genes, and the regression is present there at +4.5%. Genes
losing every isoform contribute nothing to the acted-on set, because they retain
no surviving member that could appear in a flagged locus.

Remaining differences from the original setup, none of which is established as
the cause:

- the holdout differs in size (15,122 deleted / 9,720 genes here versus 11,904 /
  7,969 there), and the baseline error differs correspondingly (1.83 versus
  1.23), so the two are not scoring the same transcript population;
- the original holdout-construction scripts were not preserved, so "origin50"
  cannot be reproduced exactly — only re-derived from its description;
- this run uses the shipped defaults (`--novel-odds-per-miss 2.0`,
  `--novel-min-misses 1`, `--novel-min-locus-reads 5`); the original sweep's
  settings are not recorded per table.

**Status: unverified, not refuted.** The accuracy effect may be specific to a
holdout construction this reproduction did not match. It should not be cited
until it is reproduced against a preserved holdout.

## Bearing on the feature

The original document's own headline is *"not 'the novel state improves
accuracy' — that is a 4th-decimal global effect on top of an already-robust
baseline"*, but rather that the model **contributes reporting**: near-perfect
precision on which loci are incomplete. That claim is the one validated here.
The secondary accuracy claim is the one that failed to reproduce, which leaves
the feature's stated primary contribution intact.

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
