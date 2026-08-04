# Evaluating the annotation-omission model (2026-07-25)

All runs: ONT cDNA, 2M simulated reads (`cdna_s1_reads_shuffled.2M.fastq`), genome read
mode (spliced genome alignment → bramble projection), `--filter-group no-filters`,
16 threads, GRCh38.p14 + `GCF.pc_lncrna.matched.gtf`.

Holdout GTFs follow TranSigner: transcripts are deleted from the set reads were
actually simulated *from*, not from the catalogue at large (deleting unexpressed
transcripts orphans no reads and measures nothing).

> **Caveat on cross-table comparison.** Each table computes metrics over the set of
> transcripts common to the arms *in that table*. The three-arm table intersects with
> the transcriptome-mode quant, so its transcript set — and therefore its absolute CCC
> — differs from the deep-sweep table. Compare within a table, not across tables.

---

## 1. The original evidence was diluted, not absent

Global CCC over ~20k transcripts, when the model touches ~50 loci, makes a 4th-decimal
delta arithmetically inevitable. Restricting to transcripts that are members of a
flagged locus:

| design | n txps | mean \|log2 err\| full → omission | change |
|---|---|---|---|
| dominant25 | 143 | 1.1203 → 1.0064 | **−10.2%** |
| origin25 | 131 | 1.0982 → 0.9585 | **−12.7%** |
| origin50 | 179 | 1.2255 → 1.0512 | **−14.2%** |
| minor25 | 30 | 0.5993 → 0.6931 | +15.7% (worse) |
| minor50 | 29 | 0.6158 → 0.7561 | +22.8% (worse) |
| locus25 | 17 | 0.3965 → 0.5632 | +42.0% (worse) |

`locus*` deletes *every* isoform of a gene, so no survivor exists to correct — flagging
there is harmful by construction. `minor*` (keep only the dominant isoform) also
regresses, on small n. **These are real limitations, not noise to be explained away.**

An intermediate restriction — all surviving transcripts at *affected genes* — is still
~10× too coarse (4,770 affected genes vs 409 acted-on loci) and shows almost nothing:

| design | group | n | mean \|log2\| full → omission |
|---|---|---|---|
| origin50 | surviving @ affected genes | 6,347 | 0.9391 → 0.9348 |
| origin50 | surviving @ unaffected genes | 5,558 | 0.5735 → 0.5703 |

## 2. The mechanism is the claimed one (over-attribution)

Signed error at acted-on loci; positive = over-attributed.

| design | median log2(est/truth) | % inflated >2× |
|---|---|---|
| dominant25 | +0.0942 → +0.1249 | 28.0% → 28.7% |
| origin25 | +0.2023 → +0.1910 | 31.3% → 26.0% |
| origin50 | **+0.4329 → +0.2909** | **39.1% → 33.0%** |

## 3. Dose-response holds to extreme incompleteness

New holdouts at 75% and 90%. Global CCC, and error reduction on acted-on loci.

| design | orphaned mass | CCC full | CCC omission | delta | flagged loci | acted-on error reduction |
|---|---|---|---|---|---|---|
| dominant50 | 21.5% | 0.85642 | 0.85722 | +0.00080 | 70 | −13.6% (n=196) |
| dominant75 | 31.2% | 0.84173 | 0.84309 | +0.00136 | 98 | −17.4% (n=270) |
| dominant90 | 40.0% | 0.83083 | 0.83243 | +0.00160 | 112 | −17.2% (n=311) |
| origin75 | 76.0% | 0.83310 | 0.83515 | +0.00205 | 78 | −16.2% (n=104) |
| origin90 | **87.6%** | 0.78673 | 0.78925 | **+0.00252** | 53 | **−19.0%** (n=37) |

Monotone in incompleteness — the effect roughly triples from 21.5% to 87.6% orphaned.

## 4. Detection is the stronger claim: high precision, very low recall

Treating "is this locus incomplete?" as classification against holdout truth. A locus is
a true positive if any member transcript's gene lost an expressed isoform.

**origin50** (7,969 genes lost ≥1 expressed isoform):

| min flagged reads | loci | precision | genes found | recall |
|---|---|---|---|---|
| 1 | 161 | 0.913 | 155 | 0.0195 |
| 2 | 104 | 0.971 | 110 | 0.0138 |
| 3 | 89 | **1.000** | 98 | 0.0123 |
| 5 | 75 | **1.000** | 82 | 0.0103 |
| 10 | 55 | **1.000** | 60 | 0.0075 |
| 50 | 27 | **1.000** | 31 | 0.0039 |

**origin25** (4,770 affected): precision 0.787 → 1.000 over the same sweep, recall
0.0180 → 0.0031. **dominant25** (1,413 affected): precision 0.522 → 0.944, recall
0.0354 → 0.0127.

Recall **saturates at 2–3.5%** even at the most permissive threshold: signal-limited,
not threshold-limited.

### Recall is not limited by the annotation

Is the deleted transcript's junction set a subset of some surviving sibling's (hence
undetectable in principle)?

| design | lost txps | single-exon | juncs ⊆ survivor | detectable ceiling |
|---|---|---|---|---|
| origin25 | 5,952 | 78 | 283 | **93.9%** |
| dominant25 | 1,413 | 2 | 104 | **92.5%** |
| origin50 | 11,904 | 206 | 349 | **95.3%** |

93–95% are distinguishable in principle; we recover 2–4%.

## 5. Where the reads actually go

Junction-evidence accounting from the projection stage:

| annotation | reads projected | junction-informative | disagree with **every** candidate |
|---|---|---|---|
| full | 1,356,492 | 900,174 (66.4%) | 7,876 (0.87%) |
| origin25 | 1,160,521 | 750,580 (64.7%) | 13,135 (1.75%) |
| origin50 | 871,411 | 525,168 (60.3%) | 16,350 (3.11%) |

Reads vanishing from projection entirely vs the full annotation: **195,971** (origin25),
**485,081** (origin50). The junction-mismatch channel the model consumes grows by only
**5,259** and **8,474** over its full-annotation floor — 37× and 57× smaller.

The 0.87% all-mismatch rate under the *correct* annotation is the false-positive floor
(alignment/projection error). Signal rises to 1.75%/3.11% — only a 2–3.5× enrichment,
which caps how aggressive the novel state can safely be.

## 6. Three-arm comparison

| design | arm | CCC | MARD | reads assigned |
|---|---|---|---|---|
| origin25 | transcriptome (no compatibility filter) | 0.86916 | 0.3217 | 1,214,436 |
| | genome: filter only | 0.87772 | 0.3090 | 1,142,652 |
| | genome: filter + model | 0.87818 | 0.3087 | 1,140,903 |
| origin50 | transcriptome | 0.83729 | 0.3640 | 957,195 |
| | genome: filter only | 0.86512 | 0.3289 | 833,605 |
| | genome: filter + model | 0.86630 | 0.3282 | 830,166 |
| origin90 | transcriptome | 0.76788 | 0.4273 | 351,302 |
| | genome: filter only | 0.82558 | 0.3661 | 282,305 |
| | genome: filter + model | 0.82764 | 0.3648 | 281,285 |

Genome-vs-transcriptome gap: +0.00856 / +0.02783 / +0.05770. **Confounded** — the arms
differ in aligner (direct transcriptome mapping vs spliced genome alignment +
projection), not only in filtering. Sections 7–8 remove that confound.

## 7. The similarity threshold is inert (negative result)

Exposed `similarity_threshold` on bramble's `ProjectionConfig` and as oarfish's hidden
`--projection-similarity-threshold`. Sweeping it on origin25:

| threshold | CCC | reads assigned |
|---|---|---|
| 0.60 (long-read default) | 0.878672 | 1,142,806 |
| 0.30 | 0.878285 | 1,142,838 |
| 0.05 | 0.878073 | 1,142,847 |

A 12× reduction admits **41 more reads out of 1.14M** (0.004%) and moves CCC by −0.0006.

**The similarity threshold accounts for essentially none of oarfish's robustness to
annotation incompleteness.** Reads orphaned by a deleted isoform do not produce
low-similarity candidates that get filtered — they produce *no viable exon-chain match
at all*, so structural matching rejects them before similarity is ever evaluated.

## 8. The junction tolerances are the real filter

Exposed `max_clip`, `max_junc_ins`, `max_junc_gap`, `max_error_exon`. Sweeping the
junction triple (`max_clip` left at its default so the intervention targets splice
disagreement specifically):

| design | tolerance | reads | CCC | MARD | vs default |
|---|---|---|---|---|---|
| origin25 | 40/40/35 (default) | 1,142,806 | 0.87867 | 0.3080 | — |
| | 120/120/120 | 1,155,507 | 0.87557 | 0.3118 | +12,701 reads, CCC −0.00311 |
| | 255/255/255 | 1,160,352 | 0.87341 | 0.3136 | +17,546 reads, CCC −0.00527 |
| origin50 | 40/40/35 | 833,792 | 0.86543 | 0.3288 | — |
| | 120/120/120 | 853,473 | 0.85889 | 0.3368 | +19,681 reads, CCC −0.00654 |
| | 255/255/255 | 864,613 | 0.85403 | 0.3419 | +30,821 reads, CCC −0.01140 |
| origin90 | 40/40/35 | 282,305 | 0.82558 | 0.3661 | — |
| | 120/120/120 | 303,290 | 0.80941 | 0.3875 | +20,984 reads, CCC −0.01617 |
| | 255/255/255 | 309,105 | 0.80402 | 0.3961 | +26,800 reads, CCC −0.02157 |

## 9. Final decomposition

| mechanism | origin25 | origin50 | origin90 |
|---|---|---|---|
| junction tolerance (40/40/35) | **+0.00527** | **+0.01140** | **+0.02157** |
| annotation-omission novel state | +0.00046 | +0.00118 | +0.00206 |
| similarity threshold (0.60) | ~0 (41 reads) | — | — |

The tolerance filter beats the novel-state model by a consistent **~10:1** at every
level of incompleteness, and both grow monotonically with it.

**Most exclusion is unbridgeable, not tunable.** At maximal permissiveness (255, the
`u8` ceiling) only 17.5k–30.8k extra reads are admitted, against ~196k (origin25) and
~485k (origin50) reads lost relative to the full annotation. So >90% of orphaned reads
have no viable exon chain in *any* surviving transcript — skipped exons are kb-scale
while tolerances cap at 255 nt. That dominant term is structural, not a filter setting.

The tolerance filter explains 62% / 41% / 37% of the genome-vs-transcriptome gap in §6;
the remainder is the aligner/mode difference plus that structural term. Note also that
relaxing tolerance can change which transcript wins for reads that *already* projected,
so this is not purely an "admit extra reads" experiment — but it is same-aligner,
same-mode, which was the confound it was designed to remove.

---

## What this means for the paper

The honest headline is **not** "the novel state improves accuracy" — that is a
4th-decimal global effect on top of an already-robust baseline. It is:

1. **Genome-projection mode is intrinsically robust to annotation incompleteness**,
   because structurally incompatible reads are excluded rather than misattributed.
   Holding CCC 0.826 with 87.6% of library mass orphaned is a strong result by itself.
2. **The annotation-omission model contributes reporting** — near-perfect precision on
   *which* loci are incomplete, a capability no competing tool offers — plus a 10–19%
   error reduction on the loci it acts on.

Recall (2–4%) is the honest weak point and should be stated as such.

## Code state

- bramble `ProjectionConfig` gained `similarity_threshold`, `max_clip`, `max_junc_ins`,
  `max_junc_gap`, `max_error_exon` (all `Option`, `None` = existing presets); passed
  through at the single `ReadEvaluator` construction site in `src/api.rs`.
- oarfish gained five matching hidden `--projection-*` flags.
- **Control:** the default path is unchanged to last-ULP. The byte-identity gate did
  *not* pass as written — 160 of 65,268 entries differ — but max absolute difference is
  **1.8e-12** on values of order 10⁶, with column sums exactly equal. That is
  accumulation-order nondeterminism in the 16-thread run; a read being admitted or
  excluded moves mass by O(1), not 10⁻¹².
- All uncommitted. The bramble **path dependency** in `oarfish/Cargo.toml` remains the
  outstanding release blocker.

## Deferred

Tag loci with non-trivial unprojectable reads → local transcript discovery/assembly →
re-quantify overlapping annotated and novel transcripts. Recorded separately; deserves
its own push.
