# oarfish 2 — first feature-phase results (2026-08-16)

Continues `oarfish2-measurement-report-2026-08-15.md`. Both features are on
`oarfish2-dev` (a87e68f), flag-gated, defaults unchanged.

## F2 — `--presence-model spike-slab` (Lead A)

Per-transcript Bernoulli presence posterior `q_t` inside the EM: the M-step
sees `q_t · counts[t]` (zero-inflated EM); `q` refreshed every 10 evaluations
after warmup from a leave-one-out log-evidence pass; fixed prior
`--presence-rho` (default 0.05). Design notes learned the hard way:

- The leave-one-out log-likelihood gain is **always ≥ 0**, so an
  empirical-Bayes `rho = mean(q)` prior only ratchets upward and suppresses
  nothing — the prior must be a fixed penalty.
- Units are interpretable: for weakly-held transcripts the evidence sum ≈
  the transcript's assigned read count, so `-logit(rho)` ≈ the soft
  effective-exclusive-read floor (default ~3), smooth and self-reinforcing.
- Unique-read invariant: any unique read drives `q → 1` (verified: 0
  suppressed transcripts with unique reads).
- Posteriors in `<prefix>.presence.tsv`; calibration on NA12878-cdna:
  q<0.1 bin is 6.0% truly present, q≈1 bin 94.3%, sharply bimodal.

### Panel B results (vs byte-identical capped baseline, logistic arm)

| sample | presence Δ | floor-0.5 Δ |
|---|---|---|
| nanosim-NA12878-cdna | **+0.0079** | +0.0063 |
| nanosim-H9-cdna | **+0.0054** | +0.0043 |
| nanosim-NA12878-drna | **+0.0053** | +0.0045 |
| nanosim-H9-drna | **+0.0046** | +0.0041 |
| tksm-RSII | −0.0036 | −0.0005 |
| tksm-SQ2 | −0.0024 | −0.0012 |
| **mean** | **+0.0029** | +0.0029 |

- **NanoSim/ONT family: +0.0058 mean — beats the count-floor stopgap on
  every sample and captures ~21% of the +0.0273 detection headroom.**
- **TKSM PacBio family regresses** (−0.002…−0.004; within the ≤0.01
  family-regression rule, but the +0.005 mean gate for a default is NOT met).
- Root cause isolated, not fixable by tuning: on SQ2, 9.6% of suppressed
  transcripts are real and carry ~12 true reads each (86,925 reads total) —
  isoforms whose reads the EM already handed to a dominant sibling, so their
  leave-one-out evidence is genuinely small. Probes: `rho` 0.05→0.35 moves
  SQ2 only −0.0024→−0.0019 (and NA12878 *gains* at 0.2: +0.0082); warmup
  10/50/200 is flat. TKSM truth is denser and more singleton-heavy (76k
  expressed, 21–27% at ≤3 reads, vs NanoSim 30k / 17.5%).

**Disposition: ships opt-in (default `none`).** Mirrors the standing PacBio
anomaly (Lead B); a per-technology default decision needs the expanded
PacBio exact-truth panel (TKSM generation configs — external ask).
Also note pre-existing caveat: bootstraps run presence-off (same class of
inconsistency as novel states; fix both together).

## F1′ — `--novel-from-failures` (Arm N)

Projection-failed reads (a segment outside every annotated exon, or every
candidate eliminated) become locus-level novel-isoform evidence instead of
vanishing: their genomically-overlapped transcript sets join the
ambiguity-component union-find, count toward `--novel-min-locus-reads`, and
inject deterministic constant mass into the locus's novel state (they are
single-candidate reads of it; library-size scaling includes them). The
genome-BAM path gains unprojectable-tracking parity with the reads path.
`unexplained.tsv` gains a `failed_reads` column.

### Holdout-20% evaluation (parity BAM, 4,856 genes losing an expressed isoform)

| | baseline | `--novel-from-failures` |
|---|---|---|
| novel loci | 397 | 1,137 |
| evidence reads | 10,070 flagged | 10,852 flagged + **57,132 failed** |
| gene-level recall | 0.93% (45) | **19.65% (954) — 21×** |
| precision (locus level) | 0.92 | 0.79–0.88 |
| annotated-transcript Spearman | 0.613438 | 0.613045 (Δ −0.0004) |

- Recall target (≥10%) **met with 2× margin**; precision target (≥0.98)
  **not met** — best 0.88 at evidence ≥5 & unexplained-fraction ≥0.02.
  Note the baseline itself scores 0.92 on this harness (its historical 1.000
  came from a differently-constructed holdout), so the feature costs ~4–13pp
  precision for 21× recall.
- Residual FPs are depth-confounded background (0.8% per-read failure rate ×
  deep loci) plus locus-attribution bleed at overlapping genes; the
  `failed_reads` and `unexplained_fraction` columns let consumers set their
  own operating point.
- Annotated quantification is untouched — the novel-state mass injection is
  isolated as designed.

**Disposition: ships opt-in under `--model-unannotated-isoforms`.** This is
the projection track's live feature per the ceiling analysis (Arms S/X cut).

## Follow-ups

1. Expanded PacBio exact-truth panel (external: TKSM generation configs) →
   revisit presence default and Lead B jointly.
2. Presence × novel-from-failures interaction arm (both on, genome mode).
3. Panel A (LongBench/SIRV) sweep for both features before any default flip.
4. Bootstrap consistency (presence + novel states) in one change.
5. Precision work for novel detection: strand-aware overlap and per-locus
   background-rate normalization are the two obvious levers.
