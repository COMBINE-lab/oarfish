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

## F4 — positional plausibility & boundary endpoints (added later on 2026-08-16)

Question posed: can a coverage-plausibility model supply the "new evidence"
for the overlap zone, and does penalizing alignments terminating at internal
exon boundaries (especially sibling-terminal ones) help?

**F4.1 measurement (juncprobe `zone-sep`, all six Panel B samples).** In the
overlap zone (no unique reads, 0 < est ≤ 5):

- Boundary-coincidence statistics are **uninformative**: endpoint-within-W-of-
  internal-boundary AUC 0.51–0.52; the sibling-terminal-conditioned variant
  0.49–0.54 even restricted to transcripts where it can fire. The archived
  per-read LR (9.9×) does not survive per-transcript aggregation in the zone.
- **Endpoint-gap surprisal does separate**: reads on false-positive
  transcripts sit at atypical positions under a length-stratified positional
  model fitted from unique reads (AUC 0.72 left-gap alone). A 5-feature
  logistic model (surprisal L/R, log length, candidate density, soft reads)
  fit on NA12878-cdna and **frozen** scores AUC 0.76–0.85 on all six samples
  (est ≤ 1 sub-zone: 0.55–0.71). FPs also skew long (length AUC 0.74).
- Limit found: presence-suppressed transcripts with est ≈ 0 have no assigned
  reads, hence no positional evidence — post-EM positional statistics cannot
  rescue the TKSM suppressed-real set (frozen-score AUC 0.52 there).

**F4.2/F4.3 integration**: `--presence-endpoint-alpha` — each assigned read
adds its responsibility-weighted log-LR of sitting where unique reads sit
(vs uniform) to the transcript's presence evidence. Panel B, logistic arm,
α = 1:

| family | presence Δ vs baseline | presence+endpoint Δ vs baseline | endpoint Δ vs presence |
|---|---|---|---|
| NanoSim (4) | +0.0058 | **+0.0062** | +0.0004 |
| TKSM (2) | −0.0029 | **−0.0020** | +0.0011 |
| all (6) | +0.0029 | **+0.0035** | +0.0007 |

Improves **every sample** (6/6) over presence-alone and shrinks the PacBio
regression (α = 2 helps TKSM further but costs ONT; α = 1 is the balanced
point). Against the pre-registered gates: TKSM-shrink met; the NanoSim
≥ +0.002-over-presence bar not met (+0.0004) — consistency without magnitude,
so the term stays **opt-in** alongside the presence flag it extends.

Net position on the ONT/NanoSim family: baseline 0.9136 → 0.9198 with
presence+endpoint, ~23% of the +0.0273 detection ceiling captured.

## F5 — pre-EM candidate-level evidence (2026-08-16; gate not met)

Question: can positional/tiling evidence computed over *candidate sets*
(θ-independent, uniform-prior responsibilities) rescue the transcripts the
EM zeroes before any post-EM evidence exists — the TKSM suppressed-real set?

`zone-sep --resp-mode uniform --include-zero` adds candidate-level coverage
**breadth** (occupied fraction of 20 coarse bins) and **occupancy entropy**,
the tiling signature: a real transcript's candidate reads span it (different
reads shared with different competitors) while a false positive's borrowed
reads concentrate in its dominant sibling's shared region.

Measured on the zeroed sets (est ≤ 0.01, ≥2 candidate reads; ~34k/45k
transcripts, ~7–10% real):

- Best single features: breadth 0.59–0.63, entropy 0.57–0.66 AUC — the
  hypothesized mechanism is real but weak. Positional surprisal *inverts* at
  candidate level (0.40–0.46): uniform responsibilities dilute every
  candidate with the sibling's reads.
- Frozen cross-sample combination: 0.672 (cdna, train) / **0.635 (SQ2)**;
  on the full SQ2 presence-suppressed set (73k of 76k now covered, vs 1.7k
  post-EM): **0.639**. Depth-stratified: cdna climbs to ~0.70 at ≥30
  candidate reads, SQ2 flat at ~0.65 at any depth.

**Gate (AUC ≥ 0.70) not met — not integrated.** The deadlock is now measured
from both sides: post-EM the shadowed transcript has no reads (no evidence);
pre-EM its profile is a mixture dominated by the sibling's reads (diluted
evidence). Absolute tiling statistics cannot break it.

The one untested formulation with a mechanistic reason to be sharper:
**differential/competitor-aware tiling** (S4 proper) — does t's candidate
set contain reads at positions its dominant competitor *cannot* explain,
computed by projecting exon structures between candidate pairs. That is a
substantially bigger build (pairwise structure projection) and is recorded
as the follow-up, not attempted here.
