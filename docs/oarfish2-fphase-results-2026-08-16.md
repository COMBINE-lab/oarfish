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

## F5b — differential (competitor-aware) shadowing evidence (2026-08-16; dead)

The follow-up formulation from F5, with a simplification that also bounds the
exon-projection variant: "reads at positions the dominant competitor cannot
explain" ⟺ "reads whose candidate set excludes the competitor" — the aligner
itself adjudicates explainability, and reads are the only observations of
positions. So `zone-sep --differential` computes, per transcript: its
dominant shadow (θ-argmax winner over its candidate reads), `d1_frac` (the
candidate mass from reads excluding that shadow), the breadth of those
non-shadow reads, the shadow's mass share, and the winner entropy
(sharing-diversity).

Measured (uniform responsibilities, cdna + SQ2):

| feature | zeroed sets | SQ2 suppressed | est(0,5] zone |
|---|---|---|---|
| d1_frac | 0.51–0.54 | 0.53 | 0.44–0.49 |
| ns_breadth | 0.51–0.54 | 0.53 | 0.48–0.54 |
| winner_entropy | 0.52–0.53 | 0.50 | 0.50–0.52 |
| shadow_share | 0.47–0.48 | 0.50 | 0.51–0.56 |

Chance-level everywhere; adding them *degrades* the frozen cross-sample
model (SQ2 zeroed 0.635 → 0.586). Interpretation: in a dense annotation an
absent transcript's borrowed reads also arrive from several expressed
neighbors, so its sharing-diversity matches a real weak transcript's — the
"one specific sibling" premise fails. Since read-set containment upper-bounds
what pairwise exon-structure projection could see, the differential approach
is closed, not deferred.

**The information boundary is now measured from three sides** for the
fully-shadowed set: post-EM positional (no assigned reads), pre-EM absolute
tiling (sibling-diluted, AUC ≤ 0.67), pre-EM differential/containment
(chance). Within one sample's read-level data, real-but-fully-shadowed and
absent are indistinguishable — which is precisely why the EM shadowed them.
Breaking this tie requires *external* evidence: cross-sample joint priors
(the same transcript observed unshadowed in a related sample), hybrid
short-read priors (`--short-quant` already exists as the entry point), or
richer annotation-independent signals (e.g. full-length flags from adapters).

## F6 — protocol-level full-length (poly(A)) evidence (2026-08-16; shipped opt-in, resurrects a "killed" idea)

Mechanism: a poly(A) tail in a read's terminal soft clip anchors the
*molecule's* 3' end at the alignment end. Among score-tied candidates this is
extra information (it lives in the clip): a candidate whose annotated 3'
terminus is flush with the read's end explains the tail; one that would
place the molecule's end mid-transcript requires a templated internal
poly(A).

**Census first (juncprobe `tail-probe`)** — two findings that rewrite the
history of this idea:

1. **NanoSim simulates no tails** (0.002% detection): the archived
   `polya-three-prime` term (killed at +0.000108) was evaluated on a panel
   where four of six samples had no signal *by construction*. Its "killed"
   verdict was a simulator artifact, not a measurement of the mechanism.
2. **TKSM per-read truth is recoverable** (FASTQ `molecule_id=` comment →
   mdf `tid=`): maps built for SQ2/RSII (14.0M reads each, 0 missing) at
   `sim-panel/tksm/read_truth/`. This unlocks read-level truth on PacBio —
   the handoff believed only transcript counts existed.

**Tie-breaking measurement** (TKSM, per-read truth): among score-tied
misassigned reads with tails, the strong configuration (truth flush ≤16nt &
winner conflict >64nt) favors the truth **1,557:0 (SQ2)** and **4,646:0
(RSII)** — zero counterexamples in 6,203 firings; the soft any-gap version
runs 72–830:1. One-directional, unlike every prior endpoint feature. The
correctable set is small (0.13–0.45% of tied misassignments), and the
per-transcript presence consumer is weak (AUC 0.50–0.59; tail reads' gaps
rarely differ between candidates).

**Integration**: the archived `polya_probability.rs` ported as opt-in
`--polya-three-prime` (per-read conditional 3'-completeness likelihood,
trained on unique tail reads, odds capped at 20; works with or without the
coverage model; skips paths without tail info; no-op on tail-free data).

| sample | logistic stack Δ | presence stack Δ |
|---|---|---|
| tksm-SQ2 (7.6% tails) | +0.00026 | +0.00032 |
| tksm-RSII (23.8% tails) | +0.00060 | +0.00050 |

Positive in all four measurable arms, magnitude ∝ tail rate; NanoSim arms
are structural no-ops. Real libraries carry far higher tail/3'-anchoring
rates than these sims (dRNA ~100%, cDNA with retained poly(A) well above
24%), so the sim-measured magnitude is a floor — but confirming that needs
real data with truth (SIRV read mode, or Panel A as a sign check).

**Killed-list amendment**: `polya-three-prime` moves from "killed" to
"shipped opt-in; prior verdict was a tails-free-panel artifact."

## F7 — multi-sample presence prior (2026-08-16; shipped opt-in, positive everywhere tested)

The external-evidence tie-breaker for the fully-shadowed set:
`--presence-prior-file <other-sample .presence.tsv>` adopts a related
sample's presence posteriors as this sample's per-transcript prior,
`rho_t = w·q_other + (1−w)·rho0` (`--presence-prior-weight`, default 0.5).
Floor design: support elsewhere protects; absence elsewhere never suppresses
below the baseline prior — so a broken sharing premise can only mute the
feature, not weaponize it.

**Panel limitation, measured first**: the simulated sibling samples do NOT
share truth (21–56% of expressed transcripts are sample-specific — each sim
profile was drawn independently), and sibling unique-read support is
near-chance (AUC 0.48–0.54) as a real-vs-FP separator on this panel. So
Panel B *understates* the mechanism; two validations bracket it:

1. **Split-half (exactly shared truth)** — NA12878-cdna split by read-name
   parity; half-1's posteriors prime half-2:
   alone 0.838836 → w=0.5 **+0.0028** → w=0.8 **+0.0042**.
   Clean upper-regime validation of the machinery.
2. **Cross-protocol (broken premise, 21–48% sample-specific)** — cdna↔drna
   within cell line, on top of the full presence+endpoint stack:

   | sample | pres+ep | +cross-prior (w=0.5) | Δ |
   |---|---|---|---|
   | NA12878-cdna | 0.895447 | 0.896595 | +0.0011 |
   | H9-cdna | 0.916759 | 0.917642 | +0.0009 |
   | NA12878-drna | 0.925941 | 0.926231 | +0.0003 |
   | H9-drna | 0.941017 | 0.941598 | +0.0006 |

   Positive 4/4 even under heavy premise violation (w=0.8 slightly better
   again: +0.0014 on NA12878-cdna). Real biological replicates share far
   more expression than these sims, so the split-half number is the better
   estimate of real-data value.

**Cumulative ONT/NanoSim position**: baseline 0.9136 → presence + endpoint +
cross-prior **0.9205** (+0.0069, ~25% of the +0.0273 detection ceiling), all
opt-in, each component positive on every sample it touches.
