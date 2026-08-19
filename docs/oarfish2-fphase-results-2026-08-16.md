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

## F8 — inference-safe uncertainty infrastructure (2026-08-17)

Response to the post-selection-inference analysis: within-experiment prior
sharing correlates replicates and biases downstream dispersion estimation
(anti-conservative); the design goal is per-sample-independent estimation
with irreducible uncertainty propagated in forms downstream tests already
handle. Three pieces:

1. **Presence/novel-aware inferential replicates.** Bootstrap replicates now
   quantify the SAME model as the point estimate: the resampled M-step
   includes novel latent states, and the presence posterior is refit per
   replicate (multiplicity-weighted leave-one-out pass). Measured on
   NA12878-cdna (30 reps, presence stack): under the old combination
   (presence point + presence-off replicates) **59% of replicates
   contradicted** suppressed point estimates — pure model mismatch. With
   aware replicates the contradiction rate halves to **27.6%**, and the
   remainder is genuine bimodal presence uncertainty (replicates refit q and
   flip at boundary transcripts), i.e. exactly what infRV-aware tests
   (swish/fishpond) are designed to discount. 95% replicate-interval
   coverage of truth on the boundary set: 0.891 vs 0.880. Known remaining
   approximation: projection-failed novel-state base mass is carried
   unresampled.
2. **`--write-identifiability` sidecar.** Per transcript: ambiguity-component
   id and size, unique reads, presence posterior, dominant shadow (the
   θ-argmax co-candidate winning its reads) and shadow share — the structure
   needed for Terminus-style group-level testing at the resolution the data
   identifies. On the parity smoke: 71,220 rows, 21,840 in multi-transcript
   components, 17,507 with a named shadow.
3. **`--presence-prior-file` re-scoped to exogenous priors.** Help text now
   warns explicitly: never build the prior from samples entering the same
   downstream test (within-condition sharing → deflated dispersion,
   anti-conservative; cross-condition → attenuated contrasts, measured at
   6–11% relative spurious-mass inflation on the presence-switch sets).
   Safe use: atlases, pilot samples, reference runs outside the tested
   design. The F7 cross-sample numbers stand as *accuracy* results; the
   feature's inferential use is exogenous-only.

## Leg 2 — SIRV real-read demonstration (2026-08-17)

Real ONT reads (E2 mix, 128-fold concentration range), Lexogen's deliberately
incorrect annotation variants: C (correct, 69 isoforms), I (insufficient, 43
— 26 real isoforms deleted across all 7 genes), O (over-annotated, 100 — 31
fabricated isoforms with certain truth-zero). Genome mode against the 7-gene
SIRV genome, full presence stack vs baseline. Assets:
`oarfish-evaluation-data/sirv-annotations/` (downloaded set1_170612a + runs).

**(a) Over-annotation → presence model, real reads, real truth-zeros.**
Baseline places **9.06% of the library (17,391 reads) on the 31 fabricated
isoforms**; the presence stack cuts this to **5.67%**, drops fabricated
calls (est>0.5) from 18 to 13, and raises Spearman on the 69 real isoforms
**0.7329 → 0.7604 (+0.0275)**. Posterior separation is crisp: fabricated
median q = 0.050 (18/31 below 0.5) vs real median q = 1.000 (2/69 below).
Control: under the correct annotation the presence stack is a perfect no-op
(0.8257 = 0.8257) — every real isoform earns q = 1 and nothing changes.

**(b) Insufficient annotation → novel-from-failures, real reads.**
Under I, **20.8% of aligned reads (41,155) fail projection outright** (all
attributable to annotated loci) — the real-data realization of the
holdout-sim recall pool. Unexplained-mass accounting: junction-flag evidence
alone reports 9.0%; `--novel-from-failures` reports **29.0%, vs a true E2
molar share of the missing isoforms of 32.9%** — 88% of the truly-missing
mass recovered and correctly attributed, on real reads. Caveat: at SIRV
density all transcripts collapse into one ambiguity component, so per-locus
resolution is trivial here (1 locus); fine-grained locus attribution was
demonstrated at genome scale in the holdout evaluation (1,137 loci).

**(c) Poly(A) census on real libraries.** Sequence-visible tails are scarce
in these public SIRV libraries: dRNA 0.98% (2017-era dRNA basecalling rarely
emits the tail — dRNA 3'-anchoring is structural, not sequence-visible),
cDNA 5.05%. The polyA term is honestly near-inert on these particular
libraries; its measured per-tail-read value (TKSM) stands, and modern
cDNA/Kinnex libraries with higher tail retention remain its target.

Net: Leg 2 delivers the two real-read headline demonstrations — presence
suppression with certain truth-zeros (+0.0275 Spearman under
over-annotation, no-op under correct annotation) and failure-driven
unexplained-mass accounting within 4pp of truth under under-annotation.

## Leg 1 — designed TKSM experiment: end-to-end inference results (2026-08-17)

Dataset built per the pre-registered spec (`scripts/leg1-dtu-sim/spec.md`):
2 conditions × 4 replicates × 10M reads, SQ2-like HiFi (identity 99.15, tail
rate 3.4%, lengths matching real SQ2), exact per-read truth (10.0M reads per
sample, 0 unmatched), truth sets: 2,000 DTE / 400+400 presence switches /
400 DTU swaps / ~52k expressed nulls. 104 GB under
`oarfish-evaluation-data/leg1-dtu-sim/`; fully regenerable from committed
scripts + seeds. Both quantification arms (presence stack / plain) with 30
presence-aware bootstraps; converted to salmon format by
`oarfish2salmon.py` (fishpond ≤2.16 has no native oarfish loader — verified
in sources; swap point marked in `run_de.R`).

### Production DE over oarfish inferential replicates (q/FDR ≤ 0.05)

| arm/method | DTE power | pres-off | pres-on | DTU | empirical FDR |
|---|---|---|---|---|---|
| pres / swish | 0.390 | 0.909 | 0.948 | 0.438 | **0.041** |
| pres / edgeR-OD | 0.704 | 0.997 | 1.000 | 0.598 | 0.070 |
| pres / edgeR-naive | 0.603 | 0.992 | 1.000 | 0.554 | 0.058 |
| plain / swish | 0.448 | 0.937 | 0.945 | 0.455 | **0.047** |
| plain / edgeR-OD | 0.704 | 0.997 | 1.000 | 0.599 | 0.068 |
| plain / edgeR-naive | 0.605 | 0.995 | 1.000 | 0.554 | 0.058 |
| shared-prior / swish | 0.393 | 0.914 | 0.945 | 0.440 | 0.043 |
| shared-prior / edgeR-OD | 0.704 | 0.997 | 1.000 | 0.598 | 0.070 |

Findings:

1. **The full uncertainty chain works end-to-end with production tools**:
   swish over oarfish presence-aware inferential replicates **controls FDR**
   (0.041–0.047 at nominal 0.05) with presence-switch power 0.91–0.95 —
   validating the replicates, the converter, and the whole stack in one shot.
2. **edgeR v4 + catchSalmon overdispersion**: highest power (DTE 0.70, DTU
   0.60, switches ≈1.00) but mildly anti-conservative on this design (FDR
   0.068–0.070, vs 0.058 naive) — the moderation buys power, not calibration,
   here.
3. **Presence arm vs plain downstream**: a small conservative shift in swish
   (FDR 0.041 vs 0.047; DTE power 0.390 vs 0.448) and no difference for
   edgeR-OD. Per-sample accuracy shows the known PacBio-family presence
   regression (−0.0021 Spearman, uniform across all 8 samples) — this
   SQ2-derived design inherits Panel B's TKSM behavior, reinforcing the
   per-technology opt-in disposition.
4. **The post-selection experiment**: within-condition prior sharing at
   w = 0.8 (each replicate primed by a same-condition sibling) does **not**
   measurably inflate FDR (swish 0.043 vs 0.041; edgeR unchanged to 3
   decimals). The theoretical anti-conservativeness is empirically bounded
   near zero under the shipped floor design (protective-only, one-directional
   priors; both conditions treated identically; swish's infRV moderation
   defending the boundary set). Together with the measured cross-condition
   contrast attenuation (6–11% on switch sets), the full picture: the
   exogenous-only recommendation stands on principle, but the feature's
   failure modes are measured small in both directions.

## F3 — PacBio physical endpoint model resurrected and gated (2026-08-17; Lead B closed)

Ported from the archive as `--coverage-model pacbio-endpoint` (bare kernel:
per-read endpoint simplex in the coverage substrate; containment guard
applied at the source — nested-short candidates neutralize the read's
endpoint term — instead of the archived post-EM clamp). Evaluated on four
exact-truth PacBio samples (TKSM RSII/SQ2 + Leg-1 A1/B1):

| sample | none | logistic | pacbio-endpoint |
|---|---|---|---|
| RSII (CLR) | 0.9328 | **0.9392** | 0.9310 |
| SQ2 (HiFi) | **0.9501** | 0.9481 | 0.9489 |
| Leg1-A1 (HiFi) | **0.9411** | 0.9401 | 0.9406 |
| Leg1-B1 (HiFi) | **0.9383** | 0.9373 | 0.9382 |

Pre-registered gate (beat BOTH none and logistic by ≥ +0.002): **not met.**
The kernel edges logistic on every HiFi sample (+0.0005…+0.0009) but never
beats `none`; on CLR it is harmful (learned mixture: 10.7% intact, 41%
3'-truncated, tolerance at clamp — endpoint geometry is noise there).

**Lead B disposition**: the reproducible "PacBio anomaly" is now explained on
exact truth — it was never a missing better coverage model. On HiFi,
positional/endpoint evidence simply does not pay; `--coverage-model none` is
the best HiFi configuration (best arm on all three HiFi samples). Practical
recommendation: per-technology coverage default of `none` for pac-bio-hifi,
`logistic` for CLR/ONT. The flag ships opt-in with this documented negative.

## Over-annotation at genome scale: mechanism isolated (2026-08-18)

A genome-scale SIRV-O analog (13,361 decoy isoforms — skip/alt-terminus/
extend perturbations at every eligible expressed gene; `make_decoys.py`) was
quantified against by the Leg-1 samples plus two single-variable arms:

| regime | mean read | decoy leak (%lib) | decoys called | presence effect |
|---|---|---|---|---|
| HiFi (99.2% id, KDE truncation) | 1,580 nt | 0.12% | ~380 | calls −10%, mass ≈0 |
| ONT-error (96% id, same molecules) | ~1,570 nt | 0.12% | ~320 | ≈ none |
| ONT-error + heavy truncation | 726 nt | **12.4%** | 3,397 | **≈ none** |
| (real dRNA SIRV-O anchor) | short/anchored | 9.1% | 18/31 | mass −37%, +0.0275 Sp |

Findings:

1. **Read completeness — not error rate — is the over-annotation
   vulnerability.** Same molecules at 96% vs 99% identity: identical 0.12%
   leak. Same pipeline with truncation forced to ~726 nt: 12.4%. Long reads
   that span discriminating features are structurally immune to annotation
   decoys; truncated reads are not. (Corollary: HiFi/Kinnex users are
   largely protected; dRNA/degraded-cDNA users are the exposed population.)
2. **Downstream robustness on HiFi**: even at 0.12% leak, over-annotation
   costs the plain pipeline ~6pp of swish DTE power (0.448 → 0.386) while
   the presence arm is unaffected (0.390 → 0.388); FDR controlled in both
   (0.043).
3. **A measured presence-model boundary**: under heavy truncation at this
   decoy density, each decoy equilibrates at ~360 assigned reads — far above
   the effective-exclusive-read floor — and the leave-one-out evidence
   genuinely supports it: presence suppresses 0% of mass-bearing decoys
   (vs 58% of fabricated SIRV isoforms, whose equilibrium mass stayed near
   the floor at spike-in depth). Individually-ambiguous truncated reads are
   the sibling-shadowing information bound in another form; suppressing
   these decoys needs evidence channels beyond assignment mass — the
   poly(A)/full-length term (absent in this sim's surviving tails),
   positional models robust to decoy pollution of the unique-read training
   set, or annotation-level structural priors. Recorded as the open problem
   the truncated-protocol population needs solved.

## Anchored 3' likelihood for dRNA (2026-08-18) — the needle-mover, mapped end to end

**Premise (four-corner measurement).** Real dRNA is strongly 3'-anchored:
67.4% of unique SIRV-dRNA reads abut the annotated 3' end in alignment
coordinates, **91.7% crediting the molecule-3'-side soft clip** (the ragged
read start is clipped); mean gap 20 nt. NanoSim dRNA does NOT reproduce this
(29.5%/34.5%, mean gap 550 nt; its cDNA looks the same) — resolving the
apparent "missing 3' signal" as simulator infidelity, disqualifying Panel B
dRNA for gating anchored/endpoint models, and retroactively explaining part
of the July endpoint-model burial.

**Feature** (`--anchor-three-prime`): every dRNA read treated as
3'-complete; per-sample empirical likelihood over a signed, clip-aware 3'
mismatch — undershoot (annotated end beyond the read after clip credit) and
overhang (clip beyond an 80 nt poly(A)/adapter allowance: the signature of a
candidate shorter than the molecule) both penalized, odds capped at 20.

**Real-data results (SIRV E2 dRNA, alignment mode):**

| | baseline | anchored (signed) |
|---|---|---|
| correct annotation, Spearman vs E2 | 0.8175 | **0.8896 (+0.072)** |
| over-annotated, Spearman (69 real) | 0.7449 | **0.8873 (+0.142)** |
| over-annotated, fabricated mass | 9.22% of lib | **0.51% (−94%)** |

Presence alone: no effect here (fabricated isoforms equilibrate far above
any evidence floor — assignment-level siphoning requires an
assignment-level likelihood). Largest single-feature gains measured in the
program, on real reads with known concentrations.

**Exact-truth boundary (anchored-truncation sims, 10M reads each).** A
dRNA-faithful simulator was built (custom 5'-side mdf truncation — tksm's
`--always-end` anchors the wrong end; fidelity: 99.9% flush). At median
600–1,300 nt anchored reads, baseline accuracy collapses (Spearman
0.34–0.42) and decoy leak *rises* to 19–22%: anchoring stacks every read on
the transcriptome's most-shared region. Leak anatomy: **alt-3' decoys
0.00%** (the aligner's scoring already eliminates overhang candidates);
the leak is entirely alt-5' (13–14%) and exon-skip (5–7%) decoys — 3'-
sharing structures that 3'-anchored fragments *definitionally* cannot
distinguish by 3' evidence, and whose 5' side is degradation-confounded.
The anchored term is honest here: +0.003, no harm.

**Disposition**: ships as the dRNA headline feature with a precisely stated
domain — large gains whenever candidate structures differ within reach of
the read (terminal variants, structurally distinct false isoforms; all of
SIRV's fabrications), structurally silent for 3'-sharing ambiguity in
heavily truncated data, where the correct response is the identifiability/
grouping machinery, not more evidence. Evaluation used real SIRV dRNA
(premise + gains) and the anchored sims (boundary); NanoSim dRNA must not
be used for either.

## Anchored-3' addendum: human calibration flips the sign (2026-08-18)

Census of **real human dRNA** (SG-NEx MCF7, HEK; RefSeq frame): only
18–25% of unique reads are 3'-flush with the annotated end (21–45%
clip-inclusive), mean gap ~0.9 kb — versus SIRV's 67%/92%/20 nt. Coherent
mechanism: proliferative lines favor proximal APA; molecules genuinely end
far upstream of RefSeq's distal annotated termini. The dRNA anchor premise
splits: "read end = molecule end" holds (poly(A) selection), but
"molecule end = annotated end" does not on human.

A sim calibrated to the human census (jitter-p0 0.3, exponential 1.2 kb
offsets; gate passed: 28.2% flush, 644 nt mean gap) gives the exact-truth
verdict: the current global-gap anchored likelihood is **harmful**
(std −0.0036, decoy arm −0.0042) — with APA, gap differences encode
polyadenylation variation, not read origin.

**Scoped disposition**: `--anchor-three-prime` is validated ONLY where 3'
termini are annotation-exact (spike-ins; curated end sets). Help text to
be updated accordingly. The SIRV result stands as the mechanism
demonstration; the human-scale version of this feature is the identified
upgrade: **per-transcript empirical 3'-end profiles** learned from each
sample's own unique reads (APA becomes per-transcript signal instead of
global noise), cross-fitted, scored per candidate. Both the calibrated sim
(A1c) and SIRV now exist as its gates. Also recorded: the four-corner
census methodology itself (juncprobe tail-probe --assume-anchored) as the
per-sample premise check any anchored model should run.

### 3'-discriminability partition of human ambiguity (A1c, RefSeq)

74.7% of reads are multi-candidate; among those, the spread of candidate
annotated-3'-gaps at the read's 3' position exceeds 100 nt for only
**21.7%** (23.7% @50 nt, 17.9% @200 nt). So the pool any 3'-end model can
address on human RefSeq is ~16% of all reads (0.747 x 0.217); ~78% of
ambiguity is 3'-SHARING — isoforms with a common terminus differing
internally/5'. This bounds the human upside of even a perfect (per-
transcript, APA-aware) 3' model well below the SIRV demonstration, where
isoforms differ mainly at their ends and termini are annotation-exact.
Combined ranking of remaining leads: per-transcript 3'-end profiles are a
scoped, bounded win (~16% addressable pool, needs APA-aware form);
presence/absence and projection tracks remain the larger levers.

### Census correction + mechanism diagnostics (2026-08-19, after review)

Challenge raised: the census conditioned on unique reads (tail-probe's flush
stat is computed only for single-candidate reads), which could deplete flush
reads since ~78% of ambiguity is 3'-sharing. Recomputed as per-read MIN
clip-aware gap over ALL candidates: MCF7 45.3->49.7% flush, HEK 20.6->19.6%;
ambiguous reads no more flush than unique (51.0% vs 49.7%). Bias real in
design, small in effect: half of MCF7 / 80% of HEK reads are flush with NO
annotated terminus.

Cause diagnostics (both samples): internal A-stretch capture DEAD
(downstream-15nt A-frac 0.25-0.27 = background); degraded-read artifact DEAD
(flushness flat across length strata); antisense DEAD (0.5% reverse).
POSITIVE: large-gap read 3' ends cluster at discrete internal sites — on
MCF7, modal +/-25nt site captures 73% of such reads per transcript
(HEK 45%) = unannotated poly(A) sites / short-isoform ends. Combined with
the SIRV control (same protocol, 91.7% flush), the protocol is fine:
"read end = molecule end" holds; "molecule end = annotated end" fails (APA).

Sign-flip mechanism: when a molecule ends at a proximal APA site, the TRUE
isoform shows a ~1kb gap while a shorter sibling whose annotated end
coincides with that site shows gap 0 — the global-gap likelihood then
confidently reassigns the read to the wrong flush candidate. Evidence
inverted, not diluted, exactly for APA-affected reads.

Recalibrated sim A1d matched to the corrected census (jitter-p0 0.5, mean
1400; gate: 53.4% min-gap flush vs MCF7 49.7%): std −0.0005, decoy −0.0020,
decoy mass +0.33pp. Softer than A1c's −0.004 but still no benefit. Verdict
robust across both calibrations: on APA-bearing human data the global-gap
anchored likelihood is neutral-to-harmful; the SIRV win is real but
annotation-exact-only. Upgrade path unchanged (per-transcript 3'-end
profiles; the 73% modal clustering is direct evidence the per-transcript
end distribution is learnable).

## Bottom-up feature census (2026-08-19)

Question inverted per review: not "which mechanisms should help" but "which
features in the data actually discriminate." Two probes: (1) exact-truth
placement census on A1d (2.5M sampled reads; truth vs best-AS winner, per
feature, overall and AS-tied); (2) real-data end-structure on MCF7 dRNA.

DEAD (measured, third dataset for the tie result):
- 94.1% of misassignments are AS-tied; among tied, EVERY per-alignment
  geometry feature is ~0 or ANTI-informative: lendiff/txpcov/txplen net
  -19%, gap5 -16%, startf -15% (short isoforms covering the read span look
  "cleaner"; truth is systematically the longer source). Naive plausibility
  likelihoods would hurt. gap3c ~0 on jittered ends (consistent w/ A1d quant).
- Transcript-level coverage-shape aggregates for FP detection: mean_qcov /
  flush_frac / end-dispersion AUC 0.51-0.52 = dead. FP separation lives in
  unique_frac (AUC 0.675) and n_reads (0.625) — the presence-track variables.

ALIVE (real data only; invisible in sims by construction):
- 3' PAS support: AATAAA/ATTAAA in [-40,-5] of read 3' ends: 53.9% (flush),
  32.6% (large-gap), 2.5% internal control. Read 3' ends are genuine pA
  sites; truth-free per-candidate feature.
- 5' structure: 26.0% of reads clip-inclusive-flush with an annotated start
  (~26x uniform-truncation expectation) — bulk full-length reads exist; and
  non-flush 5' truncation points cluster (modal +/-25nt = 48.8%/txp).

Synthesis: per-READ geometry is exhausted; the live signal is per-TRANSCRIPT
end-position structure at BOTH ends, concentrated and learnable from unique
reads. Points to a two-sided per-transcript end-profile likelihood (5' side
addresses the 78% 3'-sharing ambiguity pool that the 3' side cannot), with
the census caveat baked in: only sample-calibrated profiles are justified —
naive 5'-flush preference is anti-informative on truth. Sims need
site-structured end placement (PAS-anchored 3', clustered 5') before they
can evaluate any of this; current A1d cannot express these signals.

### End-profile evaluation verdict (2026-08-19): FAILS the quantification gate

Step 1 (held-out predictability, real MCF7/HEK): PASSED with margin
(+1.1-1.3 nats/read vs global, positive in every coverage stratum).
Step 3 results (--end-profile vs base vs anchored):

| dataset | arm | Spearman | decoy mass |
|---|---|---|---|
| SIRV E2 real (C) | base / anchor / endprof | .8175 / .8896 / .8898 | — |
| SIRV E2 real (O) | base / anchor / endprof | .7449 / .8873 / .8884 | 9.22 / 0.51 / 5.34% |
| A1d control      | base / anchor / endprof | .4123 / .4118 / .4116 | 21.9 / 22.2 / 22.5% |
| A1e site-struct  | base / anchor / endprof | .5473 / .5457 / .5440 | 16.8 / 17.1 / 17.8% |

A1e (per-transcript APA sites + 5' hotspots, census gate passed: 53.5%
flush) is the favorable case by construction — sites drawn INDEPENDENTLY
per transcript maximize sibling discriminability, more than real shared-PAS
biology would — and end-profile still loses (−0.0033 std, +1.05pp decoy
mass). Pre-registered bar (≥ +0.002) decisively missed.

Mechanism of failure (why predictive ≠ discriminative):
1. A concentrated profile helps a transcript's site-reads but PENALIZES the
   same transcript's non-site tail reads (~half of them) relative to an
   unexpressed sibling's broad global-shrunk profile — the two effects
   cancel to slightly negative.
2. Decoy forgiveness: decoys that capture a few unique reads learn
   self-excusing profiles (SIRV fabricated mass 5.34% vs anchored's 0.51%);
   shrinkage bounds but does not eliminate this.
3. The EM's abundance term already resolves ties in the direction end
   evidence would mostly confirm; where they disagree, end geometry is
   more often wrong (consistent with the feature census).

Disposition: --end-profile stays as an experimental flag with this verdict
in its provenance; NOT a headline feature. What survives from this line:
(a) --anchor-three-prime for annotation-exact termini (unchanged), (b) the
per-transcript end-profile machinery as the basis for the BANKED APA
reporting deliverable (profiles as a per-sample APA readout — the held-out
validation shows they are accurate descriptions; they are just not useful
as assignment evidence), (c) the corrected census + site-structured sim as
standing methodology. Program priority returns to presence/absence and
projection tracks.

## F4.1 end-consistency presence evidence: GATE PASSED (2026-08-19)

Global-error decomposition (Panel B logistic baseline): FP transcripts carry
0.07-0.28% of library mass but cost +0.015..+0.049 Spearman and 0.03-0.07
MARD (oracle removal); FN mass 0.2-0.6%; detected-only Spearman 0.96-0.98.
Presence is decisively the global-metric lever; magnitude error secondary
(except SIRV-C, where anchored already recovers logP 0.75->0.90).

F4.1 measurement — per-transcript end-consistency vs count baselines,
FP-vs-real AUC (est>=1, >=3 reads; zone = no unique reads):

| statistic | SIRV O real (zone) | A1e site-struct (zone) |
|---|---|---|
| flush5_frac (reads reaching OWN annotated start) | 0.974 (0.960) | 0.729 (0.748) |
| flush3_frac | 0.889 (0.860) | 0.521 (0.525) |
| mode_capture | 0.671 (0.703) | 0.530 (0.583) |
| unique_frac (baseline) | 0.729 | 0.631 |
| n_reads (baseline) | 0.524 (0.600) | 0.513 (0.523) |

Pre-registered gate (in-zone AUC >= 0.70, must add over unique_frac):
**flush5_frac passes on both datasets** — 0.960 real SIRV, 0.748 sim — and
is strongest exactly where count evidence dies. Mechanistic irony: the 5'
side, useless per-read, is the top transcript-level FP discriminator — a
real transcript earns its measured full-length excess (26% of reads reach
its start); a fabricated isoform's misassigned reads almost never reach the
FP's own start. flush3 is diluted wherever 3' ends are shared (A1e decoys)
but adds on real SIRV.

Next (F4.2): binomial evidence log-LR on (k flush-5' of n reads; per-sample
full-length rate vs homology background) into PresenceState::update_q via
the existing --presence-endpoint-alpha hook; evaluate presence vs
presence+end on Panel B + SIRV O + A1e per pre-registered F4.3 bars.

### F4.3 verdict (2026-08-19): integration FAILS despite the F4.1 AUC pass

Arms (presence vs presence + flush evidence, alpha 1; then cap-fixed alpha
1/3 on dev sets):
- Panel B: presflush −0.001..−0.002 Spearman on all 6 samples (bar was
  >= +0.002) — evidence engaged even on NanoSim (self-gate too permissive
  there) and added noise.
- SIRV O: fabricated mass 9.22% -> 9.22% (alpha 1, capped) -> 8.63%
  (alpha 3, uncapped) while real-transcript Spearman fell 0.745 -> 0.705.
- A1e: FP mass 23.65% -> 23.64%; Spearman −0.002..−0.004.

Root cause, two layers:
1. Implementation (fixed): the flush log-LR was capped at ±25 while the
   leave-one-out delta grows per-read — mass-carrying FPs were
   arithmetically immune. Uncapping (TERM_CAP 1e4) is the principled form.
2. Fundamental (the kill): FP MASS lives in phantoms that SHARE their 5'
   start with a real isoform — their flush counts look normal, so end
   evidence cannot indict them — while the AUC 0.96/0.75 was driven by the
   many low-mass phantoms that count evidence already handles. Meanwhile
   real transcripts with weak 5' representation take collateral damage.
   Mass-weighted discrimination, not ROC, is the right gate for presence
   evidence; recorded as a methodology lesson.

Consistent with the F2 finding, now measured from a second direction:
mass-carrying fabrications sit at an EM equilibrium that per-transcript
posterior evidence cannot dislodge — only assignment-level likelihoods
reach them (anchored: 0.51% vs presence-channel floor ~8.6%).

What stands after F4: (a) presence-alone captures part of the decomposed FP
headroom (Panel B logistic 0.8871 -> pres 0.8951 on NA12878-cdna; +0.005..
+0.008 across samples); (b) the remaining oracle gap (to ~0.936) is in
low-count FPs vs weakly-expressed reals — flush evidence was the candidate
and failed; (c) --presence-flush-alpha stays experimental, default 0, with
this provenance. F4 track closed.
