# Re-benchmark of the `auto` coverage stack against `logistic` (2026-08-03)

## Summary

Every arm, including `logistic`, was run on both panels specified in
`docs/coverage-auto-rebenchmark-plan-2026-08-03.md`. 448 runs, zero failures.

The headline results:

- **The demotion of rank blending and dominance pruning is confirmed.** Both
  fail decision rule 1 decisively, and they fail it on the panel with exact
  truth. Rank blending costs up to **-0.0878** Spearman on a dataset family and
  dominance pruning up to **-0.0437**. Neither should be an `auto` default.
- **`auto` should not be a default at all.** `adaptive-bare` does not beat
  `logistic` on Panel B (-0.0003, i.e. a tie), and the full `auto-new` default
  is also at parity (-0.0007) while costing **23% more wall time** than
  `logistic` on the same branch. This satisfies decision rule 4.
- **Alignment calibration and censoring pass rule 1 on the letter, but their
  effects are below the noise floor of the experiment.** The same `--model-coverage`
  differs by up to 0.0053 Spearman between `main` and the branch under test;
  the measured effects of calibration (+0.0009) and censoring (+0.0042) are of
  that order or smaller.
- **The two panels disagree completely and systematically**, and this is the
  most important structural finding. Decision rule 5 is invoked: Panel B is
  preferred, and the disagreement is reported rather than averaged away.

## Validation of the measurement path

Before any conclusion, the pipeline was checked against numbers produced by a
different harness. It reproduces every add-one ablation figure in
`docs/coverage-auto-defaults-review-2026-08-03.md` on NanoSim NA12878 1D-cDNA:

| arm | this run (Spearman / MARD) | review doc |
|---|---|---|
| `none` | 0.8392 / 0.0873 | 0.8391 / 0.0873 |
| `logistic` (dev) | 0.8867 / 0.0582 | 0.8867 / — |
| `logistic` (main) | 0.8850 | 0.8850 / 0.0590 |
| `adaptive-bare` | 0.8846 / 0.0593 | 0.8846 / 0.0593 |
| `auto-new` | 0.8824 / 0.0600 | 0.8824 / 0.0600 |
| `auto-old` | 0.7677 / 0.1576 | 0.7677 / 0.1576 |
| `+calib` | 0.8842 / 0.0594 | 0.8842 / 0.0594 |
| `+censor` | 0.8827 / 0.0599 | 0.8827 / 0.0599 |
| `+prune` | 0.8564 / 0.0713 | 0.8564 / 0.0713 |
| `+rank` | 0.7879 / 0.1426 | 0.7879 / 0.1426 |

RMSE for `none` matches to 227.4453 vs 227.4457. Kendall tau-b was
cross-validated against scipy over 200 randomised tie-heavy trials (0
mismatches).

## Method

Nine arms as specified in the plan, `--filter-group no-filters` on every run,
`--seq-tech` passed on every run (`bulk.rs` gates the PacBio physical-endpoint
kernel on it, and the original selection harness passed it too). Metrics are
Spearman, Kendall tau-b, Pearson on log1p, CCC on log1p, RMSE and MARD, with
wall time and peak RSS from `/usr/bin/time -v`.

**Identical reference transcriptome per comparison.** Each sample's evaluation
universe is the reference transcriptome shared by the run and its truth, and
the runner hashes every arm's `.quant` key set and aborts if two compared arms
differ. All 28 LongBench samples resolve to one identical 251,488-transcript
reference (`md5 9afd4f42be9d`); Panel B resolves to 177,816. This matters: the
naive union of truth and estimate keys would have scored ~134,000 transcripts
that the matched-Illumina reference never contained as truth = 0, penalising
arms for an annotation gap rather than a quantification error.

**Truth handling.** Accession versions are stripped (`\.\d+$`) before an outer
join filled with 0. Two format details the plan did not anticipate: TKSM truth
files carry a header row and versioned accessions (NanoSim's are headerless and
unversioned), and SIRV truth is a concentration table whose column is selected
by the manifest `mix`. Panel A estimates are rescaled to the truth total
(required: matched Illumina is a different protocol at different depth); Panel B
estimates are not (truth and estimate are both simulated read counts).

**SIRV E0 is equal-molar.** Every spike-in has concentration exactly 1, so truth
variance is 0 and Spearman, Kendall, Pearson and CCC are all undefined for those
six samples. They are reported as `NA`, not 0.0 — scoring them as 0.0 would have
averaged six undefined samples into every arm's panel mean and read as a large
false regression. For E0 the informative metrics are MARD and estimate
dispersion.

## Panels

| | Panel A (LongBench) | Panel A (SIRV/spike-in) | Panel B (simulations) |
|---|---|---|---|
| samples | 28 | 12 | 6 |
| truth | matched Illumina | SIRV concentrations; 2x exact counts | exact simulated counts |
| universe | 251,488 / 71,220 | 69 / 385,659 | 177,816 |
| absolute Spearman | 0.21-0.43 | 0.82-0.99 (E0 undefined) | 0.84-0.95 |

## Results: Panel B (exact truth) — per sample, Spearman

| sample | none | logistic | adaptive-bare | auto-new | auto-old | +calib | +censor | +prune | +rank |
|---|---|---|---|---|---|---|---|---|---|
| nanosim-H9-cdna | 0.8739 | **0.9084** | 0.9077 | 0.9075 | 0.8099 | 0.9076 | 0.9076 | 0.8820 | 0.8315 |
| nanosim-H9-drna | 0.9113 | **0.9363** | 0.9350 | 0.9350 | 0.8736 | 0.9350 | 0.9351 | 0.8877 | 0.9350 |
| nanosim-NA12878-cdna | 0.8392 | **0.8867** | 0.8846 | 0.8824 | 0.7677 | 0.8842 | 0.8827 | 0.8564 | 0.7879 |
| nanosim-NA12878-drna | 0.8723 | **0.9192** | 0.9168 | 0.9170 | 0.8724 | 0.9171 | 0.9168 | 0.8805 | 0.9168 |
| tksm-RSII | 0.9289 | 0.9354 | 0.9380 | 0.9378 | 0.8959 | 0.9374 | **0.9384** | 0.8966 | 0.9380 |
| tksm-SQ2 | 0.9480 | 0.9459 | 0.9481 | 0.9480 | 0.9183 | 0.9479 | **0.9482** | 0.9186 | 0.9481 |

Panel B, per sample, MARD (lower is better):

| sample | none | logistic | adaptive-bare | auto-new | auto-old | +calib | +censor | +prune | +rank |
|---|---|---|---|---|---|---|---|---|---|
| nanosim-H9-cdna | 0.0797 | **0.0568** | 0.0575 | 0.0577 | 0.1403 | 0.0575 | 0.0577 | 0.0696 | 0.1241 |
| nanosim-H9-drna | 0.0418 | **0.0318** | 0.0324 | 0.0323 | 0.0529 | 0.0323 | 0.0323 | 0.0479 | 0.0324 |
| nanosim-NA12878-cdna | 0.0873 | **0.0582** | 0.0593 | 0.0600 | 0.1576 | 0.0594 | 0.0599 | 0.0713 | 0.1426 |
| nanosim-NA12878-drna | 0.0595 | **0.0370** | 0.0382 | 0.0381 | 0.0538 | 0.0381 | 0.0382 | 0.0508 | 0.0382 |
| tksm-RSII | 0.0613 | 0.0567 | 0.0545 | 0.0547 | 0.0785 | 0.0549 | **0.0542** | 0.0778 | 0.0545 |
| tksm-SQ2 | 0.0871 | 0.0872 | 0.0846 | 0.0848 | 0.1103 | 0.0849 | **0.0845** | 0.1101 | 0.0846 |

## Results: Panel A (LongBench, matched Illumina) — mean Spearman by family

| family | n | none | logistic | adaptive-bare | auto-new | auto-old | +calib | +censor | +prune | +rank |
|---|---|---|---|---|---|---|---|---|---|---|
| ont-cdna | 8 | 0.2291 | 0.2201 | 0.2207 | 0.2216 | **0.2348** | 0.2210 | 0.2214 | 0.2217 | 0.2313 |
| ont-drna | 8 | 0.2980 | 0.2921 | 0.2915 | 0.2927 | **0.3032** | 0.2917 | 0.2924 | 0.2923 | 0.3002 |
| pac-bio | 11 | 0.3404 | 0.3370 | 0.3384 | 0.3393 | **0.3444** | 0.3386 | 0.3391 | 0.3410 | 0.3408 |
| synthetic-drna | 1 | 0.8207 | 0.8254 | 0.8258 | 0.8244 | 0.8253 | 0.8257 | 0.8245 | 0.8183 | **0.8319** |

`auto-old` — the exact pre-demotion configuration — is the best arm on 25 of 28
LongBench samples, and `logistic` is the *lowest* arm on most of them. MARD on
this panel does not discriminate at all (every arm within 0.001 of 0.45).

## Results: Panel A (SIRV and spike-ins)

Mean by family; E0 rank metrics are undefined (equal-molar truth).

| family | n | metric | none | logistic | adaptive-bare | auto-new | auto-old | +calib | +censor | +prune | +rank |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sirv-E0 | 6 | MARD | 0.2919 | 0.2916 | 0.2800 | 0.2649 | **0.2573** | 0.2774 | 0.2669 | 0.2614 | 0.2801 |
| sirv-E2 | 4 | Spearman | 0.8157 | 0.8176 | 0.8162 | 0.8528 | **0.8844** | 0.8186 | 0.8491 | 0.8574 | 0.8162 |
| independent-sim | 2 | Spearman | 0.9634 | 0.9775 | 0.9834 | 0.9835 | **0.9852** | 0.9835 | 0.9835 | 0.9852 | 0.9835 |

These 12 samples side with Panel A, not Panel B: `auto-old` is best on all three
families.

**Caveat on `independent-sim`.** These two samples have exact counts truth, but
their truth is nearly degenerate: 100,000 reads spread over ~82,000 transcripts,
with a maximum of 6 reads on any transcript and 81% of transcripts at exactly 1
read. There are only six distinct truth values, so rank metrics there measure
detection rather than abundance estimation. Panel B, by contrast, averages
320-338 reads per expressed transcript across several orders of magnitude. On
the expressed-only universe used by the original harness, this run reproduces
that harness's published values (`none` 0.6827 vs 0.667849 reported;
`logistic` 0.7043 vs `auto` 0.710010 reported).

## Why the panels disagree

This was investigated rather than assumed. One hypothesis — that rank blending
compresses the dynamic range of the estimate — was tested and **rejected**: all
arms have essentially identical global dispersion on Panel B (CV ~25.0 against a
truth CV of 25.0).

What the data do show, on NanoSim NA12878 1D-cDNA:

**1. Global Spearman on Panel B is dominated by truth-zero transcripts.** 83.0%
of Panel B's universe (147,537 of 177,816) has exactly zero true reads. The
number of those an arm gives non-zero mass to orders the arms almost perfectly:

| arm | truth-zero transcripts given non-zero mass | mass on zeros | global Spearman |
|---|---|---|---|
| logistic | 2,931 | 17,871 | 0.8867 |
| adaptive-bare | 3,079 | 22,171 | 0.8845 |
| auto-new | 3,287 | 26,453 | 0.8824 |
| +prune | 5,343 | 48,986 | 0.8563 |
| none | 8,261 | 19,999 | 0.8391 |
| +rank | 18,143 | 24,378 | 0.7879 |
| auto-old | 21,077 | 79,392 | 0.7676 |

**2. Among genuinely expressed transcripts, `auto-old` is the best arm in every
abundance stratum.** Stratified Spearman over the 30,279 expressed transcripts:

| arm | low (<=10 reads) | mid (<=59) | high |
|---|---|---|---|
| none | 0.5419 | 0.6915 | 0.9279 |
| logistic | 0.5423 | 0.7066 | 0.9467 |
| auto-new | 0.5401 | 0.7098 | 0.9484 |
| **auto-old** | **0.5576** | **0.7148** | **0.9532** |

So the demoted features do not harm abundance accuracy among expressed
transcripts — they harm the model's ability to leave unexpressed transcripts at
zero, and that is what the global metric measures. Panel A's Illumina truth is
52.9% zero versus Panel B's 83.0%, and SIRV E2 has no zeros at all, which is
consistent with those panels favouring `auto-old`. This is a coherent
explanation supported on one sample; it is not established as the sole cause.

**Both readings are legitimate, but they answer different questions.** If the
metric of record is genome-wide rank fidelity over the full annotation — which
is what every table in this evaluation series reports — `logistic` wins and the
demoted features are harmful. If the question is abundance accuracy restricted
to expressed transcripts, `auto-old` is better. The published metric is the
former, so the recommendations below follow it.

## Follow-up: is rank blending salvageable by fixing the warm-up?

The plan lists "rank blending blends toward a non-converged EM" as a defect to
fix regardless of the benchmark. That was tested directly on NanoSim NA12878
1D-cDNA by sweeping `--coverage-warmup-iterations`.

The warm-up does converge, at **17,926 evaluations** — 179x the default cap of
100.

| `+rank` config | Spearman | MARD | truth-zero FPs | warm-up |
|---|---|---|---|---|
| `logistic` (baseline) | **0.8867** | **0.0582** | **2,931** | - |
| `adaptive-bare` | 0.8846 | 0.0593 | 3,079 | - |
| warmup = 100 (**default**) | 0.7879 | 0.1426 | 18,143 | not converged |
| warmup = 1,000 | 0.8378 | 0.0906 | 8,811 | not converged |
| warmup = 5,000 | 0.8386 | 0.0890 | 8,497 | not converged |
| warmup = 20,000 | 0.8388 | 0.0887 | 8,458 | **converged @ 17,926** |
| warmup = 100,000 | 0.8388 | 0.0887 | 8,458 | **converged @ 17,926** |

Two conclusions:

**The truncation defect is real and accounts for about half the damage.**
Spearman recovers 0.7879 -> 0.8388 and false positives halve, 18,143 -> 8,458.
The false-positive count and the Spearman deficit move together, which supports
the account in the previous section. Nearly all of the recovery arrives by 1,000
iterations; full convergence adds only a further +0.0010.

**Fixing it does not rehabilitate the feature.** Fully converged, rank blending
is still 0.0479 below `logistic` and 0.0458 below `adaptive-bare` — about five
times the rule 1 tolerance — and still carries 2.9x logistic's false positives.
It also costs 224.9 s of warm-up on a run that otherwise takes ~146 s.

The residual deficit is therefore intrinsic to the design rather than to the
truncation: blending 20% of a coverage-*free* estimate into a
coverage-*corrected* one reintroduces the bias the coverage model exists to
remove. Iterating the warm-up longer cannot fix that, and the measurements show
it does not.

**Caveat.** The main EM does not converge either — `em_evaluations: 1000,
em_converged: false` at threshold 0.001 for *every* arm on this data, including
`logistic`. The problem is therefore not "converged versus truncated" in
absolute terms; it is that the warm-up was given a 10x smaller iteration budget
than the estimate it is blended against. The general non-convergence of the main
EM at the shipped cap is a separate issue worth its own investigation.

## Follow-up: can the adaptive kernel be salvaged?

The kernel has a genuine advantage that the headline metric hides, so three
salvage routes were tested. None of them yields a default.

**Finding 0: `auto` is not one kernel, it is three.** The technology dispatch
selects a different model per `--seq-tech`:

| technology | kernel |
|---|---|
| ont-cdna | `cross_fitted_adaptive_endpoint` |
| ont-drna | `ont_drna_competing_risks` |
| pac-bio, pac-bio-hifi | `pacbio_guarded_physical_endpoint` |

Every result reported as "`auto` versus X" is therefore an average over three
distinct algorithms, and they do not behave alike. This alone is a reason to
stop reporting `auto` as a single arm.

**Route 1: the expressed-transcript gain is real.** Over all 6 Panel B samples,
restricted to transcripts the truth actually expresses:

| arm | global | expressed-only | low | mid | high | FP zeros |
|---|---|---|---|---|---|---|
| `logistic` | **0.9220** | 0.9452 | 0.6187 | 0.7948 | 0.9609 | **1,812** |
| `adaptive-bare` | 0.9217 | **0.9488** | **0.6211** | **0.8005** | **0.9660** | 1,931 |
| `auto-new` | 0.9213 | 0.9501 | 0.6229 | 0.8023 | 0.9672 | 1,998 |

The bare kernel beats `logistic` on expressed transcripts in all three
abundance strata and on 5 of 6 samples (+0.0036 mean; TKSM RSII +0.0104,
SQ2 +0.0085). Its false-positive penalty is only +119 transcripts (+6.6%), and
that is enough to cancel the gain.

**Route 2: equalise the zero handling — fails.** Applying the same post-hoc
threshold to every arm (zeroing estimates below t) and recomputing global
Spearman, mean over the 6 Panel B samples:

| arm | t=0 | t=0.25 | t=0.5 | t=1 | t=2 |
|---|---|---|---|---|---|
| `logistic` | 0.9220 | 0.9244 | 0.9246 | 0.9245 | 0.9078 |
| `adaptive-bare` | 0.9217 | 0.9244 | 0.9247 | 0.9245 | 0.9081 |
| delta | -0.0003 | +0.0000 | +0.0001 | +0.0000 | +0.0002 |

The gap stays at +0.0001 at every threshold. The expressed-transcript advantage
does not transfer to the metric of record: with 83% of the universe at truth
zero, a +0.0036 gain over 17% of the transcripts is worth ~+0.0001 globally.

Incidentally, a threshold of 0.5 improves *every* arm by ~+0.0026 (`logistic`
0.9220 -> 0.9246) — more than the entire adaptive stack contributes. Cheap,
orthogonal to coverage modelling, and worth investigating separately.

**Route 3: a PacBio-conditional default — fails on the right baseline.** The
kernel beats `logistic` on **13 of 13** PacBio samples across both panels
(+0.0016 mean), which is the only cross-panel-consistent signal in the whole
evaluation. But measured against `none`:

| arm | mean d vs `none` (PacBio, n=13) | wins |
|---|---|---|
| `+rank` | +0.0011 | 6/13 |
| `auto-new` | -0.0002 | 4/13 |
| `adaptive-bare` | -0.0009 | 5/13 |
| `logistic` | **-0.0025** | **4/13** |

The 13/13 sweep is "`logistic` underperforms on PacBio", not "the kernel helps
on PacBio". No coverage model reliably beats doing nothing there. **That
`logistic` loses to `none` on PacBio (-0.0025, 4/13) is an independent finding
worth its own investigation.**

The exception is Panel B PacBio under exact truth, where `adaptive-bare` beats
`none` by +0.0091 on RSII and ties on SQ2 (+0.0001), and beats `logistic` on
both. That is the one place the kernel clears both baselines — but n=2 is far
too thin to set a default, and Panel A's 11 PacBio samples point the other way.

**Conclusion.** The adaptive kernel is not salvageable as a default on the
current evidence. Its one genuine advantage is confined to expressed transcripts
and is worth ~+0.0001 under the metric of record. The PacBio lead is real enough
to be worth testing properly, and the required experiment is well defined:
more PacBio samples with exact read-level truth, scored against `none` as well
as `logistic`.

## Where the remaining headroom is (Panel B, exact truth)

All numbers below are measured on the saved Panel B quantifications, not
projected. Several promising ideas were tested and killed; they are recorded so
they are not re-proposed.

**The problem is ambiguity, not coverage.** 89.5% (ONT cDNA) to 95.1% (PacBio
HiFi) of assigned read mass is ambiguous. Unique reads alone give Spearman
0.5559 / 0.4733; the EM lifts that to 0.8867 / 0.9459. Coverage modelling is a
second-order correction (`logistic` buys +0.0264 over `none`) on top of a
first-order ambiguity-resolution problem.

**Error decomposition for `logistic`** (mean over 6 Panel B samples):

| | Spearman | headroom |
|---|---|---|
| `logistic` as shipped | 0.9220 | - |
| + perfect **detection** (truth-zeros forced to 0) | 0.9481 | **+0.0261** |
| + perfect **abundance** (expressed set to truth) | 0.9765 | +0.0545 |

1,812 false-positive transcripts carrying only **0.156%** of the read mass cost
**0.026 Spearman** — as much as the entire benefit of coverage modelling. This
is the largest identified single lever, and it is orthogonal to the coverage
model.

**Ideas tested and rejected:**

| idea | measured result |
|---|---|
| length / effective-length correction | **oracle** gain -0.0001 Spearman, -0.0019 MARD |
| ambiguity-gated coverage correction | 0.0000 at every gate threshold |
| naive detection thresholds | +0.0026 of the +0.0261 available (10%) |
| aggressive detection (no unique read & est<5, <20) | **-0.0361, -0.1138** |

The length-dependent residual is real and large — mean log2(est/truth) runs
-0.21 (short) to -0.68 (long) on ONT cDNA, and `logistic` *worsens* the long end
(-0.68 -> -0.75) — but it is nearly rank-preserving, so correcting it perfectly
buys nothing under rank metrics. It would matter for a metric of record that
weights absolute abundance.

**Why naive detection fails, and what would work.** The two classes are
separable in bulk but overlap exactly where the false positives live:

| group | median est. reads | p90 | fraction with no unique read |
|---|---|---|---|
| truly expressed (26,196) | 30.68 | 429.45 | 0.617 |
| false positive (2,931) | 1.32 | 15.38 | **0.987** |

Of transcripts with est <= 1 and no unique read, 65.0% are truly zero — but 35%
are real. At est <= 5 the split is 40.8% / 59.2%. So a hard cutoff necessarily
trades true positives for false positives, which is why the aggressive rules
lose 0.04-0.11 Spearman. The discriminating information exists (a 23x median
separation, and 98.7% versus 61.7% on unique-read presence); it needs a
probabilistic presence/absence component, not a threshold.

**Ranked proposals, by measured headroom:**

1. **A presence/absence (zero-inflated or spike-and-slab) component in the EM.**
   Up to +0.0261 Spearman, roughly ten times what the entire adaptive stack
   delivers, and independent of which coverage kernel is used. Should emit a
   per-transcript probability of presence rather than a hard call.
2. **An empirical endpoint mixture in place of the fixed inverse-logistic
   shape.** The physical-endpoint diagnostics show large, learnable, and highly
   variable structure across libraries — TKSM SQ2 fits 2.7% intact / 54.2%
   3'-truncated, while LongBench H69 fits 26.6% intact / 26.0% 3'-truncated. A
   single fixed curve cannot express that range. This is what the PacBio kernel
   already does, and it is the one kernel that beats `logistic` consistently.
3. **A cheap immediate win: a 0.5-read floor.** Applied post-EM it gains +0.0026
   on every arm — more than the whole adaptive stack — at no cost. It captures
   only 10% of the detection headroom, so it is a stopgap, not the fix.

## Cost (Panel B, serial, idle node, 16 threads)

| arm | mean wall (s) | mean peak RSS (GB) |
|---|---|---|
| none | 112.0 | 3.21 |
| logistic | 113.9 | 3.33 |
| adaptive-bare | 139.7 | 3.42 |
| auto-new | 140.1 | 3.42 |
| auto-old | 140.7 | 3.42 |
| +prune | 146.7 | 3.42 |
| +rank | 146.1 | 3.42 |

**The reported "1.8x faster" belongs to the branch, not to `auto`.** On the same
branch, `auto-new` is **23% slower** than `logistic` (140.1 s vs 113.9 s) and
uses slightly more memory. The speedup in the review doc came from comparing dev
`auto` against *main* `logistic`; main `logistic` averages 269.9 s on Panel B,
so demote `logistic` is 2.4x faster than main and demote `auto` is 1.9x faster.
Choosing `auto` over `logistic` on the current branch costs time, it does not
save it.

## Cross-branch control

`--model-coverage` reproduces across branches, but not exactly:

| panel | max abs delta vs `main` | direction |
|---|---|---|
| B | 0.0053 | demote higher on all 6 |
| A | 0.0021 | main higher on 24/28 |

**This sets the noise floor for the whole experiment.** Any effect smaller than
~0.005 Spearman cannot be distinguished from branch-level drift in the logistic
path itself. Alignment calibration (+0.0009) and censoring (+0.0042) are both
inside that floor.

## Decision rules, applied as written

### Rule 1 — a feature is an `auto` default only if it beats `logistic`, with no dataset family regressing more than 0.01 Spearman

| feature | arm | union mean d | worst family d | verdict |
|---|---|---|---|---|
| rank blending | `+rank` | +0.0009 | **-0.0878** (nanosim-cdna) | **FAILS** |
| dominance pruning | `+prune` | +0.0004 | **-0.0437** (nanosim-drna) | **FAILS** |
| alignment calibration | `+calib` | +0.0009 | -0.0017 (nanosim-drna) | passes |
| censoring | `+censor` | +0.0042 | -0.0024 (nanosim-cdna) | passes |

### Rule 2 — rank blending and dominance pruning stay demoted unless they clear rule 1

Neither clears it. **They stay demoted.** The activation-gate re-derivation
contemplated by rule 2 is therefore not required.

### Rule 3 — alignment calibration and censoring are on probation

Both are positive on Panel A LongBench (`+calib` +0.0008, `+censor` +0.0012), so
they survive probation on the letter of the rule. But both are inside the
0.005 noise floor established by the cross-branch control, on both panels. The
honest statement is that this experiment **cannot resolve** whether either helps.

### Rule 4 — if `adaptive-bare` does not beat `logistic` on the union, `auto` should not be a default at all

| | `adaptive-bare` vs `logistic` | `auto-new` vs `logistic` |
|---|---|---|
| union (46 samples) | +0.0005 | +0.0048 |
| Panel B only | **-0.0003** | **-0.0007** |
| Panel A LongBench | +0.0006 | +0.0015 |

The union is positive only because Panel A contributes 40 of 46 samples. On the
panel with exact truth, `adaptive-bare` does not beat `logistic` — it ties, well
inside the noise floor. Under rule 5 the Panel B reading governs. **Rule 4 is
satisfied: `auto` should not be the default.**

### Rule 5 — if the panels disagree, prefer Panel B and report the disagreement

They disagree on every arm, in sign, systematically:

| arm | Panel A (all 40) | Panel B | agree? |
|---|---|---|---|
| none | +0.0055 | -0.0264 | no |
| adaptive-bare | +0.0006 | -0.0003 | no |
| auto-new | +0.0015 | -0.0007 | no |
| auto-old | **+0.0103** | **-0.0657** | no |
| +calib | +0.0008 | -0.0005 | no |
| +censor | +0.0012 | -0.0006 | no |
| +prune | +0.0019 | -0.0350 | no |
| +rank | **+0.0072** | **-0.0291** | no |

Panel A ranks `none` above `logistic` (+0.0055, 24 wins / 4 losses). A panel on
which "no coverage model" outperforms the shipping coverage model has little
power to select coverage-model features, which is a plausible account of how
these features were selected in the first place.

## Recommendations

**1. Rank blending — keep opt-in, and consider removing it outright.**
Rule 1: fails, worst family -0.0878 on NanoSim 1D-cDNA, nearly nine times the
0.01 tolerance. MARD on that family degrades from 0.0582 to 0.1426, 2.4x worse.
Confirmed on the only panel with exact truth at realistic depth. The committed
demotion stands. The follow-up above additionally shows that the feature is not
salvageable by fixing the non-converged warm-up: with a fully converged warm-up
it still fails rule 1 by ~5x while costing 224.9 s of extra EM. Since its
best-achievable configuration still loses decisively to both `logistic` and
`adaptive-bare`, "remove" is defensible and "opt-in" is the conservative floor.

**2. Dominance pruning — keep opt-in (do not revert the demotion).**
Rule 1: fails, worst family -0.0437, more than four times the tolerance, plus
the largest wall-time cost of any single feature (+5.0% over `adaptive-bare`).
The committed demotion stands.

**3. Alignment calibration — make it opt-in.**
Rule 1 passes on the letter (+0.0009 union, worst family -0.0017), but rule 3
put it on probation and the measured effect is five times smaller than the
0.0053 cross-branch reproduction spread. There is no evidence it improves on
`logistic`; there is also no evidence it harms. A default that cannot be shown
to do anything should not be a default, and it is currently `Auto` in
`prog_opts.rs`. Demote to `none` and keep it reachable.

**4. Censoring — make it opt-in.**
Same reasoning. Its +0.0042 union mean is the largest of the two survivors, but
it is driven by SIRV E2 (+0.0315 across 4 samples of 69 transcripts with 4
concentration levels) and is negative on Panel B (-0.0006). Inside the noise
floor. It is currently `Auto`; demote to `none`.

**5. Beyond the four features: `auto` should not be the default coverage model.**
Rule 4 is satisfied. `--model-coverage` (logistic) wins outright on 4 of 6
Panel B samples, ties on the other 2, has the best MARD on 4 of 6, and is 23%
faster than `auto-new` on the same branch. Ship `--model-coverage` as the
default and keep `--coverage-model auto` opt-in until it can be shown to beat
logistic on exact truth.

## Does anything contradict the committed demotion?

**No.** Panel B confirms the demotion strongly and reproduces the review doc's
numbers exactly. Nothing here recommends reverting it.

Two findings do complicate the surrounding narrative, and are stated plainly:

- **The demotion is correct but not sufficient.** The same evidence that
  demotes rank blending and dominance pruning also shows that the remaining
  `auto` stack does not beat `logistic` on exact truth. Stopping at the
  demotion leaves a default that is slower and no more accurate.
- **40 of the 46 samples in this benchmark favour the pre-demotion
  configuration.** That is not a reason to revert — those panels are proxy
  truth, equal-molar controls with undefined rank metrics, or near-degenerate
  truth with six distinct values — but it should be recorded that the
  conclusion rests on 6 samples out of 46, chosen on the pre-registered grounds
  of truth quality rather than on count.

## Deliverables and provenance

- Per-sample machine-readable results: `coverage-rebenchmark-2026-08-03.tsv`
  (448 rows; arm x sample x metric, plus wall time, peak RSS, evaluation
  universe size, reference md5, and the activated-feature diagnostics).
- Harness: `scripts/rebench-2026-08-03/` in this repo (`run_arms.py`,
  `metrics.py`, `analyze.py`, `make_deliverable.py`, `panelB_manifest.tsv`,
  `panelA_sirv_manifest.tsv`). Copied out of node-local scratch deliberately —
  the plan notes that `/scratch1` is neither shared nor backed up.
- Raw run outputs (`.quant`, `.meta_info.json`, `/usr/bin/time -v` logs) remain
  under `/scratch1/rob/long-read-ecosystem/rebench-2026-08-03/`, which is not
  durable storage.
- Binaries: `fix/demote-rank-blend-and-dominance-pruning` at `af8f3eb` and
  `origin/main`, both oarfish 0.10.3, built into separate target dirs.
- Panel B was copied from `/fs/cbcb-lab/rob/students/zahra/` — that NFS export
  is mounted directly on the benchmark host, so the `nexuscbcb01` hop in the
  plan was unnecessary (and that host rejects the available key).

### Caveats

- Panel A's absolute Spearman (0.21-0.43) is low because matched-Illumina truth
  is a cross-protocol proxy. Arm-to-arm comparison within a sample is still
  valid — every arm is scored over an identical, hash-verified reference — but
  Panel A and Panel B absolute values are not comparable to each other.
- All runs used oarfish defaults outside the arm definitions. The default EM
  accelerator is `None` and the EM reports `converged=false` at the evaluation
  cap on the LongBench samples, for every arm equally.
- Rank-metric conclusions for SIRV E0 are impossible by construction; those six
  samples contribute only MARD and dispersion.
- The mechanism proposed for the panel disagreement (zero-transcript
  discrimination) was verified on one Panel B sample, not all six.
