# Endpoint-grid geometry (2026-07-24)

## Hypothesis

The joint endpoint model in `src/util/endpoint_probability.rs` bins strand-aware *fractional*
gaps onto a uniform 20x20 grid. Because
`five_gap + three_gap = 1 - aligned_fraction <= 1` and `floor(a) + floor(b) <= floor(a + b)`, a
valid alignment can only ever land in the triangle `x + y <= GRID - 1` -- 210 of 400 cells. The
remaining 190 are structurally unreachable, yet they carry symmetric Dirichlet prior mass and are
averaged into their neighbours by the 3x3 smoother in `smooth_counts`.

The smoother is the substantive defect: it divides by the full neighbour weight while
structurally-empty neighbours contribute a zero count, so anti-diagonal cells are deflated by a
*cell-specific* factor that no scalar prior can undo.

## Where the grid is live

`src/bulk.rs` routes PacBio/PacBio-HiFi `auto` + `Full` to `PacbioPhysicalEndpoint`, a
nucleotide-scale kernel. The fractional grid is therefore reached only by **ONT cDNA `auto`**,
**ONT dRNA `auto`** (where it feeds `degradation_probability.rs`), and the explicit
`--coverage-model endpoint|hybrid|adaptive` modes. This is an ONT-facing change and cannot move
the Kinnex Figure 2 panels.

A correction to the reading of
[coverage-signal-review-2026-07-21.md](coverage-signal-review-2026-07-21.md): the
`coverage_rho_length = -0.615` and `mean_coverage_log_range = 0.092` figures reported there come
from `kinnex-wtc11/coverage-signal-export-20260721`, whose `.meta_info.json` records
`technology_kernel = pacbio_adaptive_logistic` and `effective_ablation = NoEndpoint`. Under
`NoEndpoint` the endpoint weight is zero (`hybrid_probability.rs`), so those exported
`coverage_probs` are **logistic-only**. Likewise `weighted-eq-collapse-h69-50k` ran
`pacbio_guarded_physical_endpoint`. Neither figure describes the endpoint grid, and no change to
the grid can move them.

## Data

Sixteen grid-live ONT exports were generated for this round (none existed previously):

```
for s in H69 H146 H211 H526 H1975 H2228 HCC827 SHP77; do for t in cdna drna; do
  tech=$([ $t = cdna ] && echo ont-cdna || echo ont-drna)
  oarfish --alignments longbench/bam/$s-$t-first50k.name.bam --seq-tech $tech \
    --coverage-model auto --coverage-ablation full --filter-group no-filters --threads 4 \
    --write-coverage-signals --coverage-signal-sample-rate 2 \
    --output oarfish-evaluation-data/ont-signal-export-20260724/$s-$t
done; done
```

Kernels confirmed from `.meta_info.json`: `cross_fitted_adaptive_endpoint` (8 cDNA),
`ont_drna_competing_risks` (8 dRNA).

## Findings

### 1. The reachable region is exactly the predicted triangle

Computed from the raw `starts`/`ends`/`lengths`/`strands` columns, so independent of which kernel
produced the export: exactly **210 distinct cells occupied**, **zero** records in the other 190,
across every dataset checked. 47.5% of the grid carries prior mass that can never be matched by
data, and the smoother borrows from it.

Records within one smoothing step of the anti-diagonal boundary: PacBio Kinnex 0.00%, synthetic
ONT dRNA 8.86%.

### 2. The prior grid is pinned at its lower boundary in every ONT library

`selected_prior_mass` is chosen from `PRIOR_GRID = [10, 30, 100, 300, 1000, 3000, 10000]` by
cross-fitted held-out predictive likelihood. Result: **10.0 in 16 of 16 libraries**, for both
technologies -- the smallest value offered.

This is a boundary solution: the model's own held-out criterion is asking for less prior than the
grid can express, and the search is truncated. Two consequences:

- The absorption argument for the feasible-prior variant does not apply. Concentrating the prior
  on 210 cells *raises* the per-cell share by `400/210 = 1.905x`, and the selector cannot
  compensate by stepping down because it is already at the floor. `endpoint-feasible-prior` is
  therefore predicted to be neutral-to-harmful, not neutral-to-helpful.
- Extending `PRIOR_GRID` downward is a separate, independently motivated candidate.

`capped_reads` is also high (6,992-18,626 of ~50k reads), i.e. the coverage term routinely wants
to express more evidence than `--coverage-max-bayes-factor 4` permits.

### 3. Endpoint structure is at nucleotide scale, which the grid cannot represent

Per-end gap distributions on real LongBench libraries:

| Library | 5' == 0 | 5' <= 10 nt | 5' median | 3' == 0 | 3' <= 10 nt | 3' median |
|---|---:|---:|---:|---:|---:|---:|
| ONT cDNA H69 | 34.1% | 47.7% | 17 nt | 33.5% | 47.4% | 18 nt |
| ONT cDNA H146 | 36.0% | 50.9% | 9 nt | 32.8% | 47.3% | 18 nt |
| ONT dRNA H69 | 10.6% | 20.6% | 118 nt | 41.4% | 59.6% | 3 nt |
| ONT dRNA H146 | 10.9% | 21.2% | 111 nt | 43.6% | 60.6% | 3 nt |

For a median transcript the grid's finest bin is ~90 nt wide. About half of cDNA records have
*both* endpoints inside 10 nt, and dRNA's 3' end has a median gap of **3 nt** -- the anchored end
that distinguishes dRNA -- while its 5' end is diffuse at ~115 nt. The grid resolves none of this.

The corresponding PacBio numbers (5' median 1 nt, 3' median 2 nt) explain why the nucleotide-scale
physical kernel won there, and why its `tau5`/`tau3` pin at the `MIN_TAU = 3.0` clamp.

Note that the small synthetic dRNA sample in `synthetic-drna-small/` is **not** representative on
this axis (5' median 961 nt, 3' median 380 nt); design decisions should use the real libraries.

A single monotone log warp `b(g;L) = floor(B*ln(1+g/g0)/ln(1+L/g0))` is fine near zero and
asymptotically fractional in the tail, so it would give cDNA resolution at both ends and dRNA
resolution on its anchored 3' end while leaving its diffuse 5' end in coarse bins -- without any
technology conditioning.

## Implementation

Three opt-in `CoverageAblation` variants, default behavior unchanged:

- `endpoint-feasible-smoothing` -- the smoother sums kernel weights over reachable neighbours only,
  leaves unreachable cells at zero, and rescales to preserve total mass so the correction changes
  only the *shape* of the table.
- `endpoint-feasible-prior` -- Dirichlet prior spread over the 210 reachable cells only.
- `endpoint-feasible-geometry` -- both.

A `GridConfig` is derived once in `adaptive_endpoint_probabilities` and threaded into
`smooth_counts` and a new `cell_probability` helper that replaces the four inlined copies of the
estimator (scoring, prior selection, the per-candidate loop, and the fold-variance loop). With a
default `GridConfig` the helper is exactly the historical expression.

`endpoint_probabilities` (used by `--coverage-model endpoint|hybrid`) deliberately keeps legacy
behavior, since `docs/evaluation.md` requires `endpoint` as a fixed comparison baseline.

`src/bulk.rs` keeps PacBio `auto` on the promoted physical kernel for these variants, so the
PacBio panels act as invariance controls. Without that guard every new ablation would silently
drop PacBio onto the grid and re-litigate a comparison the grid already lost.

## Verification gates

- **Refactor is a no-op**: `full` output is byte-identical to the pre-change binary on H69-cdna,
  H69-drna, H146-drna, H526-cdna.
- **PacBio invariance**: `endpoint-feasible-smoothing` and `endpoint-feasible-geometry` produce
  byte-identical `.quant` to `full` on H69-pb.
- Unit tests: the reachable region is exactly 210 cells; valid alignments never reach an
  infeasible cell; the smoother ignores unreachable neighbours; total mass is preserved; a default
  `GridConfig` reproduces the legacy estimator.

  Note the interior-cell test asserts exact agreement only for an isolated interior blob, where
  the mass-preserving rescale is the identity. On real data the rescale is a global factor
  slightly below one (boundary cells gain mass, so the pre-rescale total rises), so interior cells
  shift by that common factor. The correction is a change of *shape*, not a claim that interior
  cells are bit-identical.

### 4. Signal-level baselines for the grid-live path

`scripts/analyze_alignment_coverage_signals.py` had never been run on ONT exports and crashed on
them (`ZeroDivisionError`): candidate pruning can zero every score probability for a read. Fixed
by falling back to a uniform distribution when a normalizer is zero.

Baselines over the 8+8 ONT libraries, `full` arm:

| Metric | ONT cDNA | ONT dRNA |
|---|---:|---:|
| `mean_coverage_log_range` | 0.9871 | 0.8037 |
| `coverage_rho_length` | -0.3199 | -0.2381 |
| `coverage_winner_agreement` | 0.3076 | 0.2909 |
| `coverage_high_conf_agreement` | 0.3387 | 0.3178 |
| `combined_high_conf_agreement` | 1.0000 | 1.0000 |

So on its live path the endpoint grid carries considerably more within-class contrast than the
0.092 previously attributed to it, but it still never overturns a high-confidence alignment-score
winner.

## Results

28-case panel, `--threads 4`, arms `none full endpoint-feasible-smoothing endpoint-feasible-prior
endpoint-feasible-geometry`. Artifacts: `oarfish-evaluation-data/endpoint-feasible-20260724/panel`.

Means by technology:

| Technology | Arm | CCC | Spearman | MARD | RMSE |
|---|---|---:|---:|---:|---:|
| ONT cDNA (8) | none | 0.275396 | 0.229145 | 0.456792 | 7872.26 |
| | full | 0.282916 | 0.231478 | 0.457437 | 7799.43 |
| | feasible-smoothing | 0.282878 | 0.231554 | 0.457431 | 7800.37 |
| | feasible-geometry | 0.282897 | 0.231544 | 0.457435 | 7800.06 |
| ONT dRNA (9) | full | 0.552123 | 0.360554 | 0.434563 | 5578.33 |
| | feasible-smoothing | 0.552110 | 0.360510 | 0.434635 | 5578.54 |
| PacBio (11) | full | 0.508652 | 0.340968 | 0.443861 | 3893.18 |
| | feasible-smoothing | 0.508652 | 0.340968 | 0.443861 | 3893.18 |

Deltas versus `full`:

| Technology | Arm | dCCC | dSpearman | dMARD | dRMSE |
|---|---|---:|---:|---:|---:|
| ONT cDNA | feasible-smoothing | -0.000038 | **+0.000076 (8/8)** | -0.000006 | +0.94 |
| ONT cDNA | feasible-prior | -0.000029 | -0.000004 (3/8) | +0.000002 | +0.47 |
| ONT cDNA | feasible-geometry | -0.000018 | +0.000066 (8/8) | -0.000002 | +0.63 |
| ONT dRNA | feasible-smoothing | -0.000013 | -0.000044 | +0.000072 | +0.20 |
| ONT dRNA | feasible-prior | +0.000002 | -0.000002 | -0.000000 | -0.05 |
| PacBio | all three | **0.000000 (0/11)** | 0.000000 | 0.000000 | 0.00 |

PacBio is bit-identical on all 11 cases, confirming the `bulk.rs` guard. `coverage_seconds` is
unchanged (cDNA 0.161 vs 0.163 s; dRNA 0.540 vs 0.543 s), and peak RSS is unchanged.

**This is a measured null**, and both halves have a complete mechanistic explanation.

### Why the smoothing correction cannot matter

The fix only moves cells at or near the anti-diagonal. Measured exposure on real libraries:

| Library | on anti-diagonal | within one smoothing step |
|---|---:|---:|
| ONT cDNA H69 | 1.71% | 6.04% |
| ONT cDNA H146 | 1.53% | 5.33% |
| ONT dRNA H69 | 1.06% | 4.05% |
| ONT dRNA H146 | 0.97% | 3.75% |
| synthetic dRNA (small) | 2.43% | 8.86% |

Only ~1-2% of records sit where the bias is strongest, and the resulting per-cell correction is
then diluted by the support gate, the log-linear combination with the logistic term, and the
Bayes-factor cap. A 5th-decimal effect is the correct expectation, not a disappointment. Note the
small synthetic sample overstates this exposure by roughly 2x and should not be used to size the
effect.

### Why the prior correction cannot matter

`selected_prior_mass` is 10.0 in 16/16 ONT libraries. Spread over 400 cells that is **0.025
pseudo-counts per cell** -- the prior is effectively switched off, so redistributing it over 210
cells instead of 400 changes nothing measurable. The selection is also identical across all four
arms (10.0 in 17/28 panel cases in every arm): the 1.905x change in per-cell share is smaller than
one step of the ~3x-spaced `PRIOR_GRID`, so the selector cannot even register it.

## Decision

**Do not promote.** Neither variant meets the promotion bar in `evaluation.md`: median primary
metrics do not improve, and no technology class shows a consistent primary-metric gain. The only
consistent signal is ONT cDNA Spearman +0.000076 (8/8 libraries), which is a secondary metric and
is an order of magnitude below the +0.000252 delivered by the last promoted candidate.

Both changes are nonetheless *correct* -- the smoother should not borrow from structurally-empty
cells, and the prior should not be spread over unreachable ones -- so they are retained behind
their ablation flags as opt-in, with the default path byte-identical.

## What this rules in

The two cheap geometry corrections are exhausted. The measurements above point the remaining
leverage squarely at **resolution**, not feasibility:

- ~48% of ONT cDNA records have *both* endpoint gaps within 10 nt; the grid's finest bin is ~90 nt.
- ONT dRNA's 3' gap has a median of 3 nt with 41% exactly zero.
- This affects 100% of records, versus the ~1.7% the feasibility fix could reach.

A monotone log warp `b(g;L) = floor(B*ln(1+g/g0)/ln(1+L/g0))` with `g0 = 2` nt is fine near zero
and asymptotically fractional in the tail, so one symmetric axis serves cDNA (both ends sharp) and
dRNA (3' sharp, 5' diffuse) without technology conditioning. Under such a warp the feasibility
mask, the Jacobian, and the diffuse prior all collapse into a single object -- the integer
lattice-point count of a cell -- which is also why the feasibility work above is a prerequisite
rather than wasted effort.

The separate finding that `PRIOR_GRID` is pinned at its lower bound in 16/16 ONT libraries is an
independent, one-line candidate: the held-out selector is asking for less regularization than the
grid can express.

---

# Round 2: prior-grid extension and a resolution/measure decomposition

## Extended prior grid (`endpoint-wide-prior-grid`): rejected

`PRIOR_GRID` was extended downward to `[1, 3, 10, ...]`. The selector immediately moves to the new
floor -- 1.0 on ONT cDNA, 3.0 on ONT dRNA -- confirming the boundary was binding. It buys nothing:

| Technology | dCCC | dSpearman | dMARD | dRMSE |
|---|---:|---:|---:|---:|
| ONT cDNA (8) | +0.000006 | -0.000003 | -0.000002 | -0.05 |
| ONT dRNA (9) | -0.000005 | +0.000001 | -0.000000 | +0.05 |
| PacBio (11) | 0.000000 | 0.000000 | 0.000000 | 0.00 |

This is a 6th-decimal null, and the reason is arithmetic: at `prior = 10` spread over 400 cells the
per-cell prior is already 0.025 pseudo-counts. Lowering it to 1.0 makes it 0.0025. Both are
indistinguishable from zero against real cell counts, so the held-out selector's preference for a
smaller prior has no consequence for quantification. **The prior is not a lever on this model.**
Artifacts: `oarfish-evaluation-data/wide-prior-grid-20260724`.

## Is a finer grid worth building? A truth-anchored answer

Before implementing a warped or finer axis, the question was settled offline against ground truth.
The large synthetic dRNA BAM (`eval/parity/txp_cpp.rh.bam`) encodes the source transcript in each
read name (`NM_005885_3217_aligned_...`), so exporting its signals in alignment mode -- which,
unlike raw-read mode, preserves read names -- yields per-read truth.

`scripts`-external analysis trained cross-fitted endpoint models at several resolutions on unique
reads and scored **118,549 ambiguous reads whose true transcript is among the candidates**. Each
scheme selects its own Dirichlet prior, so finer grids are not penalised for needing more
smoothing. Comparisons use probability *per integer gap pair*, since raw cell probabilities are not
comparable across partitions.

| Scheme | cells | picks_true | delta | picks_true + measure |
|---|---:|---:|---:|---:|
| uniform-20 (current) | 400 | 0.3627 | -- | **0.3785** |
| uniform-30 | 900 | 0.3680 | +0.0053 | 0.3782 |
| uniform-40 | 1600 | 0.3650 | +0.0023 | 0.3799 |
| uniform-64 | 4096 | 0.3662 | +0.0034 | 0.3794 |
| warp-20 (g0=2) | 400 | 0.3515 | -0.0113 | 0.3791 |
| warp-40 (g0=2) | 1600 | 0.3601 | -0.0026 | **0.3829** |

Two conclusions, and they point away from the larger change:

1. **Resolution is nearly irrelevant.** Finer uniform grids gain at most +0.005 and
   non-monotonically (30 > 64 > 40 > 20), which is noise. The warp is *worse* than the current grid
   without a measure correction, because it concentrates bins where the current model already has
   adequate mass while coarsening the tail.
2. **The measure correction is the real effect, and it is resolution-free.** Dividing by the number
   of integer gap pairs a cell covers gains **+0.0158 at the current resolution** -- three times
   the best resolution gain -- and once applied, every scheme lands within 0.005 of every other
   (0.3782-0.3829). The oracle log Bayes factor moves the same way: 0.0625 -> 0.0973 for
   uniform-20, and 0.0267 -> 0.1046 for warp-20.

So a new warped-axis module would buy roughly +0.004 discrimination over the measure correction
alone, for a new length-class table, an integer lattice-area routine and a density-domain smoother.
**Not worth building.** The measure correction is a few lines on the existing grid.

This reverses the earlier decision to demote the Jacobian on a magnitude heuristic (that its
implied term is larger than the existing coverage signal, and that favouring shorter candidates
resembles the `MALAT1-256` containment failure). The heuristic reasoned about magnitude and
direction; the truth-anchored measurement shows the correction moves discrimination the *right*
way. The containment risk remains real and is what the panel must now test -- better within-class
discrimination is necessary but not sufficient, and the rejected unique-read-profile candidate is
the standing example of ranking improving while calibration collapses.

## Implementation

`endpoint-nt-measure` divides each candidate's cell probability by
`cell_lattice_area(cell, len)` before the per-read normalization, where a bin's lattice width is
`ceil((b+1)*len/GRID) - ceil(b*len/GRID)` -- exact for the binning `EndpointModel::cell` performs.
`endpoint-nt-measure-geometry` combines it with both feasible-region corrections.

Gates re-verified: default `full` remains byte-identical, and PacBio is invariant under both new
ablations.

## Results

28-case panel, arms `full endpoint-nt-measure endpoint-nt-measure-geometry`. Artifacts:
`oarfish-evaluation-data/nt-measure-20260724`.

| Technology | Arm | dCCC | dSpearman | dMARD | dRMSE |
|---|---|---:|---:|---:|---:|
| ONT dRNA (9) | nt-measure | +0.000004 (4/9) | **+0.000515 (8/9)** | **-0.000362 (8/9)** | +0.64 |
| ONT dRNA (9) | nt-measure-geometry | -0.000004 | +0.000472 (8/9) | -0.000294 (8/9) | +0.82 |
| ONT cDNA (8) | nt-measure | -0.000085 (1/8) | -0.000004 (2/8) | +0.000025 (2/8) | +2.67 |
| ONT cDNA (8) | nt-measure-geometry | -0.000097 (0/8) | +0.000049 (4/8) | +0.000025 | +2.96 |
| PacBio (11) | both | 0.000000 | 0.000000 | 0.000000 | 0.00 |

The truth-bearing `synthetic-drna` case is the largest single win: dSpearman **+0.003034**,
dCCC +0.000351, dRMSE -2.4.

This is the signature the project has learned to distrust: **rank metrics improve while calibration
does not**. ONT dRNA gains Spearman and MARD in 8 of 9 libraries with CCC flat and RMSE slightly
worse; ONT cDNA regresses CCC in 7 of 8.

### The gain is not better read assignment

Because the synthetic BAM carries per-read truth, the claim is directly testable. Measuring
true-transcript discrimination from the *final* coverage probabilities the binary actually used
(118,535 truth-labelled ambiguous reads):

| Arm | coverage picks true | combined picks true | oracle log BF |
|---|---:|---:|---:|
| `full` | 0.5213 | 0.5657 | 0.1302 |
| `endpoint-nt-measure` | **0.4779** | **0.5271** | 0.1271 |

Discrimination gets substantially *worse* (-0.043), in the same runs whose Spearman and MARD
improve. So the panel gain is not coming from assigning reads to the right transcript. It is a
systematic length reweighting -- precisely the mechanism
[coverage-signal-review-2026-07-21.md](coverage-signal-review-2026-07-21.md) identifies as able to
"change abundance without resolving read origin", and to improve an absolute length bias while
damaging fold changes.

### Why, exactly

Removing the logistic term (`--logistic-weight 0`, `--candidate-pruning none`) isolates the
endpoint model:

| Arm | coverage picks true | oracle log BF |
|---|---:|---:|
| `full`, endpoint only | 0.3800 | 0.0314 |
| `endpoint-nt-measure`, endpoint only | **0.3849** | **0.0397** |

So the measure correction *is* right on its own terms: it improves the endpoint model in isolation,
directionally reproducing the offline sweep (+0.0049 here versus +0.0158 offline, the gap being the
smoothing, gates, uniform shrinkage and Bayes cap that dilute it in situ).

But the logistic term is the dominant source of coverage discrimination in the deployed model
(0.5213 with it versus 0.3800 without), and the Jacobian's roughly `1/L^2` reweighting of the
endpoint factor fights it inside the log-linear combination. A correction that helps one factor
degrades the product.

## Decision (superseded -- see Round 6)

**Do not promote** `endpoint-nt-measure` or `endpoint-nt-measure-geometry`. It fails the
`evaluation.md` bar -- no improvement in two technology classes, cDNA regresses on a primary metric
(CCC in 7/8, MARD in 6/8) -- and the dRNA improvement it does show is traced to a length tilt
rather than to better read assignment, which is the failure mode this project has repeatedly
rejected.

Retained behind their ablation flags with the default path byte-identical.

> **This decision was reversed in Round 6 below.** It was reached by pooling the panel by
> technology, which is 27/28 comparator data. Under the truth/comparator split introduced in
> Round 5, the candidate improves every primary metric on every molecular-truth dataset while the
> comparator tier is null.

### What is now established

- The measure correction is **theoretically right and empirically right in isolation**, and wrong
  as a drop-in because it is incompatible with the logistic factor it is multiplied against. Any
  future use requires re-deriving the logistic term on the same measure, not bolting a Jacobian
  onto one factor of a product.
- **Resolution is not the bottleneck.** A warped or finer axis buys at most +0.005 discrimination
  and was correctly not built.
- **The prior is not a lever** (per-cell mass is ~0.025 pseudo-counts).
- The endpoint grid's own contribution to discrimination is modest (0.3800) next to the legacy
  per-transcript logistic term (0.5213). If a future round wants a larger effect, the logistic
  factor -- not the endpoint grid -- is where the signal actually is. Note that term is still built
  from responsibility-unweighted alignments (`oarfish_types.rs`, `add_interval(..., 1.0)`), so every
  candidate of a multi-mapping read contributes full coverage to every transcript it touches.

### Method note

Per-read truth is available without new data: the synthetic BAM encodes the source transcript in
each read name (`NM_005885_3217_aligned_...`). Alignment mode preserves read names in
`--write-coverage-signals` exports; raw-read mode does not (it writes `no_read_name_available`),
which is worth fixing if per-read validation becomes routine.

`scripts/analyze_alignment_coverage_signals.py` also crashed on every ONT export
(`ZeroDivisionError` when candidate pruning zeroes all score probabilities); fixed here.

### Evaluation trap: `--score-prob-denom` silently disables agreement calibration

`src/bulk.rs` computes `automatic_alignment_calibration` as
`... && args.score_prob_denom.is_none()`, so `--alignment-calibration auto` **abstains whenever
`--score-prob-denom` is supplied at all**. Passing the default value explicitly is therefore not a
no-op:

```
oarfish ... (no flag)                            -> quant A
oarfish ... --score-prob-denom 5.0               -> quant B  (B != A: calibration off)
oarfish ... --score-prob-denom 5.0 \
            --alignment-calibration agreement    -> quant A  (byte-identical)
```

Any sweep of the score temperature must force `--alignment-calibration agreement`, or every arm
silently also ablates the calibration promoted on 2026-07-23, and the sweep measures two changes at
once.

---

# Round 3: alignment-score temperature sweep

`--score-prob-denom` (D in `exp((score - best) / D)`, default 5.0) is the single scalar balancing
the alignment-score likelihood against coverage evidence. It had not been swept on the current
stack. Swept over D in {0.5, 1, 1.5, 2, 3, 5, 8, 12, 20} on all 28 cases with
`--alignment-calibration agreement` forced. Artifacts:
`oarfish-evaluation-data/score-temp-sweep-20260724`.

## The panel says "sharpen"; molecular truth says "flatten"

Against the LongBench Illumina comparator, every metric except Pearson improves monotonically as D
*decreases*, with no interior optimum down to D = 0.5. D = 2 versus the D = 5 default:

| Technology | Spearman | RMSE | CCC | Pearson | MARD |
|---|---|---|---|---|---|
| ONT cDNA | +0.000666 (7/8) | -13.47 (8/8) | +0.000097 (3/8) | -0.000482 (1/8) | +0.000053 (4/8) |
| ONT dRNA | +0.000636 (8/9) | -3.32 (8/9) | +0.000144 (3/9) | -0.000068 (2/9) | +0.000080 (2/9) |
| PacBio | +0.000480 (9/11) | -25.83 (10/11) | +0.001511 (7/11) | -0.000045 (6/11) | -0.000071 (9/11) |

On `synthetic-drna` -- the only case in the manifest with molecule-level ground truth -- the
direction is **exactly reversed**, also monotonically:

| D | CCC | Spearman | RMSE |
|---:|---:|---:|---:|
| 0.5 | 0.996661 | 0.841506 | 43.3 |
| 5 (default) | 0.997042 | 0.841721 | 40.8 |
| 20 | **0.997845** | 0.842182 | **34.9** |

Both objectives are boundary-seeking, in opposite directions. The current default sits between
them.

## Interpretation

The score temperature is not being optimized against accuracy; it is being optimized against
*whichever reference the panel scores with*. The LongBench comparator is matched Illumina
abundance -- a different technology with its own assignment behavior, not molecular truth -- and
sharpening the long-read likelihood evidently moves estimates toward it while moving them away from
known molecular composition.

**Recommendation: do not change the default.** More importantly, this is a property of the
evaluation harness, not of D: the 24-library Illumina-comparator mean can be driven in a direction
that molecular truth rejects. Any future single-scalar tuning should be gated on the truth-bearing
cases, not on the panel mean, and the two should be reported separately rather than pooled.

---

# Round 4: runtime and memory

Two fixed costs were identified by timing a 50k-read H69 run against GENCODE (385,659
transcripts), whose 1.6 s wall is dominated by per-reference rather than per-read work.

## `AlnInfo::prob` removed: retained (dead code), but not a memory win

The `prob: f64` field was assigned a constant at three construction sites and **never read**. It
also forced 8-byte alignment on the struct. Removing it takes `AlnInfo` from **32 to 24 bytes**
(asserted by a new unit test), avoiding 25.7 MB of allocation on a 3.2M-record sample.

Peak RSS, however, does not move. On synthetic dRNA (3 repeats, median): 347,136 KB -> 344,776 KB,
**-0.7%, within run-to-run noise** (the "before" samples span 346.6-366.4 MB, the "after" samples
are stable at 344.8 MB). Peak is set by other buffers -- the per-alignment coverage probabilities
and EM weights are each `f64` over the same records, and BGZF decompression holds its own -- so
shrinking `AlnInfo` alone does not lower the high-water mark.

Retained as dead-code removal with byte-identical output, not as a performance change. An earlier
draft of this report claimed a 10.9% saving; that compared against a single high outlier rather
than a median and was wrong.

## Deferring the seqcol digest: rejected (speed/memory trade)

`digest_from_header` runs unconditionally before a single record is parsed, and nothing consumes
its result until the run writes metadata. On GENCODE it costs **0.46 s of a 1.6 s run**. Moving it
to a worker thread (following the existing `get_digest_from_fasta` precedent) and joining at
output produced:

| Sample | reference | wall | peak RSS |
|---|---|---|---|
| H69-cdna 50k | 385k txps | 1.62 s -> **1.18 s** (-27%) | 483 MB -> **640 MB (+32%)** |
| synthetic dRNA | 33k txps | 15.5 s -> 15.4 s | +19 MB |

The digest's working set now overlaps the run's memory peak instead of being freed before it, and
that cost scales with reference size -- exactly where the speed benefit is. A 0.44 s gain for
+156 MB is not a trade this project accepts; a candidate was previously rejected over 26.7 MB.
Joining earlier does not help, because the H69 footprint is flat and reference-dominated rather
than peaking during EM.

**Reverted.** The underlying observation stands and is worth revisiting only if the digest's own
memory can be bounded (it is computed inside `seqcol_rs`, so not locally): a 27% wall reduction on
shallow samples against a large reference is otherwise available for free.

## EM convergence: examined, no action

All 28 panel cases report `converged=false` at the 1000-evaluation cap, including under SQUAREM.
This is benign. Raising `--max-em-iter` to 10000 converges at 2094 evaluations and changes **126 of
385,659 transcripts by 7.4 reads total (0.02% of the library)**, with every accuracy metric moving
in the 6th decimal. The criterion (max per-transcript relative change over entries above
`MIN_READ_THRESH`) is held above threshold by near-zero transcripts oscillating between 1.0 and
0.0 reads.

Stopping *earlier* is the real risk: L1 distance to the converged solution is 0.78% of library mass
at 100 evaluations, 0.19% at 300 and 0.017% at 1000 -- all far larger than the 5th-decimal effects
the panel is used to adjudicate. The current cap is well placed.

---

# Round 5: evaluation tooling

## Bug: raw-read mode discarded read names needed by its own export

`src/main.rs` enables read-name retention for either consumer:

```rust
.write_assignment_probs(args.write_assignment_probs.is_some() || args.write_coverage_signals)
```

but both raw-read/mapping paths derived a *second, local* flag that dropped the second clause
(`src/bulk.rs`, two sites):

```rust
let write_assignment_probs: bool = args.write_assignment_probs.is_some();
```

With only `--write-coverage-signals`, `filter_opts` therefore allocated the name vector while the
mapping workers never populated it, and every exported row was written as
`no_read_name_available`. The alignment (BAM) path has no such second flag, which is why BAM
exports carried names and raw-read exports did not.

Fixed by making both sites match `main.rs`. The independent ONT simulation now exports 99,381 rows
at 100% real read names, making per-read truth available on the raw-read path for the first time.

Immediate payoff -- a second truth-bearing dataset, which turns out not to resemble the first:

| Truth dataset | score alone | coverage alone | combined |
|---|---:|---:|---:|
| synthetic dRNA (BAM, 33k reference) | 0.4365 | 0.5213 | 0.5657 |
| independent ONT simulation (reads, GENCODE) | **0.9867** | 0.7573 | 0.9905 |

Alignment score alone resolves 98.7% of ambiguous reads on the independent simulation versus 43.7%
on the synthetic set. Any claim about "how much work the coverage model has to do" is therefore
dataset-specific.

## Truth/comparator tier split

See [evaluation.md](evaluation.md#truth-tier-versus-comparator-tier) for the protocol. Added here:

* `scripts/summarize_panel.py` -- reports per-tier means and win counts and emits a verdict,
  including `INDETERMINATE` when a run contains no truth-bearing sample. Both the rejection path
  (injected regression) and the indeterminate path are verified.
* `benchmarks/truth_tier_manifest.tsv` -- the truth-bearing panel (synthetic dRNA plus the
  independent ONT and PacBio simulations).
* `scripts/run_coverage_ablation.py` now accepts either a `bam` column or `reads` + `reference`,
  so truth-bearing panels distributed as reads run under the same driver as the BAM panel.

Applying it retrospectively to the nucleotide-measure candidate changes the reading of that
result. Split by tier rather than by technology:

| Tier | n | dCCC | dMARD | dSpearman |
|---|---:|---:|---:|---:|
| truth | 1 | +0.000351 | -0.003006 | +0.003034 |
| comparator | 27 | -0.000037 | -0.000002 | +0.000058 |

The truth tier is positive on both primary metrics while the comparator tier is null. That is the
opposite polarity to the by-technology summary that drove the original rejection, and it sits
alongside the finding that per-read discrimination *fell* 0.043 on that same truth sample. Those
are not contradictory -- abundance accuracy and read-assignment accuracy are separate axes -- but
with `n=1` in the truth tier the earlier rejection was more confident than the evidence supported.
This is precisely the ambiguity the expanded truth tier exists to resolve.

## Method note: a diagnostic that cannot work here

True-transcript discrimination (`picks_true`) is **insensitive to D by construction**: `exp((s-b)/D)`
is monotone in `s` for every `D > 0`, so the score argmax never changes. Measured across
D in {0.5, 5, 20} on 118,535 truth-labelled ambiguous reads, combined picks-true moves only
0.5656 / 0.5657 / 0.5671, and that residual comes from agreement calibration and dominance pruning,
which read score *ratios*. D changes how sharply mass is allocated in the EM, not which candidate
wins, so per-read discrimination cannot adjudicate it.

---

# Round 6: re-evaluating the nucleotide measure under the tier split

The tier split was built to stop comparator-driven gains being adopted silently. Applied to this
session's own rejected candidate, it does the opposite: it rehabilitates one.

## Truth tier, n=3, three technologies

`benchmarks/truth_tier_manifest.tsv` (synthetic dRNA, independent ONT 100k, independent PacBio
100k), `full` versus `endpoint-nt-measure`:

| Tier | n | dCCC | dMARD | dPearson | dRMSE | dSpearman |
|---|---:|---:|---:|---:|---:|---:|
| **truth** | 3 | **+0.000258 (3/3)** | -0.001075 (2/3) | +0.000253 (3/3) | -0.810 (3/3) | -0.000072 (1/3) |
| comparator | 27 | -0.000037 (4/27) | -0.000002 (9/27) | -0.000009 (9/27) | +1.094 (2/27) | +0.000058 (9/27) |

CCC, Pearson and RMSE improve on **every** molecular-truth dataset, spanning ONT dRNA, ONT cDNA and
PacBio. The comparator tier is null. Spearman is the single metric that does not improve.

## Fold changes are not damaged

The specific documented risk for a length-tilt correction is that it improves absolute calibration
while destroying condition-to-condition comparisons. Tested on the SIRV E0 -> E2 mixture at 50k:

| Arm | FC Pearson | FC Spearman | FC RMSE | Direction accuracy |
|---|---:|---:|---:|---:|
| full | 0.7497 | 0.9197 | 2.1074 | 0.986 |
| endpoint-nt-measure | 0.7493 | **0.9227** | **2.1013** | 0.986 |

Neutral to marginally better. The risk does not materialise here.

## The mechanism is bias correction, not better assignment

This remains true and is the main caveat. On both truth datasets where per-read truth is available,
the coverage term becomes a *worse* discriminator while the combined likelihood is unchanged:

| Dataset | coverage picks true | combined picks true |
|---|---|---|
| synthetic dRNA | 0.5213 -> 0.4779 | 0.5657 -> 0.5271 |
| independent ONT | 0.7572 -> 0.6700 | 0.9904 -> 0.9906 |

So the candidate does not allocate reads better; it removes a systematic length bias from the
abundance estimate. Against molecular truth that is a real gain on 3/3 datasets, and the fold-change
control is clean -- but it should be promoted as a *calibration* correction, not as an assignment
improvement, and the discrimination loss should be re-checked if the logistic/endpoint weighting
ever changes.

## Status against the promotion bar

| Requirement (`evaluation.md`) | Status |
|---|---|
| Improved median primary metrics | CCC 3/3, MARD 2/3 in the truth tier |
| Improvement in two technology classes | three (ONT dRNA, ONT cDNA, PacBio) |
| No held-out / full-length-control regression > 2% | SIRV fold change neutral; comparator tier null |
| Truth tier does not regress (Round 5 gate) | clean |
| Acceptable resource overhead | no measurable change |
| PacBio `auto` unaffected | bit-identical by the `bulk.rs` guard |

The candidate now clears the bar. It is **not** enabled by default here: `auto` still selects the
historical geometry, and promoting it is a separate decision that should be taken deliberately
rather than as a side effect of this study. The remaining reservations are the small truth tier
(n=3) and the mechanism caveat above.

---

# Round 7: the truth tier at n=12, and a flaw in the Round 5 gate

Round 6 judged `endpoint-nt-measure` on a three-sample truth tier and reported "CCC improves 3/3".
The tier now holds twelve samples -- adding the Kinnex WTC11 SIRV panel (day0 = E1, day5 = E2) and
three SIRV E0/E2 mixtures, all of which are truth-bearing and none of which were previously inside
the routine gate. The 3/3 result does not survive naively: over all twelve, CCC wins only 4/12.

## Most of the truth tier is inert for this candidate

Per-sample dCCC makes the reason obvious:

| Sample | dCCC | dMARD |
|---|---:|---:|
| synthetic-drna | +0.000351 | -0.003004 |
| independent-ont-100k | +0.000333 | -0.000158 |
| sirv-e2-drna-50k | +0.000286 | +0.000084 |
| independent-pb-100k | +0.000003 | +0.000000 |
| sirv-e0-drna-50k | 0.000000 | -0.000315 |
| sirv-e0-cdna-50k | 0.000000 | +0.000048 |
| kinnex-day0-rep1..3, day5-rep1..3 | **0.000000** | **0.000000** |

The six Kinnex rows are exactly zero *by design*: they are PacBio, and the `bulk.rs` guard keeps
PacBio `auto` on the promoted physical kernel for every endpoint-geometry ablation. Adding them was
still correct -- they will matter for any PacBio-facing change -- but they cannot move for this one.

**This exposed a flaw in the Round 5 gate.** `evaluation.md` requires improved *median* primary
metrics, and the gate computed that median across the whole tier. With half the tier structurally
inert, the median is pinned at exactly 0.000000 regardless of the candidate's merit: a real gain
and a real regression would both report as zero. `summarize_panel.py` now identifies inert samples
(exact zero on every primary metric), reports how many there are, and takes the primary-metric
median over *affected* samples only.

## Corrected reading

| Tier | inert | affected | CCC median (affected) | MARD median (affected) |
|---|---:|---:|---:|---:|
| truth | 6/12 | 6 | **+0.000145** | **-0.000079** |
| comparator | 11/27 | 16 | **-0.000062** | -0.000004 |

Among samples the change can actually reach, the truth tier improves on both primary metrics and
the comparator tier regresses on CCC -- opposite signs, which is precisely the divergence the split
was built to surface. Every truth sample with a non-trivial effect moves the right way
(+0.000351, +0.000333, +0.000286 CCC).

The candidate therefore still stands, but the honest magnitude is a **+0.000145 median CCC on six
affected truth samples**, not the "3/3" of Round 6. Round 6 overstated it by judging a
three-sample tier that happened to exclude every inert case.

## Method notes

- **Raw-read mode is non-deterministic**; see
  [evaluation.md](evaluation.md#determinism-and-a-spearman-noise-floor-in-raw-read-mode). Spearman
  moves 0.004756 between identical runs on `independent-pb-100k` while CCC/Pearson/RMSE/MARD
  reproduce exactly. The `pac-bio` truth-tier Spearman of +0.001680 above is that single noisy
  raw-read row averaged over six invariant Kinnex rows -- it is noise, not signal. BAM mode is
  byte-identical across runs.
- Judging a tier by unweighted mean or median requires knowing how many of its samples the change
  can reach. Report the inert count alongside any tier statistic.

---

# Round 8: revisiting two candidates rejected on comparator evidence

Both were flagged because their 2026-07-21 rejections rested on matched-Illumina data. Both were
re-implemented behind ablation flags and scored against the 12-sample truth tier, the 27-sample
comparator tier, and the Kinnex fold-change objective. **Both rejections stand, now on truth
evidence.**

## Continuous nested-isoform guard: rejected (corroborated)

Replaces the hard extreme-creation clamp with a continuous gate `b / (b + m)` on transcripts
co-occurring with a candidate at least 5 kb longer. Verified engaging: 4,593 protected and 503
blended on H69-pb, versus the hard clamp's zero or one.

Truth tier, per sample dCCC: **ten of twelve are exactly 0.000000**, and the two that move are
`independent-pb-100k` at **-0.005060** and `independent-ont-100k` at +0.000178 (the latter within
the raw-read Spearman/CCC noise band).

The six Kinnex rows are inert, which was not anticipated: Kinnex is PacBio, but its SIRV Set 4
reference contains no candidate pairs differing by 5 kb, so no transcript is ever *protected* and
the guard has nothing to act on. Truth-tier power depends on the mechanism under test, not merely
on matching the technology -- a caveat for the tier as a whole.

So the only truth sample that exercises this candidate regresses by 0.005 CCC, roughly sixty times
the 0.08% comparator margin the original decision was criticised for resting on. Comparator tier
agrees (CCC median over 11 affected samples -0.000431; ONT correctly inert). The original rejection
was right, and is now better supported than it was.

## Responsibility profiles: rejected (decisively)

Rebuilds the per-transcript logistic profile from coverage-free abundance responsibilities rather
than counting every alignment once, targeting the measured 5.09x ambiguous-read inflation. Verified
engaging: 5.94% of library mass moves on H69-cdna.

| Tier | dCCC | dMARD | dSpearman | dRMSE |
|---|---:|---:|---:|---:|
| truth (n=12) | **-0.014799** (4/12) | +0.014472 (1/12) | -0.025439 (0/12) | +0.034951 |
| comparator (n=27) | **-0.011287** (3/27) | +0.000310 (0/27) | -0.003183 (0/27) | +73.50 |

Regressions are one to two orders of magnitude larger than anything else measured in this session,
and both tiers agree. `independent-ont-100k` -0.0538 CCC, `independent-pb-100k` -0.0613,
Kinnex day0 -0.022 to -0.031.

### Its original selling point does not reproduce

The candidate was interesting because 2026-07-21 reported it recovering Kinnex fold change
(D Pearson 0.7526 -> 0.8025). Measured now:

| Arm | FC Pearson | FC Spearman | FC RMSE |
|---|---:|---:|---:|
| full | **0.9243** | 0.8513 | 1.3485 |
| continuous-nested-guard | 0.9243 | 0.8513 | 1.3485 |
| responsibility-profiles | **0.7668** | 0.8429 | 2.7525 |

It makes fold change substantially *worse*. The explanation is that **the baseline moved**: the
guarded physical endpoint kernel was promoted in `30a4056`, after that study. The "current auto" it
beat at 0.7526 was the pre-physical-endpoint PacBio path; today's default reaches 0.9243 without
it. The candidate was compensating for a deficiency that has since been fixed properly, and
re-applying it now damages what the fix achieved.

(These fold-change numbers use averaged replicate pairs over 69 SIRVs, not the 207-point pooled
Panel D convention of `evaluate_kinnex_fig2.py`, so they are not directly comparable to the
archived figures. The comparison *between arms* is same-code, same-data and is what matters here.)

## What this says about the revisit list

The tier split correctly identified that these two rejections rested on comparator data. It did not
follow that the rejections were wrong -- re-testing found them right, one marginally and one
emphatically. That is the intended outcome of a gate: it tells you which decisions are unsupported,
not which are incorrect.

It also adds a screening criterion that the earlier audit missed. A candidate's recorded gain is
only meaningful against the baseline it was measured on. Where the default has since improved --
as it did for PacBio in `30a4056` -- an archived gain may simply be compensating for something now
fixed. **Before revisiting any further rejected candidate, check whether the default has changed in
the relevant path since it was tested.** Applying that filter to the remaining Tier A items:

* *Joint abundance/coverage feedback* is the same family as responsibility profiles and its claimed
  gain was likewise Kinnex fold change against the pre-`30a4056` baseline. It should be considered
  closed without a re-run.
* *Equivalence-component mass conservation* (2026-07-23) postdates the physical-endpoint promotion,
  so its baseline is current and its truth-tier evidence (synthetic Pearson/CCC/RMSE all improved)
  remains untested against the wider tier. It is the one item still worth running -- though its
  0.5-1.1 s and 9-17 MiB cost is an independent objection the tier split does not address.

---

# Round 9: equivalence-component mass conservation

The last item on the revisit list, and the only one whose recorded gains postdated the
`30a4056` physical-endpoint promotion. Re-implemented as `component-mass-conservation`: union-find
over transcripts linked by shared ambiguous reads, snapshot of corrected abundances before the
abundance blend, then rescale each qualifying component so its post-blend total matches its
pre-blend total. Restricted as originally specified -- multi-transcript components only, mass at
least 0.25% of the library. Verified engaging (61 components conserved on H69-cdna) and the
union-find is unit tested, including transitive merging.

## The original result does not reproduce

| Metric (synthetic dRNA) | 2026-07-23 reported | measured now |
|---|---:|---:|
| Pearson | +0.00015 | **-0.000007** |
| CCC | +0.00016 | **-0.000006** |
| RMSE | -0.91 | **+0.036** |

Almost certainly the same baseline-drift mechanism that explained Round 8. `accuracy-refinement-2026-07-23`
evaluated *two* refinements in sequence and retained the second: agreement-calibrated alignment
likelihoods. Mass conservation was therefore measured against a baseline **without** agreement
calibration, which was then promoted on top of it. Today's `full` includes it. An implementation
difference cannot be excluded -- the 0.25% threshold could have been applied to post-blend rather
than pre-blend mass -- but the direction and the timing both point at the baseline.

## Result

| Tier | dCCC | dMARD | dPearson | dRMSE |
|---|---:|---:|---:|---:|
| truth (n=12, 6 inert) | -0.000013 (1/12) | +0.000003 (3/12) | -0.000013 | +0.003 |
| comparator (n=27) | -0.000600 (7/27) | **-0.000002 (24/27)** | -0.000392 | +11.27 |

The truth tier is a null: eight of twelve samples are exactly zero, and the largest single effect is
-0.000151 CCC on `independent-ont-100k`. The comparator tier improves MARD in 24 of 27 libraries
while regressing CCC by -0.000600 and RMSE by +11.3.

**Verdict: no promotion case.** Not a clear failure either -- it is a null with a small comparator
CCC cost and a cost in time and memory that the tier split does not excuse.

## A limitation of the Round 5 gate

`summarize_panel.py` reported `REJECT` here because the truth-tier CCC median over affected samples
was `-0.000003`. Three millionths is not a meaningful regression; it is far below the run-to-run
noise floor already documented for these samples. The gate's default `--truth-tolerance 0.0` makes
any negative median a rejection, which over-reads effects at this magnitude.

The verdict at this scale should be read as **null**, not as a demonstrated regression. The
tolerance needs calibrating against measured noise -- per metric and per sample class -- rather
than left at zero. Until then, treat a verdict whose margin is smaller than the noise floor as
"no effect detected".

## Closing the revisit list

All four Tier A candidates identified by the tier-split audit are now closed:

| Candidate | Outcome |
|---|---|
| Continuous nested-isoform guard | rejected; -0.005060 CCC on the one truth sample that exercises it |
| Joint abundance/coverage feedback | closed without re-run; same family and same superseded baseline as responsibility profiles |
| Responsibility profiles | rejected decisively on both tiers; original fold-change gain does not reproduce |
| Equivalence-component mass conservation | null; original synthetic gain does not reproduce |

The audit was still worth running: it converted four decisions resting on comparator evidence into
one rejection on truth evidence, two reproduction failures traced to a moved baseline, and one null.
None of the original decisions was overturned.

The broader lesson is that **baseline drift, not comparator bias, was the dominant confound in the
archive.** Three of the four candidates had recorded gains measured against a default that has since
improved. Any future audit of rejected work should check what the default did in the relevant code
path between the original test and now, before spending a run.

---

# Round 10: recovery probability and poly(A) — data assessment

## Truth semantics determine what is measurable

The simulation truth files total exactly 100,000 for a 100,000-read sample: they are **read counts,
not molar abundance**. Capture bias is already inside them, so a recovery probability is 1 by
construction there and any length-dependent residual is *assignment* error. Only molar-truth data
(SIRV, Kinnex) can measure recovery. This distinction decides which datasets can evaluate which
model, and is worth checking before designing any normalization term.

Residual on the simulations (aggregate share per length class, not per-transcript medians, which are
degenerate when most transcripts have a truth of one read):

| Length | ONT sim | PacBio sim |
|---|---:|---:|
| <500 | -0.060 | -0.035 |
| 500-999 | +0.001 | -0.003 |
| 1000-1999 | +0.005 | +0.005 |
| 2000-3999 | +0.020 | +0.014 |
| >=4000 | +0.034 | +0.024 |

Small, monotonic, both technologies: oarfish under-assigns short transcripts. That is an assignment
bias, and it is the same direction that `endpoint-nt-measure` corrects.

## Recovery bias is real, large, and upstream

Measured on Kinnex against SIRV molar concentrations (63 SIRVs with reads, 6 replicates, two mixes):

| Length | day0 (E1) | day5 (E2) |
|---|---:|---:|
| 161-399 | **-11.79** | **-12.88** |
| 400-699 | -4.34 | -4.31 |
| 700-999 | +0.55 | +0.71 |
| 1250-2498 | +0.96 | +0.88 |

Two SIRVs at *identical* molarity: SIRV618 (189 nt) gets **0 reads**; SIRV703 (2498 nt) gets
**48,352**. Short molecules are essentially absent below 400 nt and ~20x under-represented at
400-699 nt. This is far more severe than the "mild bias against transcripts shorter than 1.25 kb"
reported for Kinnex, and it is reproducible across all six replicates.

The same measurement on ONT dRNA SIRV is **flat** (-0.20 to +0.16), so this is technology-specific,
not universal.

**Implication.** Those reads do not exist; the loss is library preparation and size selection,
entirely upstream of oarfish. No inference recovers a molecule that was never sequenced. A recovery
term would therefore be a *reporting-stage* conversion from read counts to estimated input molarity:
estimable only where spike-ins are present, undefined exactly where the bias is worst (rho -> 0),
and unable to improve read assignment. It is a real phenomenon but a narrow lever.

## poly(A): the data exists locally, the headroom looks small

The local `hek_polyA.fastq` is synthetic -- byte-for-byte `hek_100k.fastq` with exactly 200 A's
appended to every read -- and would manufacture poly(A) evidence at positions that are not real
transcript 3' ends. SG-NEx fastq are effectively trimmed (median longest A-run 5, i.e. background).
ONT open data's RNA sets are synthetic 5-mer oligos.

But the LongBench transcriptome BAMs **retain the tail as soft-clipped sequence**: reads whose
terminal soft clip contains an A/T run of at least 12 nt are 50.7% (dRNA), 65.7% (cDNA) and 98.3%
(PacBio) of primary alignments.

The aggregate discriminative test is discouraging. Comparing the transcript-oriented 3' gap of reads
with and without such a clip:

| Sample | with poly(A) clip, 3' gap <=20 nt | without |
|---|---:|---:|
| H69 dRNA | 59.2% | 57.3% |
| H69 cDNA | 53.6% | 56.6% |
| H69 PacBio | 48.9% | 9.9% (n=739) |

On ONT the presence of a poly(A) clip barely predicts 3' completeness, most likely because dRNA
reads are *already* 3'-anchored (median 3' gap 3 nt, 41% exactly zero), so the existing endpoint
model already captures what poly(A) would tell us. On PacBio the separation is large but 98% of
reads carry the clip, leaving almost no contrast to exploit.

This does not refute the mechanism as stated -- ruling out *candidates* whose 3' end disagrees is a
per-read, per-candidate question, and only the aggregate marginal was tested here. But the headroom
is smaller than it first appears, and a sharper tail detector than "longest A/T run in the terminal
clip" would be needed before the per-candidate test is worth running.

---

# Round 11: per-candidate poly(A) test, and locating spike-in data

## poly(A): the per-candidate test is positive

Round 10 tested the wrong thing. Asking "does a poly(A) clip predict 3' completeness on average"
gave a null, but the mechanism is per-candidate: for an ambiguous read known to be 3'-complete, do
its candidates *disagree* about 3' completeness?

Detector sharpened so the tail must sit immediately adjacent to the alignment end -- forward strand,
first 30 nt of the terminal soft clip, >=80% A; reverse strand, last 30 nt of the leading clip,
>=80% T -- and candidate sets taken from name-collated groups (secondary records carry no SEQ, so
tail status is inherited from the primary).

| Sample | ambiguous reads with poly(A) | have a flush AND a far candidate | candidates rulable |
|---|---:|---:|---:|
| H69 dRNA | 14,347 | **41.1%** | 3.1 of 7.6 |
| H69 cDNA | 7,063 | **49.9%** | 2.8 of 7.2 |
| H69 PacBio | 24,123 | **36.9%** | 4.7 of 9.3 |

(flush = 3' gap <= 20 nt; far = 3' gap > 100 nt.)

So on 37-50% of poly(A)-bearing ambiguous reads, roughly 40% of the candidate set has a 3' gap over
100 nt and is incompatible with a read that demonstrably reached its poly(A). The reads-without-
poly(A) control shows a similar geometry (35.5 / 54.5 / 26.0%), as expected -- candidate geometry is
a property of the annotation, not of the tail -- so the gain comes from *knowing which reads are
3'-complete*, not from those reads being unusual.

**Why the endpoint model does not already capture this.** A large 3' gap is only mildly penalised
today, because most reads legitimately are truncated. Poly(A) converts "possibly truncated" into
"demonstrably not truncated", which licenses a far sharper penalty for exactly those reads. That is
a conditional sharpening the current marginal endpoint model cannot express.

**Risk to control.** Alternative polyadenylation and wrong annotated 3' ends would cause a
hard constraint to eliminate the true candidate. This must be a strong likelihood ratio, not a hard
filter.

## Spike-in data for rho: located, validated, and thinner than hoped

The SG-NEx ONT transcriptome BAMs are aligned to `hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa`, which
carries the spike-in panel we need:

| Family | n | length range |
|---|---:|---|
| sequin | 165 | 283 - 6,943 |
| SIRV core | 69 | 161 - 2,498 |
| long SIRV | 12 | 3,997 - 12,029 |

Combined coverage by length class: 22 (<500), 94 (500-1k), 84 (1k-2k), 33 (2k-4k), 13 (>=4k) --
a usable spread from 161 nt to 12 kb, far wider than the 161-2,498 nt our Kinnex data reaches.

Two practical findings for anyone using this:

* **Spike-ins are reference indices 0-245**, i.e. the *front* of a coordinate-sorted BAM, so a small
  range request from the file start retrieves them without downloading 752 MB.
* **Most samples contain no spike-ins at all**, and the path names do not say so. Probing the first
  record's refID is a definitive one-request test. Of six direct-RNA samples probed, only
  `SGNex_K562_directRNA_replicate4_run1` (first refID 0) and `SGNex_H9_directRNA_replicate3_run1`
  (first refID 165) carry them; A549 replicates 1/5/6, HepG2 replicate 2 and MCF7 replicate 3 do not.

Depth is the limitation. A 120 MB prefix of the K562 sample yields 572 primary sequin reads across
44 of 246 spike-in references, median 3-16 reads per reference. That is too thin to fit rho(L) from
one sample; a usable fit needs pooling across the spike-in samples and probably the cDNA libraries,
plus the Garvan sequin concentration tables, which are not in the SG-NEx S3 annotations directory.

---

# Round 12: poly(A) 3'-completeness term — implemented, unconfirmable

## Implementation

`three_prime_polya` on the `AlnRecordLike` trait keeps sequence access inside the record
abstraction. The tail must sit *immediately adjacent* to the alignment end: forward strand takes the
first 30 nt of the terminal soft clip and requires >=80% A; reverse strand takes the last 30 nt of
the leading clip and requires >=80% T, since SEQ is reference-oriented. Only the primary record
carries SEQ, so a read's tail is the maximum over its group; the primary is always retained because
it holds the best score. Raw-read mode returns 0 -- the rammap mapping does not retain the sequence
at that point -- so this is BAM-only.

In-binary detection matches the offline analysis to within a percent (dRNA 44.6%, cDNA 20.3%,
PacBio 85.8%), confirming the port is faithful.

`src/util/polya_probability.rs` learns a log-spaced histogram of 3' gap in nucleotides from
*single-candidate* tail-bearing reads, compares bins as densities (divided by bin width), scores
ambiguous tail-bearing reads against it, and caps the odds at 20x so annotation error cannot
eliminate a true candidate. Reads with no detected tail are untouched: absence of a tail is not
evidence of truncation. Everything is sample-learned.

The term is strong where it applies -- mean log odds range **2.05 (dRNA) / 2.32 (cDNA)** against the
existing coverage term's 0.80 / 0.99, so two to three times the entire current coverage signal.

## Result

| Tier | dCCC | dMARD | dSpearman | dRMSE |
|---|---:|---:|---:|---:|
| comparator (n=27) | **+0.002457** (19/27) | **-0.000336 (27/27)** | **+0.003401 (27/27)** | -37.15 (18/27) |
| truth (n=12) | -0.000008 (7/12) | -0.000003 (7/12) | -0.000446 (2/12) | -0.000087 |

The comparator gain is the largest and most consistent of the entire study -- 27/27 on both MARD and
Spearman, improving in all three technologies (PacBio CCC +0.00495, dRNA +0.00132, cDNA +0.00017).

## The truth tier cannot evaluate this candidate

Tail exposure across the whole truth tier is roughly **700 reads**:

| Sample | tail reads | why |
|---|---:|---|
| independent-ont-100k, independent-pb-100k, all 3 SIRV | **0** | raw-read mode; detector returns 0 by construction |
| kinnex-day0/day5 x3 | 10-24 each (of ~1M) | tails trimmed upstream of the provided BAMs |
| synthetic-drna | 605 (of 1.36M) | simulated without tails |

So "truth tier clean" here means "truth tier measured nothing". The verdict is uninformative, not
supportive. This is a **coverage gap in the truth tier**: it cannot evaluate any BAM-only signal,
because five of its twelve samples are raw-read and the rest do not retain tails.

Nor is a truth-anchored discrimination test available: the only per-read-truth BAM
(`synthetic-drna`) has 605 tail reads, and the simulations that do have per-read truth run in
raw-read mode.

The SG-NEx spike-in samples were checked as a candidate fix. `SGNex_K562_directRNA_replicate4_run1`
retains tails in only **5.6%** of primary reads by this criterion (against 44.6% for LongBench
dRNA), and carries just 572 primary sequin reads in a 120 MB prefix. Not enough on either axis.

## Decision: do not promote

A large, perfectly consistent comparator gain with zero truth evidence is the exact configuration
this study has learned to distrust -- it arose three times (unique-read profiles, the nucleotide
measure, the score temperature) and truth contradicted the comparator every time.

There is also a concrete mechanism by which this gain could be spurious. Illumina poly(A)-selected
libraries are themselves 3'-biased. Sharpening how we model 3' completeness could move our estimates
*toward Illumina's own bias* without making them more correct -- the same failure mode as the score
temperature, where the comparator wanted `D -> 0.5` and molecular truth wanted `D -> 20`.

Retained behind `--coverage-ablation polya-three-prime`, default byte-identical.

## What would settle it

Truth-bearing, BAM-mode data that retains poly(A). None of our current sources qualifies. The
cleanest options, in order of cost:

1. An untrimmed dRNA or cDNA run over a spike-in sample, aligned to the transcriptome and
   name-collated. The tails must survive the basecalling and alignment pipeline, which most public
   deposits do not preserve.
2. Dorado `pt:i` tags, which estimate the tail from raw signal rather than basecalls and are
   therefore robust to trimming. This would also let raw-read mode participate. Requires uBAM
   deposits, which SRA's fastq normalization discards.
3. Simulation with realistic tails *and* realistic alternative polyadenylation. Note that naive
   simulation would place every tail exactly at the annotated 3' end and inflate the result: the
   measured `diffuse_fraction` of 0.32-0.38 says roughly a third of real unambiguous tail-bearing
   reads sit more than 256 nt from their annotated 3' end.

---

# Round 13: spike-in data with retained tails — built, run, and underpowered

## The data

Surveying SG-NEx transcriptome BAMs for samples that have *both* spike-ins and untrimmed poly(A)
found them among the cDNA libraries, not direct RNA:

| Sample | spike-ins | tail rate |
|---|---|---:|
| SGNex_HEYA8_cDNA_replicate1_run4 | yes | **39.6%** |
| SGNex_HEYA8_cDNA_replicate3_run3 | yes | 38.9% |
| SGNex_H9_cDNA_replicate3_run4 | yes | 30.9% |
| SGNex_K562_directRNA_replicate4_run1 | yes | 5.6% |
| A549 / HepG2 / MCF7 direct RNA and directcDNA | mostly none | 5-14% |

`SGNex_HEYA8_cDNA_replicate1_run4` gives **158,269 primary sequin reads at 30.3% tail retention**
across 77 spike-in references (283-6,943 nt), with five sibling replicates available for pooling.

Extraction: spike-ins occupy reference ids 0-245, so a 150 MB prefix of the coordinate-sorted BAM
contains their alignments. `extract_spikein.py` rewrites them name-collated, keeping exactly the
first 246 references so record refIDs stay valid, and emitting a minimal BGZF stream. oarfish parses
the result directly (247,089 records / 158,746 reads; EM converged in 51 evaluations). Sequin
concentrations come from `rnasequin_isoforms_2.4.tsv`; the mixture is identifiable empirically
(log-log Pearson 0.5875 against MIX_A versus 0.4917 against MIX_B, so SG-NEx HEYA8 is Mix A).

## The detector is validated

| | LongBench (real human) | sequins |
|---|---:|---:|
| `flush_fraction` | 0.54-0.59 | **0.9927** |
| `diffuse_fraction` | 0.32-0.38 | **0.0031** |

On synthetic constructs with defined 3' termini, 99.3% of unambiguous tail-bearing reads are flush
and 0.3% are diffuse. The detector therefore does **not** generate false positives; the 32-38%
diffuse rate on human data is real alternative polyadenylation and mis-annotated 3' ends. That is
worth knowing independently of this candidate.

## Result: null on molar truth

| Arm | Pearson | Spearman | CCC | RMSE | MARD |
|---|---:|---:|---:|---:|---:|
| full | 0.50507 | 0.40958 | 0.50323 | 254.061 | 0.70165 |
| polya-three-prime | 0.50506 | 0.40961 | 0.50322 | 254.067 | 0.70261 |

Deltas are 5th-decimal with mixed signs, despite the term firing hard (mean log odds range 1.74 on
17,915 scored reads).

## Why: spike-ins under-represent the ambiguity that matters

| | avg candidates per ambiguous read | discriminating | rulable per read |
|---|---:|---:|---:|
| sequins | **2.8** | 27.9% | 1.4 |
| LongBench dRNA | **7.6** | 41.1% | 3.1 |

Sequin equivalence classes are 2.7x simpler than human ones. Worse, dominance pruning had already
removed 50,316 of 245,485 candidates in this sample, so the candidates poly(A) would eliminate were
largely eliminated already.

**So the sequin panel cannot adjudicate this candidate.** A null there does not refute the
comparator gain, because the regime the term operates in -- large, structurally diverse candidate
sets -- barely exists in the spike-in data.

This generalises beyond poly(A), and is the main lesson of this round: **synthetic spike-ins are
excellent for calibration questions (recovery, length bias, absolute abundance) and weak for
assignment questions**, because their ambiguity structure is far simpler than a real transcriptome.
A truth tier built only from spike-ins and simulations will systematically under-test read
assignment. Choosing a truth instrument requires matching its *ambiguity structure* to the mechanism
under test, not merely having ground truth.

## Status

`polya-three-prime` remains unpromoted and opt-in, default byte-identical. The evidence now is:
detector validated; mechanism engages strongly; large consistent comparator gain (27/27 MARD and
Spearman); no truth instrument capable of confirming or refuting it. Settling it needs truth-bearing
data with *human-like* isoform ambiguity -- which in practice means per-read truth on a real
transcriptome, i.e. simulation with realistic APA, rather than spike-ins.

The same data does serve the rho(L) question well, since recovery is exactly the calibration
question spike-ins are good at: 77 references spanning 283-6,943 nt with known molarity, six
replicates, in a real cDNA library.

---

# Round 14: rho(L) from sequin spike-ins

## A truncation artifact, and how it was caught

The first attempt used the same 150 MB BAM prefix built for the poly(A) test and produced an
erratic, non-monotonic curve:

| length | log2 rho (INVALID) |
|---|---:|
| 0-500 | +0.300 |
| 500-750 | +1.184 |
| 750-1000 | -1.232 |
| 1000-1500 | +0.265 |
| 1500-2500 | -0.601 |
| 2500+ | -0.358 |

Adjacent bins swinging 1.2 log2 is not a recovery function, and quantification accuracy against
MIX_A was poor (Pearson 0.505). Cause: **the prefix reached only reference id 84 of 245**, because
BAM records here carry SEQ and QUAL and average roughly 600 bytes, so the spike-in block is far
larger than the 150 MB assumed. 161 of 246 spike-in references therefore had zero reads *by
truncation*, and their obs/exp ratios were meaningless.

A prefix is safe for asking "does this sample contain spike-ins and do its reads retain tails"
(both were answered correctly) but not for quantification, which needs every reference's reads. The
check is cheap: track the highest reference id reached and require it to exceed the spike-in block.

The invalid numbers are recorded here rather than deleted, because the failure mode -- a
coordinate-sorted prefix silently starving later references -- is easy to repeat.

## Getting the data right

The BAM is **33.29 GB**, not the 752 MB quoted earlier (that figure came from the A549 sample and
was never checked for this one). Downloading it whole is unnecessary: the `.bai` index gives per
reference virtual offsets, and parsing it shows the spike-in block (references 0-245) ends at
**704.3 MB** -- a 47x reduction, and *derived* rather than guessed, which is exactly what the first
attempt got wrong.

A second false alarm is worth recording. The extraction reached reference id 242, not 245, and the
coverage check flagged truncation. The index shows references 243, 244 and 245 have **zero chunks**,
i.e. no alignments anywhere in the file; 25 spike-in references are empty in total, and
246 - 25 = 221 is exactly the number extracted. Coverage was complete. The check should compare
against the highest reference that the *index* says has alignments, not against the last reference id.

With complete coverage, quantification against MIX_A rises from Pearson **0.5051** (truncated) to
**0.8825**, confirming the earlier numbers were an artifact and these are sound.

## Result: no useful length-dependent recovery bias in ONT cDNA

158 sequins, 276,748 assigned reads.

| Length | n | obs share | exp share | log2 rho |
|---|---:|---:|---:|---:|
| 0-500 | 9 | 0.0182 | 0.0213 | -0.228 |
| 500-750 | 26 | 0.2838 | 0.1973 | +0.524 |
| 750-1000 | 34 | 0.3430 | 0.3543 | -0.047 |
| 1000-1500 | 36 | 0.1838 | 0.1693 | +0.118 |
| 1500-2500 | 40 | 0.1563 | 0.2477 | -0.664 |
| 2500+ | 13 | 0.0149 | 0.0100 | +0.570 |

Non-monotonic. The weighted regression on log2 length over the 94 sequins with at least 20 reads
gives a slope of **-0.234 log2 per doubling** with a **residual standard deviation of 0.973** -- the
per-sequin scatter is roughly four times the trend across the whole length range. There is at most a
weak tendency for shorter transcripts to be *better* recovered, and it explains little.

## Recovery bias is protocol-specific, not a general long-read property

| Protocol | Instrument | Result |
|---|---|---|
| PacBio Kinnex | SIRV, 6 replicates | **-11.8 to +0.95 log2, monotonic**; 189 nt gets 0 reads where 2498 nt gets 48,352 at equal molarity |
| ONT direct RNA | SIRV | flat (-0.20 to +0.16) |
| ONT cDNA | sequins, 158 refs, 283-6,943 nt | slope -0.23/doubling, swamped by 0.97 scatter |

**A general rho(L) term is therefore not warranted.** For ONT there is no meaningful length-dependent
recovery to correct. For Kinnex the bias is real and enormous, but it is upstream size selection --
the reads do not exist -- so rho -> 0 below 400 nt and no inference recovers them. That was the
conclusion reached provisionally in Round 10 on much weaker evidence; it now rests on a wide-range
spike-in measurement in a real library.

One observation worth keeping for later: the ~1 log2 per-sequin residual scatter that length does
*not* explain implies roughly two-fold recovery differences between transcripts of similar length.
Whatever drives that -- GC content, secondary structure, sequence composition -- is a larger effect
than length, and is a covariate nobody here is modelling.

## poly(A) on sequins, re-run on corrected data

The Round 13 poly(A) evaluation used the truncated extraction, so it was re-run:

| Arm | Pearson | Spearman | CCC | RMSE | MARD |
|---|---:|---:|---:|---:|---:|
| full | 0.88252 | 0.96318 | 0.88252 | 118.124 | 0.37944 |
| polya-three-prime | 0.88245 | 0.96322 | 0.88245 | 118.168 | 0.37810 |

Still a null, now on a sound baseline. The Round 13 conclusion stands -- and stands for the right
reason: the spike-in panel's equivalence classes are too simple to test the mechanism, not because
the data were broken.

---

# Round 15: sequence-driven recovery variation — does not replicate

The 1 log2 per-sequin recovery scatter that length does not explain motivated testing sequence
covariates. On HEYA8 cDNA rep1 sequins (n=94, read-weighted), 5' GC looked like a strong lead:

| Predictor | weighted R2 | coefficient |
|---|---:|---:|
| log2 length alone | 0.031 | -0.234 per doubling |
| **GC of 5' 200 nt alone** | **0.218** | **-4.30 log2 per unit GC** |
| both | 0.299 | length -0.39, GC5 -4.89 |

Correlations were 5'-asymmetric -- GC of the 5' 200 nt at -0.51 against GC of the 3' 200 nt at
-0.12 -- which is what reverse-transcription processivity would predict, since RT initiates at the
3' poly(A) and must traverse the molecule to reach the 5' end. `internal_A6_runs` at +0.34 was
consistent with internal priming producing easier-to-complete products.

## It reverses on an independent molecule set

| Molecules | Sample | R2 GC5 | GC5 coefficient |
|---|---|---:|---:|
| sequins (160) | HEYA8 cDNA rep1 | 0.218 | **-4.30** |
| SIRVs (69) | HEYA8 cDNA rep1 | 0.068-0.142 | **+7.59 to +9.20** |
| SIRVs (69) | H9 cDNA rep3 | 0.001-0.116 | **+1.27 to +6.56** |

The sign flips, in the *same library*, on a different spike-in panel. H9 could not be used for the
sequin replication at all: it carries only 122 sequin reads (34 references, none above 20 reads)
against 135,987 SIRV reads, so it was spiked with SIRVs and not sequins.

## Conclusion: unusable as a covariate

A coefficient that reverses between spike-in panels cannot be used to model real transcripts,
whichever panel is closer to the truth. The sequin correlation is most plausibly a property of that
panel's sequence design rather than a transferable recovery effect.

The SIRV test is not a clean replication either, and this should be stated: SIRVs are deliberately
overlapping isoforms, so their observed/expected ratios are dominated by *assignment* error, which
could mask or invert a genuine recovery effect. Neither panel is a clean instrument for
sequence-driven recovery -- sequins because their design may correlate GC with something else,
SIRVs because their ambiguity swamps recovery. The stark sign reversal is nonetheless not what a
shared underlying effect would produce.

What would resolve it: sequins in a second, independent library (HEYA8 replicate 2 or 3 both carry
them), which separates library-level reproducibility from molecule-set specificity. That test is
worth roughly one 700 MB download if the question is revisited. Even a positive result there would
only establish the effect *for sequins*, not that it transfers to a real transcriptome.

**Status: lead closed.** The ~1 log2 unexplained recovery scatter is real and larger than the length
effect, but no covariate tested here explains it in a way that transfers.

---

# Round 16: annotation-omission modeling — the signal is not in transcriptome space

Annotation omission is the largest documented error source in long-read quantification, and oarfish
models none of it: reads from unannotated isoforms are force-assigned to annotated ones. The natural
signal seemed to be *unexplained clipping* -- read sequence that could have aligned to the assigned
transcript but did not -- which `censoring_probability::unexplained_clip` already computes.

## Validation: it does not work

Using the synthetic dRNA BAM, whose read names encode the true source transcript, the model's actual
decision was simulated per read: how much worse is a read's best fit when its true transcript is
removed from the candidate set?

Over **221,273** reads that have both their true transcript and at least one alternative:

| Quantity | Result |
|---|---|
| best unexplained clip **with** the true transcript | median 95 nt |
| best unexplained clip **without** it | median 97 nt |
| per-read increase on removal | **median +0, p75 +0, p90 +0** |
| reads where removal raises the clip by >30 nt | **4.3%** |
| reads where removal changes nothing at all | **93.0%** |

For 93% of reads, deleting the true source transcript does not degrade the fit by a single
nucleotide.

## Why: transcriptome-space alignment discards splice structure

Isoforms of a gene share most of their sequence. A read from isoform A, with A removed, aligns to
sibling B across the shared exons with an identical clipping profile. The evidence that would betray
the omission is the read's *junction combination*, and transcriptome-space alignment cannot
represent it: the read is aligned contiguously to a linear transcript, so a disagreeing splice
structure simply looks like a normal alignment.

This is a property of the alignment space, not of the statistic chosen. No transcriptome-space
covariate can recover it, which retires the whole family of transcriptome-mode approaches to this
problem.

## The feature belongs in genome mode

In genome mode oarfish spliced-aligns to the genome and projects through bramble, so junctions
survive. Two pieces already exist:

* **Projection-loss counters** (`bulk.rs`, `n_genome_mapped` vs `n_projected`): reads that align to
  the genome but yield no transcriptome projection. On the parity dataset this is only 1,366,728 vs
  1,356,492, i.e. **0.75%** -- the extreme cases, where no transcript is compatible at all.
* **Junction-mismatch scoring**: bramble already counts internal junction mismatches, which is what
  the hidden `--junc-miss-discount` (default 1.0, off) discounts similarity by.

The larger and more interesting pool is not projection failure but reads that *do* project while
disagreeing with every candidate's junction structure. Those are direct evidence of an unannotated
isoform, and they are currently absorbed silently.

**Enabling change:** `ProjectedAlnRecord` carries only the discounted `similarity`; bramble does not
return the junction-mismatch *count*. Exposing that count per projection is the prerequisite. With
it, a read whose every candidate has at least one internal junction mismatch can be routed to a
per-locus novel state, and the accumulated mass reported as a per-gene "unexplained" diagnostic.

## Consequences for the oarfish 2 plan

1. Annotation-omission modeling is a **genome-mode feature**, and the paper must say so. That is
   defensible -- it is a real advantage of projection-based quantification over transcriptome-only
   tools, and it gives genome mode a reason to exist beyond convenience.
2. The proposed holdout demonstration (delete a fraction of expressed isoforms, show competitors
   silently redistribute the mass) remains the right experiment, but must be run in **genome mode**.
   In that setting it should work by construction: deleting a transcript deletes its junction
   combination, so its reads mismatch every survivor.
3. The transcriptome-mode path should be dropped from the feature rather than attempted.

Cost of the validation that produced this: one analysis over an existing BAM. It would have been an
expensive thing to discover after building the model.

---

# Round 17: junction evidence exposed, and the omission signal validated

## bramble change

Local bramble synced to `origin/main` (9 commits, including a fix for non-deterministic primary
alignment on score ties). `bramble-rs` already computed internal junction agreement -- it is what
`--junc-miss-discount` penalizes -- but only the *discounted similarity* reached callers.

`ProjectedAlignment` now also carries:

* `junc_misses: i32` -- internal boundaries where the read disagrees with this transcript.
* `junc_hits: i32` -- boundaries where it agrees.

Both are needed. A single-exon read, or one contained within an exon, spans no internal boundary and
has `junc_hits == 0 && junc_misses == 0`; it is uninformative about splice structure and must be
distinguishable from a read that spans four boundaries and matches them all. Reporting misses alone
would conflate the two.

oarfish is pointed at the local checkout by path while this is developed; the crate is published
later. Note the counts feed `similarity_score` through `junc_miss_discount`, so a consumer using
both should avoid double-counting.

Two unit tests lock the fields. Note bramble's `tests/evaluate.rs` does **not** compile at upstream
HEAD -- it references `ReadEvaluationConfig` fields renamed by the "updated user-facing parameters"
commit (`max_ins`, `small_exon_size`). Nine errors, all pre-existing and unrelated to this change;
the library tests pass (9/9 with the additions). Worth fixing before publishing.

## Validation: the signal responds to annotation completeness

oarfish now reports, per run, how many reads span at least one internal boundary and how many of
those disagree with *every* candidate transcript. Measured on the parity genome BAM (1,356,492
aligned reads) against progressively degraded annotations:

| Annotation | junction-informative reads | disagree with every candidate |
|---|---:|---:|
| complete (71,220 transcripts) | 899,200 | 7,876 (**0.88%**) |
| 10% of transcripts removed | 846,804 | 11,802 (**1.39%**) |
| 25% of transcripts removed | 717,258 | 16,691 (**2.33%**) |

The parity reads are simulated *from this annotation*, so with the complete GTF there are no
unannotated isoforms by construction and **0.88% is a false-positive rate** -- the cost of junction
tolerances and alignment noise. Removing annotation raises the detected mass monotonically, 2.6x at
25% holdout.

## Interpreting the sensitivity

Deleting 25% of transcripts raises detected mass by only 1.45 points, so most reads from a deleted
isoform are *not* flagged. That is correct behaviour rather than a shortfall: when isoform A is
removed, most of its reads remain fully compatible with sibling B's junction structure, and a read
compatible with an annotated isoform is not evidence of an annotation gap. The feature detects
*unexplainable* reads, not *misassigned* ones -- and only the former justifies a novel state.

This also explains why the transcriptome-space test in Round 16 found nothing: there, *no* read is
ever unexplainable, because splice structure is invisible and any sibling absorbs it silently. In
genome space a residual 1-2% is recoverable, with a low false-positive floor.

## Remaining work

1. Route the flagged mass to a per-locus novel state in the EM rather than only counting it, so it
   stops being misattributed to annotated isoforms.
2. Attribute it per gene/locus and emit a per-gene "unexplained mass" diagnostic -- the user-facing
   deliverable.
3. Quantify the accuracy benefit under holdout: does absorbing this mass improve estimates for the
   *surviving* isoforms, which is the claim that matters.
4. Fix `tests/evaluate.rs` upstream and publish the bramble crate; revert oarfish's path dependency.

---

# Round 18: the novel latent state, and a dose-response

## Model

Reads flagged in Round 17 now receive an extra EM candidate: a per-locus *novel* latent state
representing a hypothetical unannotated isoform that would explain the read's splice structure.

* **Locus** = the read's ambiguity component (transcripts linked by shared ambiguous reads), reusing
  `ambiguity_components`. No gene annotation is required.
* **Weight** = the best annotated candidate's weight scaled by `novel_odds_per_miss ^ misses`
  (default 2.0 per unmatched junction). Each unmatched internal junction is one doubling of evidence
  that the true isoform is absent.
* Novel states occupy indices past the annotated transcripts in the same counts vector, so consumers
  indexing by reference id are unaffected and the states are simply absent from the output.

Two details mattered for correctness. Single-candidate reads normally short-circuit the M-step
because their weight cancels; that is no longer true once a novel alternative competes, so the
short-circuit is now conditional. And exact equivalence-class collapsing had to include the novel id
and weight in its hash, or reads with different novel states would be merged.

The parallel M-step does not implement the novel term, so the serial path is forced when the feature
is active.

## Dose-response on the surviving transcripts

The claim that matters is not "we detect unannotated mass" but "absorbing it improves the estimates
for the isoforms that *are* annotated". Parity genome BAM, accuracy against molecular truth
restricted to transcripts surviving each holdout:

| Holdout | surviving n | dPearson | dCCC | dRMSE | dMARD |
|---|---:|---:|---:|---:|---:|
| **0% (complete)** | 23,809 | **-0.00000** | **+0.00000** | **+0.00** | +0.00104 |
| 10% | 21,468 | +0.00029 | +0.00036 | -0.55 | +0.00075 |
| 25% | 17,898 | +0.00069 | +0.00077 | -0.72 | +0.00077 |

Two properties, both required:

1. **Neutral when the annotation is complete.** Pearson, CCC and RMSE are unchanged to five decimals
   at 0% holdout. A model that fires on complete annotations would be unusable regardless of its
   benefit elsewhere.
2. **Benefit grows monotonically with incompleteness**, roughly doubling from 10% to 25% holdout.

At 25% holdout the novel states divert 2,422 reads off surviving transcripts.

## The blemish: MARD

MARD is consistently *worse* by about +0.001, and worst at 0% holdout where nothing should change.
That points at the 0.88% false-positive reads: diverting correctly-annotated reads costs MARD, which
weights low-abundance transcripts heavily, more than it costs the correlation metrics. Since MARD is
a primary metric under `evaluation.md`, this must be resolved before promotion. Two levers are
untested: `--novel-odds-per-miss` (default 2.0), and requiring more junction evidence -- for example
at least two unmatched junctions, or a minimum `junc_hits` -- before a read is eligible.

## Status and remaining work

Implemented behind `--coverage-ablation annotation-omission`; transcriptome-mode output is
byte-identical; 53 tests pass.

1. Tune eligibility and odds to remove the MARD regression.
2. Emit per-locus novel mass as a user-facing per-gene "unexplained mass" report -- the deliverable
   that distinguishes this from a silent internal correction.
3. Validate on real data with genuinely incomplete annotation, not only simulated holdout. The
   dose-response above is a controlled experiment on simulated reads; it demonstrates the mechanism,
   not the real-world magnitude.
4. Fix bramble's `tests/evaluate.rs`, publish the crate, revert oarfish's path dependency.

---

# Round 19: bramble test repair, and tuning out the MARD regression

## bramble `tests/evaluate.rs` fixed

The breakage came from `bce5396` ("updated user-facing parameters for clarity"), which renamed
`max_ins` -> `max_junc_ins` and folded `small_exon_size` into the pre-existing `max_error_exon`
without updating the integration test. The asserted *values* were all still correct (short read
0/0/0, long read 40/40/35), so only the field names needed updating. bramble's suite is now green:
9 + 3 + 1 + 7 tests, zero compile errors.

## The MARD regression: diagnosis

The damage was concentrated exactly where diverting a read is most costly in relative terms:

| Truth abundance | mean change in absolute relative difference |
|---|---:|
| <= 5 reads | **+0.00207** |
| > 50 reads | +0.00001 |

Losing one read from a 3-read transcript is a large relative error; from a 500-read transcript it is
nothing. At complete annotation every such diversion is a false positive.

## Two levers that did not work, and why

**Per-read stringency fails.** Requiring at least two unmatched junctions removes the MARD
regression and the entire benefit together, because only 85 of 7,876 flagged reads have more than
one mismatch. Over 98% of the signal is single-mismatch reads -- the same regime where alignment
noise lives. Requiring the read to also *match* a junction changed nothing (7,861 of 7,876 already
do).

**Lowering the odds fails.** Benefit scales strongly with `novel_odds_per_miss` (CCC +0.00002 at 1.1
to +0.00077 at 2.0) while MARD damage stays nearly constant (+0.00073 to +0.00077). The cost is not
how much mass moves; it is that any mass moves off low-abundance transcripts.

## The lever that works: filter loci, not reads

A locus accumulating one or two disagreeing reads is noise; a genuinely unannotated isoform
accumulates many. `--novel-min-locus-reads` drops under-supported loci before the EM sees them:

| min locus reads | complete: dMARD | 25% holdout: dCCC | 25% holdout: dMARD |
|---|---:|---:|---:|
| 1 | +0.00104 | +0.00077 | +0.00077 |
| 3 | +0.00024 | +0.00077 | -0.00016 |
| **5 (new default)** | **+0.00012** | **+0.00077** | **-0.00029** |
| 10 | +0.00012 | +0.00077 | -0.00035 |

The correlation metrics are completely insensitive to this threshold while the false-positive MARD
damage disappears, which confirms the two were driven by different read populations: benefit by loci
with concentrated evidence, damage by loci with scattered single reads.

## Tuned dose-response

| Holdout | dPearson | dCCC | dRMSE | dMARD |
|---|---:|---:|---:|---:|
| 0% (complete) | -0.00000 | +0.00000 | +0.00 | +0.00012 |
| 10% | +0.00029 | +0.00036 | -0.55 | **-0.00023** |
| 25% | +0.00069 | +0.00077 | -0.72 | **-0.00029** |

All four metrics now improve at both holdout levels, and the model remains neutral on a complete
annotation. This is the first candidate in this study to improve every primary and secondary metric
on truth data without a compensating regression.

Caveat unchanged: this is simulated data with a *random* holdout. Real annotations are incomplete in
structured ways -- whole unannotated loci, systematically truncated 5' ends, unannotated retained
introns -- which a uniform random deletion does not mimic. The magnitude here demonstrates the
mechanism, not the real-world effect size.

---

# Round 20: incompleteness *modes*, real data, and the report

## Simulation design: holdout must target origin transcripts

Following TranSigner, holdout is applied to the **origin transcripts** -- those reads were actually
simulated from -- rather than to the catalogue at large. Deleting unexpressed transcripts orphans no
reads. `make_holdouts.py` builds four modes that stress different mechanisms:

| Mode | What is deleted | Orphaned mass at 25% |
|---|---|---:|
| `origin` | random expressed transcripts | 25.0% |
| `locus` | every isoform of a fraction of expressed genes | 23.8% |
| `minor` | all but the dominant isoform | 3.4% |
| `dominant` | only the dominant isoform | 10.1% |

## The benefit is concentrated, and the modes explain why

| Design | orphaned | flagged | dPearson | dCCC | dRMSE | dMARD |
|---|---:|---:|---:|---:|---:|---:|
| **dominant25** | 10.1% | 12,207 | **+0.00044** | **+0.00058** | **-0.76** | **-0.00038** |
| origin50 | 53.2% | 16,130 | +0.00048 | +0.00056 | +0.04 | -0.00071 |
| origin25 | 25.0% | 12,504 | +0.00002 | +0.00013 | -0.08 | -0.00030 |
| locus25 | 23.8% | 4,521 | +0.00000 | +0.00000 | +0.00 | +0.00007 |
| minor25 | 3.4% | 6,363 | -0.00000 | -0.00000 | +0.00 | +0.00009 |
| minor50 | 8.1% | 8,697 | +0.00000 | +0.00000 | -0.01 | +0.00009 |

* **`dominant` is where it pays.** A highly expressed isoform missing from the annotation dumps its
  reads onto lowly expressed siblings, where the relative distortion is largest. Only 10% of the
  library is orphaned yet the benefit is the largest measured.
* **`locus` yields nothing**, as predicted in Round 16: with every isoform of the gene deleted there
  is no sibling to disagree with, so no junction mismatch is generated. Those reads are lost to
  projection failure instead, which this model does not address.
* **`minor` yields nothing** because the orphaned mass is small and lands on the dominant sibling.
* The earlier +0.00077 from deleting 25% of *all* transcripts was the most favourable case, because
  it also removed unexpressed siblings that would otherwise absorb orphaned reads silently.

Realistic magnitude is therefore around **+0.0006 CCC**, in the specific regime of a missing
dominant isoform, not the +0.0008 first measured.

## Real data

Real HEK293 ONT direct-RNA reads, genome mode, against the same RefSeq annotation:

| Quantity | Simulated (reads from this annotation) | Real reads |
|---|---:|---:|
| projection loss | 0.75% | **18.6%** |
| junction-informative reads disagreeing with every candidate | 0.88% | **3.28%** |

Real data shows 3.7x more unexplained junction evidence than reads simulated from the annotation,
and 25x more projection loss. Both are consistent with genuine annotation incompleteness -- but
**real reads also carry real error**, and ONT direct RNA at 5-10% error will generate spurious
junction mismatches that simulated reads do not. The elevation is therefore an upper bound on
genuine incompleteness, not a measurement of it. Separating the two needs a matched error-model
simulation, which is the obvious next experiment.

## Report, and the gene-label problem

`<prefix>.unexplained.tsv` gives one row per locus:
`locus, gene, n_transcripts, flagged_reads, unexplained_mass, annotated_mass, unexplained_fraction,
transcripts`.

The locus is an **ambiguity component**, not a gene. That is deliberate: bramble exposes transcript
names only, oarfish carries no gene ids in genome mode, and -- more fundamentally -- a genuinely
novel locus has no gene by construction, so keying on genes would make the most interesting cases
unreportable. A gene is emitted when one can be recovered from the reference name (GENCODE's
pipe-delimited header), and `.` otherwise; on RefSeq every row is `.` and the member transcript list
is the identifier. Example from a real run: a 12-transcript locus with 75 flagged reads, 58.1
unexplained against 227.9 annotated -- 20% of that locus's mass unexplained.

## Two bugs found

* **Locus remap:** loci were filtered with `retain` but their surviving *values* still held
  pre-filter ids, so the report indexed position 1040 of a 347-element table. Crashed only on the
  designs where enough loci were dropped.
* **Mid-run rebuild (again):** the first two design runs executed against the previous binary
  because it was rebuilt while they were queued. All designs were rerun on one binary. This is the
  second time in this study; benchmark runs and builds must not overlap.

## Genome *read* mode did not exercise the model

The first real-data run reported junction diagnostics but created no novel states: the raw-read
genome path builds groups directly and never calls `add_projected_group`, so the per-read junction
fields were never populated. Fixed by carrying `(min_junc_misses, max_junc_hits)` per read through
the worker channel. Without this the feature silently did nothing on the most common genome-mode
entry point -- and the counters still printed, which is exactly the kind of partial wiring that
looks like it works.
