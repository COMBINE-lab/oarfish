# Degradation-kernel challenger study: 2026-07-20

This study tests whether flexible degradation distributions improve on the
frozen three-component competing-risk model. The incumbent uses a constant
per-base degradation hazard. Challengers were deliberately bounded:

- `piecewise2`: constant plus two regularized directional shapes;
- `piecewise3`: constant plus four regularized directional/central shapes;
- `beta`: length-dependent mean retained fraction and concentration 2, 5, or 10.

Piecewise models embed the constant model and penalize squared log-rate
departures. All models retain the same intact and length-invariant technical
components, cross-fitting, degradation posterior, evidence gate, reliability
gate, and Bayes-factor cap. No model was tuned after examining abundance truth.

## Endpoint fit and stability

Across all 15 PRJEB53210 degradation libraries:

| Kernel | Held-out log density | Mean hazard fold SD | Selected shapes |
|---|---:|---:|---|
| constant | -2.22634 | 0 | constant |
| piecewise2 | -2.21159 | 0 | 0.67, 1.33 |
| piecewise3 | -2.20325 | 0 | 0.67, 1.00, 1.33 |
| beta | **-2.11577** | 0.0098 | concentration 2 or 5 |

All flexible models fit held-out endpoints better, and piecewise hazards are
stable across folds. This does not translate monotonically into quantification
accuracy: beta is the clearest counterexample.

## LongBench

Values are means across eight 50,000-read ONT direct-RNA libraries.

| Kernel | CCC | Pearson | RMSE | MARD | Coverage time |
|---|---:|---:|---:|---:|---:|
| **constant** | **0.497474** | **0.712606** | **6257.00** | **0.449311** | **0.100 s** |
| piecewise2 | 0.497413 | 0.712542 | 6257.46 | 0.449320 | 0.203 s |
| piecewise3 | 0.497388 | 0.712519 | 6257.73 | 0.449321 | 0.313 s |
| beta | 0.496288 | 0.711223 | 6267.30 | 0.449427 | 2.104 s |

Piecewise regressions are very small but consistent in the aggregate. Beta is
both less accurate and roughly 21 times slower than constant at this depth.

At 250,000 HCC827 reads, piecewise2 is stable and slightly better than constant:

| Kernel | CCC | Pearson | RMSE | Coverage time |
|---|---:|---:|---:|---:|
| constant | 0.562581 | 0.739437 | 6869.19 | 0.426 s |
| piecewise2 | **0.562617** | **0.739453** | **6868.54** | 0.919 s |

The gain is only 0.000035 CCC and does not offset the complete-panel result or
the twofold runtime cost.

## Degradation trajectories

| Kernel | TS10 Pearson | TS12 7.2 Pearson | TS12 7.3 Pearson |
|---|---:|---:|---:|
| constant | 0.96500 | **0.94876** | **0.94013** |
| piecewise2 | **0.96681** | 0.94833 | 0.93930 |
| piecewise3 | 0.96671 | 0.94835 | 0.93930 |
| beta | 0.96439 | 0.94831 | 0.93885 |

Piecewise hazards improve TS10 but regress both TS12 endpoints. There is no
uniform robustness winner beyond the incumbent.

## Frozen truth and large-sample runtime

| Kernel | Synthetic CCC | Synthetic RMSE | SIRV CCC | Source CCC | Synthetic coverage time |
|---|---:|---:|---:|---:|---:|
| **constant** | **0.996926** | **24.092** | 0.671305 | **0.983387** | **2.77 s** |
| piecewise2 | 0.996907 | 24.167 | **0.671326** | **0.983387** | 7.55 s |
| piecewise3 | 0.996907 | 24.166 | 0.671255 | **0.983387** | 12.69 s |
| beta | 0.994516 | 32.118 | 0.640705 | 0.983387 | 63.97 s |

The piecewise differences are small, but neither improves the large synthetic
truth. Beta decisively fails truth and runtime gates despite its much better
held-out endpoint density.

## Parameter-recovery checks

Deterministic simulations verify that the regularized piecewise2 challenger
collapses to constant hazard for constant data and recovers the direction of a
simulated 0.67/1.33 regional hazard. Thus its negative benchmark result is not
caused by a nonfunctional implementation.

## Decision

Keep `constant` as the default and validated kernel. Retain only the bounded
piecewise2 challenger behind the explicitly experimental
`--degradation-kernel constant|piecewise2` option. The piecewise3 and beta code
was removed after evaluation to avoid carrying unjustified complexity. Do not
automatically select a kernel from endpoint likelihood: the beta experiment
proved that this would select a worse quantification model. Piecewise2 is the
only credible future challenger, but its current gains are dataset-specific and
too small to justify replacement.

Artifacts are under `oarfish-evaluation-data/degradation/kernel-*` and
`oarfish-evaluation-data/longbench/kernel-*`.
