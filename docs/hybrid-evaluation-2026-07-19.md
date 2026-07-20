# Hybrid coverage evaluation: 2026-07-19

The hybrid model combines normalized coverage likelihoods in log space:

`log p = logistic_weight * log(p_logistic) + endpoint_weight * gate * log(p_endpoint)`

The per-read support gate is `s / (s + endpoint_support_scale)`, where `s` is
the mean number of unambiguous training observations in the candidates'
endpoint cells. This makes unsupported endpoint evidence neutral. Combination
uses a probability floor, stable softmax, and per-read normalization.

## Weight sweep

All runs use support scale 25. The large synthetic runs reuse the fixed BAM and
truth from the preceding coverage report. SIRV E0 runs use the identical first
50,000 reads of SRR6058584.

| Logistic weight | Endpoint weight | Synthetic CCC | Synthetic MARD | Synthetic RMSE | SIRV CCC | SIRV MARD | SIRV RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.50 | 0.986983 | 0.143080 | 49.2501 | 0.704631 | **0.144445** | **0.456486** |
| 0.75 | 0.25 | 0.994548 | **0.142315** | 31.9998 | 0.678820 | 0.150798 | 0.484979 |
| 1.00 | 0.25 | 0.995937 | 0.142579 | 27.6751 | 0.680976 | 0.149686 | 0.482583 |
| **1.00** | **0.50** | **0.996132** | 0.143295 | **26.9672** | 0.702601 | 0.144598 | 0.458714 |

For context, standalone logistic has synthetic CCC 0.994990 and RMSE 30.7349;
standalone endpoint has SIRV CCC 0.732652 and RMSE 0.425909. The selected
`1.0/0.5` hybrid improves synthetic CCC and RMSE beyond logistic while retaining
68% of endpoint's CCC and RMSE improvement over no coverage on SIRV E0. Its
synthetic MARD is 0.5% worse than standalone logistic but remains 15.3% better
than no coverage.

On the 1,634-read source-labelled direct-RNA set, the selected hybrid gives CCC
0.986243, MARD 0.0000871, RMSE 0.010849, and Spearman 0.993215. This is nearly
identical to standalone endpoint and better than no coverage on every metric.

## Decision

- Add `hybrid` as an experimental opt-in coverage model.
- Use logistic weight 1.0, endpoint weight 0.5, and support scale 25 as its
  defaults. Keep the overall oarfish default at `--coverage-model none`.
- Do not tune more parameters on these same datasets. The next test must use a
  held-out public mixture or technology, ideally after obtaining verified SIRV
  E2 transcript concentrations.

Artifacts are under
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/hybrid-sweep`.
