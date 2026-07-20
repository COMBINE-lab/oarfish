# Adaptive coverage evaluation: 2026-07-19

## Model and safeguards

The experimental `adaptive` model retains the selected hybrid's logistic and
endpoint signals but makes the endpoint contribution data-adaptive:

1. assign reads deterministically to five folds;
2. train each read's endpoint distribution without that fold's unambiguous
   reads;
3. smooth the 20 by 20 endpoint table locally with a separable 1-2-1 kernel;
4. select the symmetric Dirichlet prior mass from
   `10, 30, 100, 300, 1000, 3000, 10000` by held-out predictive likelihood;
5. gate the combined evidence by both endpoint support and agreement between
   the normalized logistic and endpoint distributions; and
6. shrink the result toward uniform and cap its per-read Bayes factor (default
   4) before EM inference.

These choices address leakage, sparse endpoint cells, sample-specific model
disagreement, and catastrophic overcorrection. The overall default remains
`--coverage-model none`.

## Benchmarks

All comparisons use the same mapper/filter settings within a dataset and
rescale estimates to the known truth library size. Metrics are computed by
`scripts/evaluate_quant.py`. The public data and generated subsets live outside
the repository under
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data`.

- **Synthetic dRNA:** fixed 1,356,482-read BAM with transcript truth. This was
  used during model development and is not held out.
- **ONT direct RNA SIRV E0:** the first 50,000 reads of SRR6058584, plus fixed
  first-10,000 and first-25,000 prefixes for depth sensitivity.
- **Held-out ONT cDNA SIRV E0:** the first 50,000 reads of ERR3588903. This
  protocol/dataset was selected before inspecting adaptive results. The source
  FASTQ MD5 is `995bb92361e88e7d77f8e6fc631503ba`.

### Accuracy at full evaluated depth

| Dataset | Model | CCC | MARD | RMSE |
|---|---|---:|---:|---:|
| synthetic dRNA | none | 0.956952 | 0.169178 | 89.4590 |
| synthetic dRNA | logistic | 0.994990 | **0.142519** | 30.7349 |
| synthetic dRNA | endpoint | 0.955198 | 0.153293 | 91.4660 |
| synthetic dRNA | hybrid | 0.996132 | 0.143295 | 26.9672 |
| synthetic dRNA | adaptive | **0.996935** | 0.143259 | **24.0550** |
| ONT dRNA SIRV E0 50k | none | 0.639192 | 0.166302 | 0.529723 |
| ONT dRNA SIRV E0 50k | logistic | 0.632240 | 0.166512 | 0.537734 |
| ONT dRNA SIRV E0 50k | endpoint | **0.732652** | **0.141257** | **0.425909** |
| ONT dRNA SIRV E0 50k | hybrid | 0.702601 | 0.144598 | 0.458714 |
| ONT dRNA SIRV E0 50k | adaptive | 0.671348 | 0.154340 | 0.493311 |
| held-out ONT cDNA E0 50k | none | 0.273643 | 0.268102 | 1.148663 |
| held-out ONT cDNA E0 50k | logistic | 0.273379 | 0.266203 | 1.149426 |
| held-out ONT cDNA E0 50k | endpoint | **0.277680** | **0.257236** | **1.137111** |
| held-out ONT cDNA E0 50k | hybrid | 0.275404 | 0.259503 | 1.143595 |
| held-out ONT cDNA E0 50k | adaptive | 0.274599 | 0.264218 | 1.145909 |

Adaptive is the best synthetic model by CCC and RMSE and improves all three
reported metrics over no coverage on both public datasets. It intentionally
captures less of the endpoint-only gain on SIRV: its safeguards favor a small,
transferable correction over a strong correction when model agreement or
support is limited. Endpoint remains strongest on these two SIRV E0 samples,
but that observation is not yet evidence that its unbounded behavior
generalizes beyond SIRV.

### Direct-RNA depth sensitivity

| Reads | Model | CCC | MARD | RMSE |
|---:|---|---:|---:|---:|
| 10,000 | none | 0.644593 | 0.167276 | 0.523536 |
| 10,000 | endpoint | 0.718655 | 0.144364 | 0.441149 |
| 10,000 | hybrid | 0.687806 | 0.149402 | 0.475013 |
| 10,000 | adaptive | 0.656112 | 0.162923 | 0.510442 |
| 25,000 | none | 0.641661 | 0.166783 | 0.526890 |
| 25,000 | endpoint | 0.728911 | 0.141126 | 0.429977 |
| 25,000 | hybrid | 0.696602 | 0.145511 | 0.465308 |
| 25,000 | adaptive | 0.661218 | 0.159550 | 0.504678 |
| 50,000 | none | 0.639192 | 0.166302 | 0.529723 |
| 50,000 | endpoint | 0.732652 | 0.141257 | 0.425909 |
| 50,000 | hybrid | 0.702601 | 0.144598 | 0.458714 |
| 50,000 | adaptive | 0.671348 | 0.154340 | 0.493311 |

The adaptive correction grows smoothly with evidence: mean reliability is
0.402 at 10k, 0.520 at 25k, and 0.604 at 50k. It improves on no coverage at
every depth without the abrupt low-support behavior that an ungated empirical
table can produce.

## Runtime and learned diagnostics

Coverage time is wall-clock time around model construction/application; EM is
timed separately in `.meta_info.json`.

| Dataset | Model | Coverage time | EM time | Mean reliability | BF-capped reads |
|---|---|---:|---:|---:|---:|
| dRNA SIRV 10k | none | <0.001 ms | 30.0 ms | - | - |
| dRNA SIRV 10k | hybrid | 1.63 ms | 26.9 ms | - | - |
| dRNA SIRV 10k | adaptive | 2.07 ms | 27.9 ms | 0.402 | 1 |
| dRNA SIRV 25k | none | <0.001 ms | 72.7 ms | - | - |
| dRNA SIRV 25k | hybrid | 3.55 ms | 65.9 ms | - | - |
| dRNA SIRV 25k | adaptive | 4.86 ms | 66.1 ms | 0.520 | 418 |
| dRNA SIRV 50k | adaptive | 9.81 ms | 146.2 ms | 0.604 | 2,180 |
| held-out cDNA 50k | none | <0.001 ms | 80.2 ms | - | - |
| held-out cDNA 50k | hybrid | 11.26 ms | 70.7 ms | - | - |
| held-out cDNA 50k | adaptive | 14.79 ms | 73.0 ms | 0.731 | 17,051 |
| synthetic dRNA 1.36m | adaptive | 404.39 ms | 11.34 s | 0.957 | 26,469 |

The adaptive overhead is linear and small relative to mapping and EM: 15 ms
versus 9.85 s of mapping on held-out cDNA, and 0.40 s versus 11.34 s of EM on
the large fixed-alignment synthetic benchmark. All three evaluated datasets
selected prior mass 10, the least-regularized edge of the frozen candidate
grid. This is a useful warning: future work should test a broader prior grid on
new training data rather than retuning it on these evaluation sets.

## Decision and next evaluation

- Keep `adaptive` experimental and opt-in. Its cross-dataset direction is
  consistently safe, but three mixtures and two ONT library types are not
  sufficient to claim universal parameter transfer.
- Retain the frozen defaults: five folds, support scale 25, logistic/endpoint
  weights 1.0/0.5, and Bayes-factor cap 4.
- Before considering a default change, evaluate biological-replicate
  consistency and differential-expression calibration, PacBio Iso-Seq, a
  non-SIRV mixture with transcript-level truth, and degradation/RNA-integrity
  strata. Report per-gene isoform accuracy in addition to global metrics.
- Learn any next-stage hyperparameters with nested sample-level validation:
  whole samples, not reads, must be the outer holdout unit.
