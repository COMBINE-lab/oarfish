# Coverage-default benchmark expansion: 2026-07-20

## Question and decision

This evaluation asks whether the experimental fixed-abundance blend should
replace the technology-aware `auto` coverage model. It adds independent
truth-bearing simulations, known-concentration SIRV mixtures, a depth series,
a held-out ONT cDNA library, and a full public PacBio SIRV library. Runtime and
peak RSS were measured with `/usr/bin/time -v`; model phase times come from
oarfish metadata.

**Decision: retain `auto` as the recommended coverage mode.** The fixed blend
remains an experimental rank-preserving alternative. It performs very well on
the original development simulation, but the new independent ONT and PacBio
simulations favor `auto` on every accuracy metric. The blend also removes much
of `auto`'s correction on equal-molar SIRVs and is less consistent between the
E0 and E2 mixtures.

The global CLI default is still `none`; changing that separate default is a
release-policy decision. Based on accuracy, documentation should recommend
explicitly running `--coverage-model auto --seq-tech ...` for bulk data.

## New data

- **Independent ONT simulation:** 100,000 reads from
  `ultra_sim_ONT_homo_2M.04.fq.gz`, simulated at 4% error.
- **Independent PacBio-profile simulation:** 100,000 reads from
  `ultra_sim_PB_homo_2M.005.fq.gz`, simulated at 0.5% error.
  Both are from the public [lr-kallisto simulation dataset](https://zenodo.org/records/11201284),
  were generated independently of this coverage work, and encode source
  transcript IDs in read names. They were mapped to the 385,659-sequence
  GENCODE v47 transcriptome. The evaluated subsets contain 81,977 and 82,023
  truth transcripts respectively.
- **ONT direct-RNA SIRV E0/E2:** public SRR6058584 and SRR6058583 reads at
  10k, 25k, 50k, and full depth. Lexogen's supplied workbook provides the 69
  transcript concentrations; E0 is equal-molar and E2 spans multiple levels.
- **Held-out ONT cDNA SIRV E0:** the existing ERR3588903 50k subset.
- **PacBio UHRR + SIRV E0:** all 6,775,127 full-length reads from PacBio's
  [public UHR Iso-Seq dataset](https://downloads.pacbcloud.com/public/dataset/UHR_IsoSeq/).
  Of these, 376,539 mapped to the SIRV transcript reference.

All downloaded data, subsets, outputs, logs, and timing files are under
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/coverage-broader-20260720`.

## Independent simulations

| Technology | Model | Pearson | Spearman | CCC | RMSE | MARD |
|---|---|---:|---:|---:|---:|---:|
| ONT 4% | none | 0.840008 | 0.667849 | 0.829567 | 0.307278 | 0.059144 |
| ONT 4% | **auto** | **0.903328** | **0.710010** | **0.898958** | **0.228732** | **0.037378** |
| ONT 4% | abundance blend | 0.848308 | 0.681482 | 0.838947 | 0.297116 | 0.055271 |
| PacBio 0.5% | none | 0.846556 | 0.676405 | 0.836901 | 0.300457 | 0.054880 |
| PacBio 0.5% | **auto** | **0.872953** | **0.702883** | **0.865794** | **0.268963** | **0.042434** |
| PacBio 0.5% | abundance blend | 0.853485 | 0.689380 | 0.844694 | 0.291913 | 0.051356 |

`auto` wins all five metrics on both held-out simulations. The blend improves
over no coverage, including Spearman, but retains only a minority of the full
correction. This reverses the conclusion from the original development
simulation, where the blend slightly dominated `auto`, and is the strongest
reason not to promote the fixed blend.

## SIRV absolute abundance and depth

For equal-molar E0, correlations and CCC are undefined/uninformative because
truth has zero variance. RMSE, MARD, and the coefficient of variation (CV) of
normalized estimates measure deviation from the expected uniform profile;
lower is better.

| Dataset | Model | RMSE / estimate CV | MARD |
|---|---|---:|---:|
| dRNA E0 10k | none / auto / blend | 0.6258 / **0.6050** / 0.6195 | 0.2570 / **0.2463** / 0.2545 |
| dRNA E0 25k | none / auto / blend | 0.6311 / **0.5959** / 0.6239 | 0.2555 / **0.2385** / 0.2519 |
| dRNA E0 50k | none / auto / blend | 0.6355 / **0.5776** / 0.6253 | 0.2516 / **0.2243** / 0.2469 |
| dRNA E0 full (184k) | none / auto / blend | 0.6292 / **0.5436** / 0.6093 | 0.2513 / **0.2111** / 0.2420 |
| cDNA E0 50k | none / auto / blend | 1.5389 / **1.5350** / 1.5367 | 0.5039 / **0.4954** / 0.5006 |
| PacBio E0 full | none / auto / blend | 1.2073 / **1.2068** / 1.2104 | 0.4231 / **0.4218** / 0.4246 |

The direct-RNA correction strengthens smoothly with depth. The blend moves
toward `auto` as abundance rises but consistently gives back most of the E0
improvement. On full PacBio, `auto` is a small improvement while the blend is
slightly worse than no coverage.

For the non-uniform E2 direct-RNA mixture:

| Depth | Model | Pearson | Spearman | CCC | RMSE | MARD |
|---:|---|---:|---:|---:|---:|---:|
| 10k | none | 0.891443 | 0.771141 | 0.889405 | 0.687783 | 0.318297 |
| 10k | auto | 0.892835 | 0.770420 | 0.890410 | 0.686896 | 0.317586 |
| 10k | blend | **0.892853** | **0.774245** | **0.890432** | **0.686800** | **0.317473** |
| 50k | none | 0.877134 | 0.755903 | 0.875127 | 0.730845 | 0.278627 |
| 50k | **auto** | **0.879573** | 0.755903 | **0.876925** | **0.729458** | 0.274831 |
| 50k | blend | 0.879555 | **0.758371** | 0.876907 | 0.729515 | **0.274794** |
| full | none | 0.873441 | 0.762610 | 0.871741 | 0.738615 | 0.270990 |
| full | **auto** | **0.877284** | **0.763159** | **0.874771** | **0.735058** | **0.263920** |
| full | blend | 0.877275 | **0.763159** | 0.874761 | 0.735095 | 0.263947 |

Here the blend and `auto` are nearly identical at moderate/full depth. The
blend's benefit is therefore sample-dependent: it suppresses useful correction
on uniform E0 but barely changes the high-abundance E2 result.

## E0-to-E2 fold-change recovery

Mixture estimates were normalized within sample before computing transcript
log2 fold changes. At 25k, `auto` has the best Pearson, Spearman, RMSE, and
direction accuracy. At 10k and 50k, no one model wins every metric, but the
blend has substantially worse fold-change RMSE.

| Depth | Model | FC Pearson | FC Spearman | FC RMSE | Direction accuracy |
|---:|---|---:|---:|---:|---:|
| 10k | none | **0.8486** | **0.8698** | **1.4059** | 0.8696 |
| 10k | auto | 0.8362 | 0.8647 | 1.4723 | **0.8986** |
| 10k | blend | 0.7724 | 0.8660 | 2.1210 | 0.8696 |
| 25k | none | 0.8726 | 0.8697 | 1.2586 | 0.8841 |
| 25k | **auto** | **0.8756** | **0.8882** | **1.2429** | **0.9275** |
| 25k | blend | 0.7989 | 0.8707 | 1.8414 | 0.8841 |
| 50k | none | **0.8791** | **0.8636** | **1.2358** | 0.9130 |
| 50k | auto | 0.8634 | 0.8593 | 1.2707 | **0.9420** |
| 50k | blend | 0.7804 | 0.8627 | 1.8261 | 0.8986 |

The likely cause is that a fixed abundance gate applies different correction
strengths to the same transcript in samples with different abundance
distributions. This is undesirable for differential analysis and rules out the
current blend as a general default.

## Performance cost

| Dataset | Model | Wall s | Coverage s | Peak RSS |
|---|---|---:|---:|---:|
| independent ONT 100k | none | 107.68 | ~0 | 3.58 GiB |
| independent ONT 100k | auto | 107.43 | 0.103 | 3.56 GiB |
| independent ONT 100k | blend | 107.88 | 1.279 | 3.66 GiB |
| independent PacBio 100k | none | 102.03 | ~0 | 2.78 GiB |
| independent PacBio 100k | auto | 103.45 | 0.102 | 2.82 GiB |
| independent PacBio 100k | blend | 106.21 | 1.352 | 2.86 GiB |
| full dRNA E2 | none | 42.21 | ~0 | 132.8 MiB |
| full dRNA E2 | auto | 42.42 | 0.287 | 131.5 MiB |
| full dRNA E2 | blend | 41.40 | 0.528 | 120.5 MiB |
| full PacBio E0 | none | 132.43 | ~0 | 138.0 MiB |
| full PacBio E0 | auto | 132.56 | 0.089 | 127.5 MiB |
| full PacBio E0 | blend | 132.97 | 0.338 | 139.7 MiB |

Mapping dominates end-to-end runtime. `auto` adds about 0.1--0.3 seconds on
these panels. The blend adds an additional warmup EM pass: about 1.2 seconds on
the large-reference simulations and 0.25--0.35 seconds on the spike-ins. Peak
RSS differences are small relative to mapping and are not consistently signed.

## Recommendation

1. Recommend `--coverage-model auto --seq-tech <technology>` for routine bulk
   quantification.
2. Keep `abundance-blend` experimental and do not expose it as the default.
   It can still be useful for analyses explicitly prioritizing rank recovery
   on data resembling the original development simulation.
3. Keep `none` as an accessible conservative control. Whether the CLI default
   should move from `none` to `auto` should be handled as an explicit release
   change, with migration notes.
4. If blend work continues, optimize a *sample-consistent* objective including
   fold-change preservation. A fixed per-sample abundance threshold is not a
   sufficient gate.

## Reproduction

Repository inputs and drivers:

- `benchmarks/coverage_default_public_manifest.tsv`
- `scripts/extract_sirv_concentrations.py`
- `scripts/run_coverage_default_benchmark.py`
- `scripts/evaluate_sirv_mixture.py`

The public simulation files have MD5 values supplied by Zenodo
(`d8f87a9517ed1d4816ec8fac370480a4` and
`f311fb893418a37366a2b5197077a666`).
