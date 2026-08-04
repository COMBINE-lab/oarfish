# Inference and coverage evaluation

Every inference change is evaluated against fixed alignments before mapping is
included.  A stage report records the oarfish commit, reference and input
checksums, command, thread count, random seed, elapsed time, peak RSS, EM
evaluations, convergence status, and the JSON emitted by
`scripts/evaluate_quant.py`.

The evaluator rescales estimated counts to the truth library size before
count-scale metrics are calculated; use `--raw-counts` only when the two files
represent the same library.  Use `--strip-version` when the reference uses
versioned accessions but the truth does not.

## Required comparison matrix

For EM acceleration, run `none`, `squarem`, and `daarem` with coverage disabled,
then repeat with the accepted coverage model and with bootstrap replicates.  Use
a tightly converged `none` run as the numerical reference.  Report count-mass
error, maximum relative error for transcripts with at least one expected read,
M-step evaluations, wall time, and peak RSS.

For coverage, always include `none`, repaired `logistic`, and `endpoint`.  Hold
the BAM and all filtering options fixed.  Primary metrics are MARD, CCC, and
within-gene isoform-fraction error; secondary metrics are RMSE, Pearson,
Spearman, expressed-transcript precision/recall, wall time, and peak RSS.

## Truth tier versus comparator tier

Panel results must be reported in two tiers, split by the provenance of the truth column, and
never pooled into a single mean:

* **Truth tier** (`truth_type = counts`): molecular ground truth or known spike-in concentration —
  the synthetic dRNA set, the independent ONT/PacBio simulations, and the SIRV mixtures.
* **Comparator tier** (`truth_type = illumina`): matched short-read abundance. This is a
  *different technology with its own assignment behaviour*, not ground truth.

The reason is empirical. `benchmarks/coverage_ablation_manifest.tsv` is 27 comparator libraries
and one truth-bearing sample, so a pooled mean is ~96% comparator — and on three separate
occasions a candidate that improved the comparator was contradicted by molecular truth:

| Round | Comparator | Molecular truth |
|---|---|---|
| Unique-read-only profiles (2026-07-21) | better Spearman/MARD in all 24 libraries | synthetic CCC 0.99704 -> 0.94994 |
| Endpoint nucleotide measure (2026-07-24) | ONT dRNA Spearman +0.000515 (8/9) | true-transcript discrimination -0.043 |
| Score temperature (2026-07-24) | prefers `D -> 0.5`, monotonically | prefers `D -> 20`, monotonically |

**Promotion additionally requires that no primary metric regress in the truth tier.** A gain that
exists only in the comparator tier is evidence about agreement with Illumina, not about accuracy.

Report both tiers with:

```
python3 scripts/summarize_panel.py <results.tsv ...> \
    --manifest benchmarks/coverage_ablation_manifest.tsv \
    --manifest benchmarks/truth_tier_manifest.tsv \
    --control full --candidate <arm>
```

It prints per-tier means and win counts and emits a verdict, including `INDETERMINATE` when a run
contains no truth-bearing sample. `benchmarks/truth_tier_manifest.tsv` holds the truth-bearing
panel; `scripts/run_coverage_ablation.py` accepts either a `bam` column (alignment mode) or
`reads` + `reference` (raw-read mode), so both tiers run under one driver.

### Determinism, and a Spearman noise floor in raw-read mode

Alignment (BAM) mode is fully deterministic: repeated identical runs are byte-identical, so every
metric on the 28-case comparator panel and on the BAM rows of the truth tier is exactly
reproducible.

**Raw-read mode is not.** Parallel mapping breaks ties differently between runs. The effect on
abundance is negligible — on `independent-pb-100k`, 59 of 385,659 transcripts differ and total
movement rounds to 0.0 reads — but Spearman is rank-based and hypersensitive to reordering among
near-zero transcripts. Measured across two identical `--threads 8` runs:

| Metric | run 1 | run 2 | \|delta\| |
|---|---:|---:|---:|
| pearson | 0.962394 | 0.962394 | 0.000000 |
| ccc | 0.961812 | 0.961812 | 0.000000 |
| rmse | 0.136305 | 0.136306 | 0.000000 |
| mard | 0.010690 | 0.010690 | 0.000000 |
| **spearman** | 0.753194 | 0.748438 | **0.004756** |

So on the five raw-read rows of the truth tier (`independent-*`, `sirv-*`), **Spearman deltas below
about 0.005 are noise**, which is larger than nearly every effect this project measures. Judge
those rows on CCC, Pearson, RMSE and MARD; treat their Spearman as uninformative unless it moves by
more than the floor, or re-run with repeats. The BAM rows (`synthetic-drna`, `kinnex-*`) have no
such caveat.

Two cautions:

* The truth tier is itself heterogeneous. Alignment score alone picks the true transcript for
  98.7% of ambiguous reads on the independent ONT simulation but only 43.7% on the synthetic dRNA
  set, so the two disagree about how much work coverage has to do. Do not treat any single case as
  decisive.
* Aggregate abundance accuracy and per-read assignment accuracy are distinct axes and can move in
  opposite directions on the *same* dataset — the nucleotide-measure candidate improved synthetic
  CCC while lowering true-transcript discrimination on that same sample. When a candidate changes
  a length-dependent term, check both.

## Stage report template

1. **Hypothesis:** the new information the model is expected to capture.
2. **Data:** manifest, checksums, truth source, technology, and whether the data
   are tuning or held out.
3. **Commands:** complete invocations and environment information.
4. **Results:** machine-readable metrics plus stratification by abundance,
   transcript length, multimapping degree, isoform similarity, and completeness.
5. **Failures:** components with the largest positive and negative changes.
6. **Decision:** promote, revise, or reject.  Promotion requires improved median
   primary metrics, improvement in two technology classes, no held-out or
   full-length-control regression above 2%, and acceptable resource overhead.

## Implemented checkpoints

### Shared EM driver and acceleration

The serial, parallel, and bootstrap paths now share an absolute-relative-change
criterion and a strict M-step budget.  SQUAREM and DAAREM have synthetic
fixed-point tests for convergence, feasibility, and budget compliance.  Public
accuracy and resource measurements remain required before changing the default.

### Repaired logistic baseline

Coverage uses zero-based half-open intervals internally, includes the terminal
bin, applies posterior/interval weights, and no longer mutates raw bins when
adding pseudocounts.  Boundary and weighted-terminal-bin behavior is unit tested.

### Joint endpoint model

The first novel candidate learns a regularized joint distribution of normalized
5-prime and 3-prime gaps in four transcript-length strata.  Only single-candidate
reads train it.  A uniform prior equivalent to 1,000 observations makes sparse
strata neutral.  Endpoint-cell separation and empty-training neutrality are unit
tested.  The mixture, predictive, and segment stages must not be promoted until
this model has a completed public stage report.

The first implementation checkpoint, including failures and the public SIRV
manifest, is recorded in
[evaluation-results-2026-07-18.md](evaluation-results-2026-07-18.md).

The first truth-based comparison of all three coverage models is recorded in
[coverage-evaluation-2026-07-19.md](coverage-evaluation-2026-07-19.md).

The subsequent support-gated hybrid implementation and weight sweep are
recorded in [hybrid-evaluation-2026-07-19.md](hybrid-evaluation-2026-07-19.md).
