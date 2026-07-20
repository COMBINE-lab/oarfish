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
