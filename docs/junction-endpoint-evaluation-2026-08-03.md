# Splice-junction endpoint evidence: evaluation (2026-08-03)

## Idea

An alignment whose transcript-coordinate boundary falls at an *internal* splice
junction of the transcript it is aligned to is suspicious: a molecule generated
by that transcript has no reason to stop at an interior exon boundary. Such an
endpoint is usually better explained by a sibling isoform that *terminates*
there. Down-weighting those alignments should improve isoform disambiguation.

## Verdict

**The signal is real, large, and reproducible. It does not convert into
quantification accuracy.** Net effect at the best setting is **+0.00026 Spearman**
on exact truth, for +5.5% wall time and a required GTF. Recommendation: keep
opt-in (the current default is off); do not enable by default.

## What was measured

### 1. The raw signal is unambiguous

NanoSim read names carry the true source transcript, so the likelihood ratio is
directly measurable. All 34,899,437 alignments of NanoSim NA12878 1D-cDNA
(9.15 M true-source, 24.3 M decoy):

| window | 3' end LR (decoy/true) | 5' start LR |
|---|---|---|
| exact (0 nt) | **9.90x** | **7.14x** |
| 2 nt | 4.49x | 5.30x |
| 5 nt | 3.13x | 3.42x |
| 10 nt | 2.48x | 2.34x |
| 20 nt | 2.05x | 1.67x |

An alignment ending exactly on an internal junction is ~10x more likely to
belong to a decoy than to the true source. Independently implemented in Python
and Rust with agreeing results.

### 2. It moves posterior mass to the correct transcript

Applying `exp(-lambda)` to junction-terminated candidates in oarfish's own
per-read posteriors (`--write-assignment-probs`) and renormalising within each
read. Parameters fit on NA12878-cDNA (w=5, lambda=6), then held fixed:

| sample | baseline mass on true source | delta | per movable read |
|---|---|---|---|
| na12878-cdna *(fit)* | 0.7363 | +0.00725 | +0.0247 |
| na12878-drna *(held out)* | 0.7987 | +0.00651 | +0.0289 |
| h9-cdna *(held out)* | 0.7560 | +0.01076 | +0.0359 |
| h9-drna *(held out)* | 0.8340 | +0.00602 | +0.0335 |

No overfitting: held-out samples match or beat the fitted one.

### 3. It does not convert into accuracy

Implemented in oarfish and run end to end. Panel B, exact truth, w=5,
lambda=0.5:

| sample | logistic | + junction | delta Spearman | delta MARD |
|---|---|---|---|---|
| na12878-cdna | 0.886717 | 0.887282 | +0.000565 | -0.000217 |
| na12878-drna | 0.919238 | 0.919404 | +0.000165 | -0.000107 |
| h9-cdna | 0.908404 | 0.908670 | +0.000265 | -0.000121 |
| h9-drna | 0.936348 | 0.936513 | +0.000166 | -0.000023 |
| tksm-RSII | 0.935371 | 0.935606 | +0.000235 | -0.000114 |
| tksm-SQ2 | 0.945946 | 0.946093 | +0.000147 | -0.000107 |
| **mean** | | | **+0.000257** | **-0.000115** |

Experimental data (LongBench, 28 real ONT/PacBio cell-line samples, matched
Illumina comparator, GENCODE v47): mean **+0.000067**, 20 wins / 7 losses.
By technology: ont-cdna +0.000078 (7/8), ont-drna +0.000088 (6/9),
pac-bio +0.000041 (7/11). The `synthetic-drna` sample uses a different
reference the GTF does not cover and shows exactly 0.000000 delta, confirming
the feature is inert when the annotation does not apply.

Real data therefore shows a *smaller* effect than simulation, contradicting the
expectation that greater real-world isoform complexity would amplify it.

### 4. Cost

| | wall | peak RSS |
|---|---|---|
| logistic | 1:01.99 | 2.87 GB |
| + junction endpoint | 1:05.40 | 2.88 GB |

+5.5% wall time (dominated by the GTF scan), no meaningful memory change, and a
GTF becomes a required input.

## Why a 10x signal buys 0.0003 Spearman

Three measured reasons:

1. **Low volume.** Only 14-21% of alignments are junction-terminated, and the
   sharp high-LR core (exact coincidence) is ~1% of decoy alignments.
2. **It fires where the EM is already nearly right.** Moving a read from 0.74 to
   0.76 posterior on the true source rarely changes a transcript's rank.
3. **It penalises true transcripts too.** Genuinely expressed transcripts carry
   a mean junction-terminated fraction of 0.12 (versus 0.18 for absent ones),
   so the term removes real evidence as well as spurious evidence. The
   per-transcript separation is only 1.54x, far weaker than the 9.9x
   per-alignment ratio.

This is the same pattern as the length-bias result in
`docs/coverage-auto-rebenchmark-results-2026-08-03.md`: a large, real,
well-characterised effect that the metric of record is nearly insensitive to.

## On the pre-registered criterion

Before seeing results the criterion was: *a gain on >=5 of 6 Panel B samples, no
dataset family regressing more than 0.005, and it must beat `logistic`.*

**As literally worded, the feature passes**: 6/6 samples gain, nothing regresses,
MARD improves on all six, and the comparison is same-binary flag-on/flag-off on
identical input, so it is deterministic rather than noise. The sign is certainly
real.

But the pre-registration **omitted a minimum effect size**, and the same note
projected "~+0.009 if it captures even a third of the detection headroom". The
observed +0.00026 is roughly 35x smaller than that projection. The criterion was
therefore satisfied on consistency while failing the implied magnitude by a wide
margin. That gap is a defect in how the criterion was written, and is recorded
here rather than resolved by moving the threshold after the fact.

## Methodological finding: proxy ladders over-predict

Each cheaper stage over-predicted the next:

| stage | signal |
|---|---|
| likelihood ratio | 9.90x |
| posterior mass moved per movable read | +0.0247 |
| offline EM (degraded candidate set) | +0.0021 Spearman |
| **real oarfish EM** | **+0.00026 Spearman** |

Two specific traps:

- The reweight test optimum (**lambda=6**) was *catastrophic* in the real EM
  (**-0.0167** Spearman). One-shot per-read renormalisation does not model the
  EM's iterative compounding; the real optimum is lambda~0.5, twelve times
  milder.
- The offline EM rescoring reconstructed conditionals as `p_rt / theta_t` from a
  posterior dump, which silently drops every candidate on a transcript the EM
  zeroed -- 58% of candidates here. Its baseline (0.8751) sits well below the
  real system (0.8867), so its deltas are not comparable to a benchmark.

## Can annotation help at all? An information ceiling

The question of whether *this particular formulation* failed, versus whether
annotation-derived evidence is fundamentally unable to help, is decidable. For
each read the EM misassigns (true source holding <50% of the posterior), compare
the alignment score of the **true source** against that of the transcript the EM
actually favoured:

| | NA12878-cdna | H9-cdna | NA12878-drna |
|---|---|---|---|
| misassigned reads | 1,813,676 | 2,096,794 | 1,080,645 |
| **TIED alignment score** | **77.04%** | **74.42%** | **66.77%** |
| true source scores HIGHER | 22.56% | 25.28% | 32.89% |
| true source scores LOWER | 0.38% | 0.28% | 0.31% |
| junction signal favours truth, within the TIED set | **0.30%** | **0.20%** | **0.15%** |

**About three quarters of all misassignment is information-theoretically
irreducible from the alignment.** The true source and the winning decoy score
*identically*. An identical score over the read's span means the two transcript
sequences are indistinguishable there, which in turn means their exon structure
over that span is identical -- so no feature derived from the annotation and the
alignment can separate them, regardless of formulation. Measured directly: the
junction-endpoint signal favours the truth in 0.15-0.30% of misassigned reads
within that bucket.

The remaining ~27% are cases where the true source **already scores higher** and
is nonetheless outvoted by the abundance prior. There the alignment evidence
already points the right way, so an annotation feature is largely re-deriving
information the model has and is choosing to override; consistent with this, the
junction signal flags only 6.5-9.7% of misassignments in that bucket.

This is the answer to "is it the idea or the information?": **it is the
information.** Splice-pattern compatibility, junction-aware effective length, or
any other annotation-derived per-alignment term faces the same ceiling, because
they are all functions of the same alignment over the same sequence.

What the decomposition does suggest, and it is not annotation-based: the ~27%
bucket is a *calibration* problem between alignment likelihood and abundance
prior, not a missing-feature problem. oarfish already exposes
`--score-prob-denom` for exactly that trade-off. That is a scalar recalibration
of existing machinery rather than new surface area, and it addresses roughly a
quarter of the error rather than a fraction of a percent.

## Implementation

Opt-in and isolated, so removal is mechanical. `grep -rn "junction-endpoint" src/`
finds every site:

- `src/util/junction_endpoint.rs` -- all logic (GTF parse, transcript-coordinate
  junction offsets, per-alignment penalty)
- `src/util.rs` -- one `pub mod` line
- `src/prog_opts.rs` -- three CLI arguments, in one tagged block
- `src/bulk.rs` -- one tagged block applying the penalty to
  `coverage_probabilities` immediately before the EM

Verified no-op when the flag is absent: Spearman 0.886717 and MARD 0.058173
reproduce to six decimals, maximum per-transcript count difference 3.3e-10
(thread-scheduling noise only).

Annotation reconciliation is reported at INFO and in `meta_info.json`. GENCODE
v47 reconciles perfectly against its own transcript FASTA (0 length
mismatches, 356,707 transcripts with internal junctions); RefSeq
GCF_000001405.40 has 4,925 length mismatches (2.8%) that are excluded from the
term rather than risk mis-registered coordinates, plus 5,378 single-exon
transcripts for which it is a no-op.

## Artifacts

- `scripts/junction-endpoint-2026-08-03/` -- Phase 0 measurement scripts and raw
  output
- `/scratch1/rob/long-read-ecosystem/juncprobe/` -- Rust harness (`phase0`,
  `reweight`, `em` modes); ~10 s for a full 35 M-record BAM
