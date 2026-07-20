# Coverage-model evaluation: 2026-07-19

This report compares `none`, repaired `logistic`, and experimental `endpoint`.
All accuracy metrics were calculated by `scripts/evaluate_quant.py`. Estimated
counts were scaled to the truth library size for the tables below.

During this evaluation we found and fixed an endpoint integration bug: the
model learned probabilities, but the EM switch controlling their use remained
off. The results below were regenerated after the fix and are the first valid
endpoint results.

## Large synthetic benchmark

Input was the fixed transcriptome alignment
`eval/parity/txp_cpp.rh.bam`: 1,356,482 retained reads and 3,214,282 alignments.
Truth was `eval/subset_ground_truth.tsv`, matched without accession versions.
Using one BAM removes mapping variation between coverage models.

| Model | CCC | MARD | RMSE | Pearson | Spearman | EM evals | Converged |
|---|---:|---:|---:|---:|---:|---:|:---:|
| none | 0.956952 | 0.169178 | 89.4590 | 0.957010 | 0.819877 | 1,000 | no |
| logistic | **0.994990** | **0.142519** | **30.7349** | **0.994995** | 0.824666 | 1,000 | no |
| endpoint | 0.955198 | 0.153293 | 91.4663 | 0.955218 | **0.825089** | 1,000 | no |

The repaired logistic model is decisively best here: relative to no coverage,
RMSE drops 65.6% and MARD drops 15.8%. Endpoint reduces MARD by 9.4% and has the
best rank correlation, but slightly regresses CCC, Pearson, and RMSE. This
suggests endpoint evidence helps many low/moderate-abundance assignments while
making a smaller number of high-count errors worse.

## Small source-labelled synthetic direct-RNA set

Input was `drna_test/long_reads.fastq` (1,634 reads). The source transcript is
encoded in each read name; `scripts/truth_from_fastq_names.py` generated truth.

| Model | CCC | MARD | RMSE | Pearson | Spearman |
|---|---:|---:|---:|---:|---:|
| none | 0.981515 | 0.0001168 | 0.012606 | 0.981683 | 0.990737 |
| logistic | 0.981661 | 0.0001157 | 0.012555 | 0.981826 | 0.990737 |
| endpoint | **0.986311** | **0.0000848** | **0.010822** | **0.986403** | **0.993215** |

Endpoint is best on every reported metric, although this is a small and nearly
full-length dataset and all models are already highly accurate.

## Public ONT direct-RNA SIRV E0

The official Lexogen GTF and genomic FASTA were converted to a 171-sequence
transcript FASTA by `scripts/extract_gtf_transcripts.py`. The E0 truth contains
79 equal-molar SIRV transcripts. Results use the same first 50,000 reads of
SRR6058584 for every model; 29,396 reads were retained.

| Model | CCC | MARD | RMSE | Pearson | Spearman | EM evals |
|---|---:|---:|---:|---:|---:|---:|
| none | 0.639192 | 0.166302 | 0.529723 | 0.685358 | 0.787281 | 430 |
| logistic | 0.632240 | 0.166512 | 0.537734 | 0.679886 | 0.787281 | 415 |
| endpoint | **0.732652** | **0.141257** | **0.425909** | **0.760328** | **0.831883** | 420 |

Endpoint improves all five metrics: CCC rises 14.6%, RMSE falls 19.6%, and MARD
falls 15.1%. Logistic is slightly worse than no coverage on this public set.

The same 50,000-read comparison was completed for SIRV E2 (SRR6058583), but a
verified transcript-level E2 concentration table was not available locally and
could not be obtained from the official resources inspected. Those outputs are
retained for future scoring rather than evaluated against inferred or circular
truth.

## Interpretation and next decision

- The repaired logistic model is very strong on the large synthetic benchmark,
  but its small SIRV regression argues against making it the default.
- The joint endpoint model generalizes to the public E0 data and the small
  source-labelled synthetic set. Its mixed result on the large synthetic set
  indicates that a single global endpoint distribution is not sufficient.
- The next endpoint iteration should target robustness rather than simply more
  flexibility: cross-fitted training, adaptive shrinkage by stratum/cell, and a
  learned mixture with a neutral component are the most directly motivated
  changes. Evaluation should specifically stratify the large synthetic errors
  by abundance and endpoint cell to identify the high-count regressions.
- Keep all coverage modes opt-in until a second public technology/dataset and
  the E2 truth table are scored.

## Artifacts

Outputs are under:

- `/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/synthetic-coverage`
- `/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/synthetic-drna-small`
- `/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/public-sirv`
