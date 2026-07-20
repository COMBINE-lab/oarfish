# Evaluation checkpoint: EM acceleration and endpoint model

> **Superseded endpoint result:** the 2026-07-19 evaluation found that endpoint
> probabilities were trained but not enabled in the EM. Endpoint conclusions in
> this checkpoint are invalid; see
> [coverage-evaluation-2026-07-19.md](coverage-evaluation-2026-07-19.md).

Date: 2026-07-18. This is an implementation checkpoint, not a promotion
report. Commands used the release build from the working tree and fixed input
alignments or identical raw-read mapping options between model variants.

## Bundled alignment smoke test

Input: `eval/parity/small_txp.bam` (4,198 aligned reads; 9,579 retained
alignments), `--filter-group no-filters`, four threads, coverage disabled.

| EM method | M-step evaluations | Converged at 0.001 |
|---|---:|:---:|
| ordinary EM | 1,000 | no |
| SQUAREM | 223 | yes |
| DAAREM | 406 | yes |

The three output totals agree to floating-point precision (4,198). Their
truth-comparison metrics are indistinguishable at the displayed precision. The
available truth file describes a much deeper library than this small BAM, so it
is useful for numerical parity after library-size scaling but is not treated as
an accuracy benchmark.

## Local direct-RNA checkpoint

Input: `drna_test/hek_100k.fastq`, GENCODE v47 transcriptome, ONT direct-RNA
preset. Ordinary and repaired-logistic runs completed on 95,324 mapped reads.
The full endpoint run was repeatedly terminated by the execution environment
during raw-read mapping, before endpoint-model training. A 2,000-read subset
completed and trained on 781 unambiguous reads. This is recorded as an
unresolved resource/execution failure; it is not evidence for or against the
model's accuracy.

## Public SIRV checkpoint

Public ONT direct-RNA SIRV mixtures were downloaded to
`/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/public-sirv`:

| Mixture | Run | FASTQ MD5 | Compressed bytes |
|---|---|---|---:|
| E2 | SRR6058583 | `1146f3fbd028d1b6502cb4dd38b8453e` | 222,870,088 |
| E0 | SRR6058584 | `ff98d2f899dcf52774b55310fc70b268` | 120,341,594 |

The matching Lexogen reference FASTA MD5 is
`69cf02d2cd69ff118967e291c00df21a`. Transcriptome-mode runs against its 114
reference records completed for all three coverage modes. On E0, using the
known equal isoform mixture aggregated to the seven SIRV loci, ordinary EM and
endpoint were identical (CCC 0.948587, Pearson 0.978646, RMSE 0.940893);
logistic changed only in the sixth decimal place (CCC 0.948589, RMSE 0.940877).

This is only a locus-level sanity check: the distributed FASTA contains SIRV
loci plus ERCC records, while isoforms are represented in the matching GTF.
Genome-projection runs were terminated during mapping in this environment.
Therefore these results do not satisfy the transcript-level promotion gate.

## Decision

- Keep `--em-accel none` and `--coverage-model none` as defaults.
- Accept SQUAREM and DAAREM as opt-in inference features; both reached the same
  numerical solution with fewer map evaluations in the bundled smoke test.
- Retain `endpoint` as an experimental opt-in model. Do not implement or
  promote endpoint mixtures, predictive coverage, or segment models until a
  held-out, transcript-truth evaluation and the local full-data failure are
  resolved.
