# LongBench coverage evaluation: 2026-07-19

This report is the expanded biological holdout requested by the coverage-model
evaluation plan. All model defaults were frozen before this panel was examined.
The public [LongBench dataset](https://registry.opendata.aws/longbench/) contains
matched ONT PCR-cDNA, ONT direct-RNA, PacBio Kinnex, and Illumina measurements
from eight lung-cancer cell lines.

## Design

For each cell line and long-read technology, the first 50,000 FASTQ records were
aligned once to GENCODE v47 transcripts with minimap2 (`map-ont` for ONT and
`map-hifi` for PacBio). The same name-sorted BAM was quantified with `none`,
`adaptive`, and `auto`, with alignment filters disabled. Estimates were compared
with the matched public Illumina `NumReads` matrix after removing transcript
versions and rescaling the long-read estimates to the Illumina library total.
All metrics use the 251,488 shared transcript identifiers. This is an external
biological comparator, not molecule-level ground truth; protocol and GENCODE
version differences remain potential confounders.

The screen comprises 24 libraries, 72 quantifications, and 1.2 million long
reads. A predeclared 250,000-read confirmation tier then tested the largest
gains, the clearest regression, and the largest `adaptive`/`auto` separation.
The confirmation tier contains five libraries and 15 quantifications.

## Complete 50,000-read screen

Values are unweighted means across eight cell lines.

| Technology | Model | CCC | Pearson | RMSE | MARD | Coverage time | EM time |
|---|---|---:|---:|---:|---:|---:|---:|
| ONT cDNA | none | 0.27547 | 0.42684 | 7871.5 | **0.45681** | <0.001 ms | 1.062 s |
| ONT cDNA | adaptive/auto | **0.28217** | **0.43501** | **7797.1** | 0.45724 | 38 ms | 1.072 s |
| ONT dRNA | none | 0.49051 | 0.70388 | 6302.8 | **0.44849** | <0.001 ms | 1.063 s |
| ONT dRNA | adaptive | **0.49751** | **0.71265** | **6257.0** | 0.44925 | 39 ms | 1.077 s |
| ONT dRNA | auto | 0.49605 | 0.71101 | 6271.3 | 0.44938 | 112 ms | 1.080 s |
| PacBio | none | 0.51868 | 0.57926 | 4053.0 | **0.44799** | <0.001 ms | 1.151 s |
| PacBio | adaptive/auto | **0.52701** | **0.58924** | **4018.9** | 0.44853 | 44 ms | 1.170 s |

Coverage correction improves mean CCC, Pearson, and RMSE for all three
technologies. The CCC and Pearson changes are positive in every ONT library.
MARD is flat to slightly worse, showing that the gain is concentrated in the
moderate/high-count transcripts that dominate correlation and squared error;
coverage correction should not be claimed to improve every accuracy endpoint.

For four direct-RNA libraries the learned degradation hazard is the grid null
and `auto` reduces exactly to `adaptive`. Four receive a modest correction.
All eight remain better than `none`, but `auto` averages 0.00145 lower CCC than
`adaptive`, so the degradation layer is conservative but not yet beneficial on
this cross-platform comparator.

## 250,000-read confirmation

| Library | Auto - none CCC | Auto - none Pearson | Auto - none RMSE | Auto - adaptive CCC |
|---|---:|---:|---:|---:|
| PacBio H69 | +0.05190 | +0.04480 | -287.2 | 0 |
| PacBio H146 | +0.00910 | +0.00999 | -55.5 | 0 |
| PacBio H526 | +0.03404 | +0.04661 | -192.6 | 0 |
| ONT cDNA HCC827 | +0.01706 | +0.01639 | -223.9 | 0 |
| ONT dRNA HCC827 | +0.00424 | +0.00494 | -41.0 | -0.00428 |

The apparent PacBio H146 regression at 50,000 reads reverses at 250,000 reads;
the two large PacBio gains and the cDNA gain persist. The direct-RNA result also
confirms both sides of the screen: automatic coverage correction remains better
than no positional model, while its learned degradation adjustment remains
slightly worse than the simpler adaptive kernel. At 250,000 reads HCC827 learns
hazard 0.1/kb, intact fraction 0.616, and full correction weight.

Mean 250,000-read adaptive coverage time is 0.147 s for PacBio versus 5.51 s of
EM (2.7% of EM time), and 0.134 s versus 5.22 s for cDNA (2.6%). Direct-RNA
`auto` costs 0.488 s versus 4.62 s of EM (10.6%); degradation fitting is the
only material coverage-model runtime cost observed here. Mapping and BAM
parsing are excluded from these metadata timings.

## Decision

- The cross-fitted adaptive endpoint kernel now has positive evidence across
  all eight cell lines and all three technologies, including depth confirmation.
- Keep `auto` as the preferred experimental interface: it is identical to
  adaptive for cDNA/PacBio and never regressed below `none` in the ONT screen.
- Do not yet make coverage correction the global default. Matched Illumina is
  not transcript truth, low-count MARD does not improve, and direct-RNA
  degradation weighting still leaves measurable accuracy on the table.
- The planned degradation-calibration follow-up is complete. Internal scalar
  calibration failed the degradation-trajectory gate, while a three-component
  competing-risk kernel removed the LongBench penalty and retained robustness;
  see the [follow-up evaluation](competing-risk-degradation-evaluation-2026-07-20.md).

Machine-readable results are
`oarfish-evaluation-data/longbench/panel-50k-results.tsv` and
`panel-250k-results.tsv`; the sample manifest is
`benchmarks/longbench_coverage_panel.tsv`.
