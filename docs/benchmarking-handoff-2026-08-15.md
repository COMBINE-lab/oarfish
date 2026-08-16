# oarfish coverage benchmarking — handoff (2026-08-15)

Written to be picked up cold. It covers what the evaluation panels are, where
the data lives, what has been decided and why, and the accuracy leads that were
measured but not pursued.

Everything below was verified present on 2026-08-15.

---

## 1. Where things stand

`dev` is at `ac0248c` and in sync with `origin/dev`. The coverage-model program
that grew after the logistic model was measured against it on data with exact
read-level truth, and retired: `--coverage-model` now accepts only
`none|logistic`.

| change | evidence |
|---|---|
| removed rank blending, dominance pruning | -0.0878 / -0.0437 Spearman, worst family |
| removed alignment calibration, censoring (both had shipped **on**) | -0.000143 / -0.000233 on exact truth |
| removed `auto\|adaptive\|endpoint\|hybrid\|degradation` + 4 modules | `endpoint` -0.0428 (worse than no coverage model at all); `auto` -0.0004 at +25.9% wall |
| `--score-prob-denom` default 5 → 3 | improved 32 of 34 samples |
| added `--model-unannotated-isoforms` (genome mode) | detection precision 1.000; -10.7% error at acted-on loci |

Retired code is preserved, each branch carrying a "do not merge without new
evidence" note:

- `archive/coverage-extras-2026-08-03` — the four extras
- `archive/coverage-kernels-2026-08-03` — the five non-logistic kernels
- `archive/junction-endpoint-2026-08-03` — splice-junction endpoint term
- `wip/coverage-and-annotation-2026-07-25` — in-flight work; **three of its four
  threads were deliberately not ported** (endpoint grid geometry, further
  coverage experiments, evaluation-tier methodology)

> **History caveat.** `fef285b` (the WIP commit) is an ancestor of `dev` via the
> port merge, so those three unported threads read as "already merged" although
> their content was resolved away. `git merge wip/...` is a no-op. The archive
> branches are their real home.

---

## 2. The panels

The single most important fact for anyone benchmarking here: **the two panels
disagree in sign on essentially every arm.** Panel A ranks `none` above
`logistic`; Panel B does not. Any conclusion drawn from Panel A alone is
unreliable, and that is how the retired features were selected in the first
place.

### Panel B — exact read-level truth. Use this one.

Six transcriptome BAMs, 177,816 RefSeq identifiers, realistic depth
(320–338 reads per expressed transcript). NanoSim read names encode the
generating transcript (`NM_015477_2260_aligned_...`), so **per-read truth is
available**, not just per-transcript counts.

| sample | BAM (under `oarfish-evaluation-data/sim-panel/`) | `--seq-tech` | extra |
|---|---|---|---|
| nanosim-NA12878-cdna | `nanosim/NA12878_1DcDNA/cdna_T.bam` (11G) | `ont-cdna` | |
| nanosim-NA12878-drna | `nanosim/NA12878_directRNA/drna_T.bam` (11G) | `ont-drna` | `-d fw` |
| nanosim-H9-cdna | `nanosim/H9_1DcDNA/cdna_T.bam` (15G) | `ont-cdna` | |
| nanosim-H9-drna | `nanosim/H9_directRNA/drna_T.bam` (18G) | `ont-drna` | `-d fw` |
| tksm-RSII | `tksm/RSII_T.bam` (14G) | `pac-bio` | |
| tksm-SQ2 | `tksm/SQ2_T.bam` (13G) | `pac-bio-hifi` | |

Truth: `sim-panel/nanosim/ground_truth/*.csv`, `sim-panel/tksm/ground_truth/*.csv`.
Manifest: `rebench-2026-08-03/panelB_manifest.tsv`.

Upstream copies (NFS, mounted on this host — no ssh needed):
`/fs/cbcb-lab/rob/students/zahra/{nanosim_data,tksm_new_model}`.

**Truth-format traps.** NanoSim truth is headerless and unversioned. TKSM truth
has a `Transcript_ID/Count` header *and* versioned accessions. Strip
`\.\d+$` on both sides before joining.

### Panel A — proxy truth. Treat as a weak comparator.

- **28 LongBench samples**, `benchmarks/coverage_ablation_manifest.tsv`
  (8 cell lines x ONT-cDNA/ONT-dRNA/PacBio, plus 3 PacBio 250k). Truth is
  *matched Illumina*, quantified against a 252k annotation while the BAMs carry
  386k — score only over the shared reference (251,488), never the union.
  Absolute Spearman is 0.21–0.43; on this panel `none` beats `logistic` by
  +0.0055.
- **12 SIRV / spike-in samples**, `rebench-2026-08-03/panelA_sirv_manifest.tsv`,
  read mode. **SIRV E0 is equal-molar**, so truth variance is zero and every
  rank metric is undefined — report `NA`, not `0.0`. Only the 4 E2 samples
  support Spearman. The 2 `independent-sim` samples have exact counts but a
  near-degenerate truth (max 6 reads/transcript, 81% singletons).

### Genome-mode inputs

- `eval/parity/genome.bam` (2G) — spliced genome BAM, **NanoSim read names, so
  per-read truth**; 1,356,492 reads project under the full annotation
- `eval/GCF.pc_lncrna.matched.gtf`, `test_data/GCF_000001405.40_GRCh38.p14_genomic.fna`
- `test_data/cdna_ground_truth.csv` — the expressed set, for building holdouts
- GTFs for junction work: `rebench-2026-08-03/junc/refseq.gtf` (RefSeq, matches
  Panel B) and `junc/gencode.v47.gtf` (matches LongBench; reconciles perfectly,
  0 length mismatches, vs 4,925 for RefSeq)

> `/scratch1` is node-local and **not backed up**. Panel B is re-copyable from
> the NFS paths above; the LongBench BAMs and `eval/parity/genome.bam` are not
> known to exist elsewhere.

---

## 3. Harness

`rebench-2026-08-03/` (also committed under `scripts/rebench-2026-08-03/` and
`scripts/omission-validation-2026-08-04/`):

| file | purpose |
|---|---|
| `run_arms.py` | run N arms x M samples; alignment **and** read mode; captures wall/RSS |
| `metrics.py` | Spearman, Kendall tau-b, Pearson/CCC on log1p, RMSE, MARD |
| `analyze.py` | per-sample and per-family tables, deltas vs `logistic` |
| `make_holdout.py`, `score_omission.py` | annotation-holdout construction and detection scoring |
| `coverage-rebenchmark-2026-08-03.tsv` | 448-row per-sample deliverable |

`juncprobe/` is a Rust harness (built binary present) with `phase0`, `reweight`,
`em` and `ceiling` modes — ~10 s over a 35M-record BAM, versus ~10 min for the
equivalent Python. Use it for anything that must touch every alignment.

### Scoring rules that have already cost time

1. **Score over a hash-verified identical reference transcriptome.** `run_arms.py`
   md5s each arm's `.quant` key set and aborts if two compared arms differ.
2. **Renormalise to library size before comparing counts.** A holdout that drops
   35% of assigned reads deflates every surviving transcript ~3x. This inverted
   the sign of a whole validation once; the tell was a baseline signed error of
   -1.70 where the reference analysis had +0.43. **If your baseline's sign
   disagrees with the reference, stop and reconcile before interpreting.**
3. **`--seq-tech` is required even in alignment mode** — `bulk.rs` gated the
   PacBio kernel on it, and the original harness passed it.
4. **Zero truth variance means undefined, not zero** (SIRV E0).
5. Arm-vs-arm on one binary is deterministic. The ~0.005 "noise floor" quoted in
   the older docs is a *cross-branch* figure for the same nominal model; it does
   not apply to flag-on/flag-off on one build.

---

## 4. Decision rules (pre-registered; keep using them)

From `docs/coverage-auto-rebenchmark-plan-2026-08-03.md`:

1. A feature earns a default only by beating **`logistic`** — not `none` — on the
   union of both panels, with no dataset family regressing more than 0.01
   Spearman.
4. If the bare kernel does not beat `logistic`, the kernel should not be a
   default either.
5. **If the panels disagree, prefer Panel B and report the disagreement.**

Two amendments learned the hard way:

- **Pre-register a minimum effect size.** The junction-endpoint term passed the
  consistency criterion as literally written (6/6 samples improving) while
  delivering +0.00026 — about 1/35th of the projected magnitude. Consistency
  without magnitude is not a bar.
- **Rule 5 is not mechanical.** When Panel A is *strongly and uniformly*
  opposed rather than merely noisy, that is a conflict to investigate, not a
  tie to break by fiat. See the `responsibility-profiles` entry below.

---

## 5. Where the remaining headroom is

Measured on Panel B, `logistic`:

| | Spearman | headroom |
|---|---|---|
| as shipped | 0.9220 | — |
| + perfect **detection** (truth-zeros forced to 0) | 0.9481 | **+0.0261** |
| + perfect **abundance** (expressed set to truth) | 0.9765 | +0.0545 |

**1,812 false-positive transcripts carrying 0.156% of the read mass cost 0.026
Spearman** — as much as the entire benefit of coverage modelling.

### Lead A — a presence/absence component in the EM. Largest known.

Worth up to **+0.0261**, roughly ten times anything the retired stack delivered,
and orthogonal to which coverage model is used. Naive thresholds capture only
10% of it (+0.0026 at a 0.5-read floor) and aggressive ones lose 0.04–0.11,
because the classes overlap exactly where the false positives live: of
transcripts with est ≤ 1 and no unique read, 65% are truly zero but **35% are
real**.

The discriminating information exists — false positives have median 1.32
estimated reads versus 30.68 for true positives (23x), and 98.7% have no unique
read versus 61.7%. It needs a probabilistic component (zero-inflated or
spike-and-slab) emitting a per-transcript presence probability, not a cutoff.

**A 0.5-read floor is available today** for +0.0026 at no cost, as a stopgap.

### Lead B — `logistic` underperforms on PacBio. Repeatedly observed.

On TKSM SQ2, `none` (0.9501) beats `logistic` (0.9479). Across both panels
`logistic` loses to `none` on PacBio (-0.0025, 4/13), and the retired `auto`
and `hybrid` kernels both beat it on the two TKSM samples. No kernel beat `none`
consistently, so nothing was kept — but this is the most reproducible anomaly in
the series and the likeliest place a *new* model finds room. n=2 on exact truth
is the limiting factor: **more PacBio samples with exact read-level truth is the
single highest-value data acquisition.**

### Lead C — `responsibility-profiles`, unresolved and worth revisiting

On `wip/coverage-and-annotation-2026-07-25` as a `--coverage-ablation`. It
rebuilds the logistic coverage profile from cross-fitted abundance
responsibilities instead of counting every alignment once — the current
construction lets an ambiguous read deposit a full read of coverage on every
candidate it touches (**measured 5.09x inflation on H69 cDNA**). A real defect.

Measured on `logistic`:

| dataset | truth | delta | consistency |
|---|---|---|---|
| Panel B (NanoSim/TKSM) | exact | **+0.0023** | 6/6 |
| independent-sim | exact counts, **different simulator** | -0.0013 | 0/2 |
| SIRV-E2 | **real reads, known concentrations** | -0.0004 | 1/4 |
| Panel A LongBench | matched Illumina | -0.0030 | 0/28 |

Not adopted: one dataset family of four supports it, and the two that are
*real with real truth* do not. But the defect it fixes is genuine, so a
different correction for the same inflation may still be worth having. Note it
is **not** redundant with `--score-prob-denom 3` — the gain was larger at D=3
than at D=5.

### Ideas already tested and killed — do not re-propose without new evidence

| idea | result |
|---|---|
| length / effective-length correction | **oracle** gain -0.0001 Spearman; the bias is real but rank-preserving |
| ambiguity-gated coverage correction | 0.0000 at every gate threshold |
| splice-junction endpoint term | 9.9x likelihood ratio, +0.025 posterior mass/read, **+0.00026** accuracy |
| endpoint grid geometry corrections | +0.000186; the nucleotide-measure variant is **harmful** (-0.005) |
| `component-mass-conservation` | byte-identical on 6/6 — never engages |
| `polya-three-prime` | +0.000108 |
| `continuous-nested-guard` | helps only the broken pre-demotion baseline |

**The junction-endpoint result generalises and is worth internalising.** An
information-ceiling analysis found that **~73% of misassigned reads have the true
source and the winning decoy scoring *identically*.** Equal alignment score over
the read's span implies indistinguishable sequence, hence indistinguishable exon
structure — so *no* annotation-derived per-alignment feature can separate them,
regardless of formulation. The other ~27% already score higher for the true
source and lose to the abundance prior, where annotation is redundant. Reproduce
with `juncprobe ceiling` before investing in any per-alignment signal.

---

## 6. Open items

- **`--model-unannotated-isoforms` recall is ~1%** (89 of 9,720 genes) at
  precision 1.000. Precision is not the weak point; recall is. The gates are
  `--novel-min-misses`, `--novel-min-locus-reads` (default 5) and
  `--novel-require-hits`.
- **The main EM does not converge at the shipped cap** — 1,000 evaluations,
  `converged: false`, at threshold 0.001, for *every* arm including `logistic`,
  on 177k-transcript data. Surfaced incidentally; affects every configuration
  and has never been investigated.
- **`--coverage-model auto` now hard-errors.** Nothing in `.github/`, `justfile`
  or `scripts/*.sh` used it, but external workflows and the paper's methods were
  not checkable from here.
- **Panel B accuracy validation for the unannotated-isoform feature is
  genome-mode only**; the transcriptome-BAM panel cannot exercise it.
- The original `origin50`-style holdout scripts from the July evaluation were
  **not preserved**, so its holdouts can only be re-derived from their
  description, not reproduced exactly. `make_holdout.py` re-derives them.
