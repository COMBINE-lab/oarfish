# Modeling unannotated isoforms (`--model-unannotated-isoforms`)

**Status: experimental, off by default.** Genome (projection) mode only.

When the annotation is incomplete, reads from an isoform that is missing from the
GTF have nowhere correct to go. oarfish's default behaviour is already
conservative — a read whose splice structure matches no annotated transcript
fails projection and is dropped, contributing zero rather than being
misattributed — but that means its mass disappears with nothing recording why.

This flag turns those reads from silent drops into reported evidence.

```bash
oarfish --genome GRCh38.fa --reads reads.fq --annotation annotation.gtf \
        --seq-tech ont-cdna --model-unannotated-isoforms \
        --output sample
```

## What it does

1. **Per-locus novel latent state.** A read whose splice structure disagrees with
   *every* annotated candidate (≥ `--novel-min-misses` junction mismatches against
   all of them) gets an extra latent "novel isoform" candidate, shared by all such
   reads at the same locus. The locus is the read's *ambiguity component* — the
   connected component of transcripts linked by shared ambiguous reads — so **no
   gene annotation is required**, which matters precisely at the novel loci that
   have no gene label.
2. **Locus attribution for unprojectable reads.** Reads that project onto nothing
   are matched against the annotation by plain genomic *overlap*, ignoring splice
   compatibility, so their mass can be attributed to a locus. These reads still
   never enter the EM — no annotated estimate is perturbed — they are only
   counted.
3. **Two reports** (below).

## Why a single flag

The obvious worry with a feature like this is a pile of knobs the user must tune
to get a good result. Measurement says otherwise: **the defaults are already the
best setting for everything except the model itself.** See
[`annotation-omission-evaluation-2026-07-25.md`](annotation-omission-evaluation-2026-07-25.md).

| knob | best value | evidence |
|---|---|---|
| projection junction tolerance | **default (40/40/35)** | relaxing costs −0.005 to −0.022 CCC; it is doing real work |
| projection similarity threshold | **default (0.60)** | inert — a 12× reduction admits 41 reads out of 1.14M |
| novel-state odds/thresholds | **defaults** | `--novel-min-locus-reads 5` removes the MARD regression |

So there is nothing to bundle beyond "turn the model on", and this flag is simply
the supported way to do that. `--coverage-ablation annotation-omission` remains an
equivalent entry point reserved for the ablation harness.

## Reports

### `<output>.unexplained.tsv` — per-locus unexplained mass

One row per locus that acquired a novel state and attracted nonzero mass.

| column | meaning |
|---|---|
| `locus` | ambiguity-component id (arbitrary but stable within a run) |
| `gene` | gene label if the reference name carries one, else `.` |
| `n_transcripts` | annotated transcripts in the locus |
| `flagged_reads` | reads disagreeing with every annotated candidate |
| `unprojectable_reads` | dropped reads overlapping the locus (lower bound; see below) |
| `unexplained_mass` | EM mass assigned to the novel state |
| `annotated_mass` | EM mass on the locus's annotated transcripts |
| `unexplained_fraction` | `unexplained / (unexplained + annotated)` |
| `transcripts` | up to 8 member transcript names |

A high `unexplained_fraction` is the signal to act on: it says a substantial share
of the locus's reads fit no annotated isoform.

### `<output>.unprojectable.tsv` — per-transcript overlap of dropped reads

One row per annotated transcript whose exons overlap at least one unprojectable
read. This is the **complete** view: the locus report only covers loci that
produced junction-mismatch evidence, which is a small minority. On a 50%-holdout
test the locus report listed 75 loci while this report listed 15,604 transcripts.

Two caveats when reading it:

- Counts are **not normalized by expression**, so highly expressed genes dominate.
  Join against the `.quant` file to get a rate if you want to rank loci.
- A read overlapping several transcripts increments each of them. In
  `unexplained.tsv` the per-locus figure is therefore the **maximum** over member
  transcripts — a lower bound on distinct reads — rather than the sum, which would
  double-count.

## Interpreting the log line

```
unprojectable reads: 495,283 (112,009 attributable to an annotated locus, 383,274 intergenic)
```

"Intergenic" means the read overlaps no *surviving* annotated exon. With a partial
annotation (e.g. protein-coding + lncRNA only) this is typically the larger share
and is expected; it is not evidence of a missing isoform at an annotated locus.

## Known limitations

- **Recall is low.** On holdout simulations the locus flag has excellent precision
  (0.87–1.00, and 1.00 at `--novel-min-locus-reads` ≥ 3 for random holdouts) but
  recall saturates around 2–4% of genes that lost an isoform. Most orphaned reads
  never reach the junction-evidence channel because they fail projection first.
  Treat a flagged locus as trustworthy and an unflagged one as uninformative.
- **Whole-locus deletions regress.** When *every* isoform of a gene is missing
  there is no surviving sibling to correct, and flagging measurably hurts
  (`locus25`: mean |log2| error 0.3965 → 0.5632). The model assumes the locus is
  partially annotated.
- **Signal-to-noise is modest.** Under a fully correct annotation 0.87% of
  junction-informative reads still disagree with every candidate (alignment and
  projection error). The signal rises only to 1.75%/3.11% at 25%/50% holdout, so
  the novel state is deliberately capped by `--novel-odds-per-miss`.
- Global accuracy gains are small (+0.0005 to +0.0021 CCC); the benefit is
  concentrated at acted-on loci (10–19% error reduction) and in the reports.

## Tuning

| option | default | effect |
|---|---|---|
| `--novel-odds-per-miss` | 2.0 | odds multiplier per unmatched junction favouring the novel state |
| `--novel-min-misses` | 1 | junction mismatches required before a read is flagged |
| `--novel-min-locus-reads` | 5 | flagged reads a locus needs before it gets a novel state |
| `--novel-require-hits` | off | additionally require the read to *agree* at ≥1 junction |

Raising `--novel-min-locus-reads` trades recall for precision; the sweep in the
evaluation doc shows precision reaching 1.000 by 3–5 on random holdouts.
