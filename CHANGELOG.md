# Changelog

All notable changes to `oarfish` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> Release tooling note: `dist` (cargo-dist) extracts the section whose header
> matches the release tag and uses it as the GitHub Release body. Cut releases with
> `./bump_and_publish.sh <version>`, which renames the `## Unreleased` heading below
> to the version (e.g. `## 0.10.0 - 2026-06-06`) and starts a new one above it.

## Unreleased

### Added

- **A warning when the EM stops at its iteration limit before converging**, in the log and
  in the `warnings` array of `meta_info.json`, for every bulk EM: `em_not_converged` for the
  default equivalence-class EM and the per-read EMs, `truncation_em_not_converged` for the
  truncation EM, each with its own advice. Its `fraction` is the EM's last convergence
  statistic and its `threshold` is `--convergence-thresh`. Estimates do not change. The
  default EM can stop just short of the threshold at 1,000 iterations on inputs with many
  transcripts, so some runs that gave no warning in 0.12.1 now give this one.
- **`em` in `meta_info.json`**: which EM ran (`standard` or `truncation`), its iterations and
  limit, whether it converged, and its last convergence statistic, for every bulk run. The log
  also reports the iteration count (`EM finished in ... (N iterations)`).
- **spearfish 0.3.0: alignments are weighed by spearfish's quantification score.** With
  `--mapper spearfish`, oarfish weighs a read's hits by `Hit::quant_score`, and from a RAD file
  by the `qscore` alignment tag when the file has it (the `score` tag otherwise, so files from
  earlier spearfish versions quantify as before). The quantification score is the mapping score
  plus an end fit for reads without a poly(A) tail, which favours the candidates whose annotated
  3' end the read reaches, and spearfish's orientation evidence for antisense pairs. Which hits
  a read has is unchanged; only their weights change. The score filters see the same score.
- **`--end-fit`, `--end-fit-tol` and `--end-fit-from chain|alignment`** (with `--mapper
  spearfish`): the end fit's points (0: off), its tolerance in compressed bases, and how where
  a read ends is measured. Each defaults to the preset's value. With `--model-truncation`, the
  end fit measures every candidate by alignment (spans are extended anyway) unless
  `--end-fit-from chain` is given.
- **`--orient-evidence quant|mapping|off`** (with `--mapper spearfish`): what spearfish's
  orientation evidence for antisense pairs does. `quant` (the preset's, so the default) moves
  only the quantification score; `mapping` moves the mapping score and decides the orientation
  (spearfish's `--no-orient-quant`); `off` turns it off, as spearfish 0.2.1 mapped (spearfish's
  `--no-orient-same-exact --orient-tail-after=80`). With `--end-fit 0 --orient-evidence off`,
  alignments are weighed by their mapping scores alone, and quantification is as with
  spearfish 0.2.1. `--end-fit 0` alone leaves the orientation evidence on. Recorded as
  `spearfish_settings.orient_evidence`.
- **A warning for spearfish indexes without gene groups** (`index_without_gene_groups`, in the
  log at warn level, also with `--quiet`, and in `meta_info.json`): when every transcript of the
  index is its own gene, isoforms compete one by one for spearfish's per-gene candidate places
  and ties among them go by FASTA order. The warning says how to fix it: build the index with a
  transcript-to-gene table or with gene inference.
- **`spearfish_settings.index_genes` in `meta_info.json`**, for every spearfish run: where the
  index's gene labels came from (`source`: `fasta` or `t2g`, both provided; `inferred`; `none`;
  or `unknown`), and its `genes`, `transcripts` and `has_gene_groups`, from the index manifest.
  For RAD input `source` is `unknown`: spearfish 0.3.0 RAD headers do not record it.
- **`--t2g FILE` and `--infer-genes`** (with `--mapper spearfish`, when oarfish builds the index
  from `--annotated`/`--novel`): gene labels from a transcript-to-gene table, or inferred from
  shared sequence, as `spearfish index --t2g` / `--infer-genes` give them. The index manifest
  records the source; an index oarfish builds in memory now has it too.
- **`spearfish_settings` in `meta_info.json`**, and a line in the log at the start of the run:
  with `--mapper spearfish`, and for RAD input, the settings the alignments were made with:
  `source` (`in_process` or `rad_header`), `preset`, `extend_ends` and `extend_ends_reason`
  (`preset`, `model_truncation`, `preset_off`, `rad_header`), `strands`, `quant_scores`, and
  the end fit (`end_fit`, `end_fit_tol`, `end_fit_from`), `orient_quant` and
  `orient_evidence`; for a RAD file
  also `mapper` and `mapper_version`, with the end fit read from the mapper's parameters. Spans
  decide how much of a read counts as aligned, and so which reads `--min-aligned-fraction`
  drops.
- **Trimming guidance for the aligned-fraction filter** (docs site and README): on untrimmed
  or older direct cDNA reads, `--min-aligned-fraction` can drop a large share of the reads, and
  trimming adapters and primers (keeping the poly(A) tails) reduces that. For gene-level work
  `--filter-group no-filters` or a lower `-m` may suit; for transcript-level estimates the
  default agreed better with short-read estimates in our tests. The `aligned_fraction_filter`
  warning now says the same and links to it. The filter's default is unchanged.

### Changed

- **spearfish 0.3.0** (`spearfish`, `long-read-rad` 0.3). Mapping output is as with 0.2.1;
  default quantification changes where the quantification score differs from the mapping score.
  `--mapper spearfish` runs get a quantification score in every preset except the `-fast`
  ones, and RAD files written by `spearfish map` 0.3.0 carry `qscore`.
- **The truncation EM's convergence test leaves out transcripts with 0.01 reads or fewer**
  (salmon's cutoff; it was 1e-5 reads), and still counts relative changes in either
  direction. Abundances decaying toward zero kept the test from passing for 2,300-2,700
  iterations, so every `--model-truncation` run stopped at the 1,000-iteration limit. It now
  converges after about 1,050-1,250 iterations. Estimates move by less than 0.02% of the reads,
  and spike-in and simulation accuracy is the same to four digits.
- **With `--model-truncation`, `--max-em-iter` defaults to 5,000** (1,000 otherwise, as
  before), since the truncation EM needs a little more than 1,000 iterations to converge. An
  explicit `--max-em-iter` is used as given. The extra iterations cost little: the truncation
  EM took 0.03-0.13 s longer on our test sets.
- **The `--model-truncation` help text says that, with `--mapper spearfish`, it extends
  alignment spans to the read ends** for the `ont-drna` and `pac-bio-hifi` presets (as it did
  in 0.12.1; the behaviour is unchanged), and so changes which reads `-m` drops.

## 0.12.1 - 2026-11-13

### Added

- **Warnings about the run, in the log and in `meta_info.json`.** When a bulk run finishes,
  oarfish checks the read counts of the discard table and warns when
  - more than 10% of the mapped reads were dropped by the aligned-fraction filter
    (`aligned_fraction_filter`; reads that still carry adapters or primers are the usual
    cause),
  - the strand filter (`-d fw|rc`) left more than 10% of the mapped reads without an
    alignment (`strand_filter`),
  - the other alignment filters did (`alignment_filters`),
  - more than 30% of the reads did not map under `--map-strand fw|rc`
    (`map_strand_unmapped`; the reads are probably not oriented), or
  - more than half of the reads did not map (`low_mapping_rate`).

  Each warning is logged, and recorded in a new `warnings` array of `meta_info.json` with its
  code, the share of reads, the threshold and the message. The checks use counts that were
  kept already; estimates are not affected.
- **`discard_ori_reads` in the discard table**: the reads all of whose alignments the strand
  filter removed (they are also counted in `no_valid_aln`, as before).
- **A guide on recommended practices by sequencing technology** (docs site): trimming of ONT
  cDNA reads (adapters and primers, not the poly(A) tails), when to use `--map-strand fw` and
  `-d fw` for PacBio reads, and what each warning means.

### Changed

- In the genome modes, a read with no usable projected alignment is now counted in
  `no_valid_aln` of the discard table (it was not counted anywhere). Estimates do not change.

- **`--mapper spearfish` embeds spearfish 0.2.1** (was 0.2.0; spearfish's changelog has the
  details):
  - a read that ends where a transcript ends, followed by a base or two that align to nothing
    (adapter, primer, basecalling junk), no longer prefers a longer isoform that those bases
    can be aligned into as mismatches. Such a read now ties between the two isoforms, and the
    EM decides. Estimates change slightly for genes with isoforms that are a prefix of one
    another;
  - mapping takes a little less CPU time, with the same mappings otherwise (the index data a
    read's lookups need is asked for ahead of use).
- **The per-thread arenas are on by default with `--mapper spearfish`**, as they are in
  spearfish's own command line from 0.2.1. They were opt-in (`--bump-alloc`) in 0.12.0. A run
  takes a little less time, for a few MB more memory a thread; what is mapped and the
  estimates do not change. `--bump-alloc=off` (or `none`) turns them off, and
  `--bump-alloc=<size>` still sets the largest allocation they serve. Other inputs and
  backends are not affected.

### Fixed

- **`--write-assignment-probs` reports a probability for each transcript, not for each
  alignment.** A read with several alignments to one transcript had an entry for each of them
  in the `.prob` file (the same transcript listed more than once, each with part of its
  probability). Those are now added, so that every transcript is listed once with the
  probability that the read comes from it, and `--display-thresh` applies to that sum.
- **`ambig_info.tsv` counts reads, not alignments.** A read with several alignments to one
  transcript was counted once for each in that transcript's total, and was counted as
  ambiguous even when it aligned to no other transcript. It now counts once, and as unique in
  that case. Estimates do not change.
- Typos in the help text (`--min-aligned-fraction`, `--write-assignment-probs`), in two log
  messages and in the README.

## 0.12.0 - 2026-11-12

### Added

- **`--mapper spearfish` writes its index with `--index-out`** (also with `--only-index`), as a
  spearfish `.sfi` with its build manifest next to it as `<path>.sfi.json`. Previously a saved
  index had to come from the separate `spearfish index`. Quantifying from a saved index
  (`--index transcripts.sfi`) skips reading the FASTA and building, and on a human transcriptome
  lowers peak memory by about 2.6 GB (TKSM SQ2, 14M reads, 16 threads: 79 s to 72 s and 6.3 GB
  to 3.7 GB) with the same estimates. The index carries the reference's seqcol digests in oarfish's
  footer, as rammap indices do, so the run records the same reference signature as building in
  process. An `.sfi` given to the rammap backend is now rejected with a hint to pass
  `--mapper spearfish`.
- **`--map-strand <fw|rc|both>` chooses the read orientation spearfish maps** (with
  `--mapper spearfish`). The default is the preset's: both orientations, and the read as given
  for `ont-drna`. For reads already oriented to the transcript strand (Iso-Seq FLNC reads, cDNA
  reads oriented by their primers), `--map-strand fw` maps the one orientation, which saves the
  work of the other. It is not a filter: a read of the other strand can still have a weak hit
  in the mapped orientation, so it is for oriented reads only; `--strand-filter` remains the
  way to keep one orientation's alignments of reads mapped in both. Given without
  `--mapper spearfish`, it is an error.
- **`--bump-alloc` turns on per-thread arenas for the spearfish backend** (with
  `--mapper spearfish`). The mapper's small allocations (up to 8 KB) then come from memory
  that belongs to each mapping thread instead of from the general allocator: a little less
  mapping time, for a few MB of memory a thread, with the same mappings. Off by default. A
  size in place of 8 KB is given as `--bump-alloc=1KB` (at most 32 KB), and
  `--bump-alloc=off` or `none` is the same as leaving the option out. The arenas are the [`nudgy`](https://crates.io/crates/nudgy) crate, the allocator behind
  `spearfish map --bump-alloc`: oarfish's global allocator is now `nudgy::Bump` around
  mimalloc (or around the system allocator without the `mimalloc` feature), and without the
  option it forwards every call, so nothing else changes.

### Changed

- **`--mapper spearfish` embeds spearfish 0.2** (was 0.1), whose default presets map
  differently (spearfish's changelog has the details):
  - a read with only a chance match to the transcriptome (reads from unannotated loci,
    contaminants) is left unmapped, where it used to be placed on whichever transcript
    shared a seed with it, so fewer reads count as mapped;
  - ONT cDNA and direct RNA reads whose first chain leaves them in doubt get more work
    (re-chaining at spearfish's rescue level, a whole-transcriptome table of rescue seeds,
    antisense transcripts of the same locus), which places more reads correctly at error
    rates of 5% and above. ONT cDNA reads take 1.1 to 1.15 times the mapping CPU time, and
    mapping holds about 1.5 GB more of the index in memory for a human transcriptome;
  - HiFi reads that end in a poly(A) tail in one orientation are mapped in that orientation
    first, and in the other only when the first leaves most of the read uncovered.

  An index written with `--index-out` now holds the tables for the second item, and is a
  little over twice the size it was (2.8 GB for a human transcriptome). An index written by
  an earlier oarfish, or by spearfish 0.1, still works: the tables are then made at each
  start.
- **spearfish hits are weighted by `exp(-d / 2)`**, where `d` is how far a hit scores below
  its read's best; it was `exp(-d)`. With `--mapper spearfish`, `--score-prob-denom` (in
  tenths of a spearfish score point) therefore defaults to 20, where it was 10. This is the
  scale spearfish 0.2 suggests for its scores (`spearfish::SCORE_DENOM`): oarfish now takes
  its in-process default from there, as it already took a RAD file's from the file, so both
  routes weigh hits alike. The previous scale suited simulated samples in which every
  transcript is expressed; when most are not, as in real samples, it gave unexpressed
  relatives of the expressed transcripts too many of the reads they share.
  `--score-prob-denom 10` restores it.
- **Read mode uses upstream `rammap-core` 1.2.0** instead of our temporary fork
  `rammap-core-temp`, now that the per-thread DP cache cap (jwanglab/rammap#9) is released
  upstream. Upstream redefines `is_supplementary` as the SAM FLAG 0x800 decision (every
  non-secondary mapping after a read's first); in the fork it marked z-drop chain splits only.
  Read mode's supplementary filter therefore now matches alignment mode, which filters on
  0x800. On 2M simulated ONT cDNA reads this drops 5 more mappings in transcriptome read mode,
  and genome read mode is unchanged. Estimates are otherwise identical, and rammap indices
  saved by 0.11 load unchanged.

### Performance

- **The EMs run on a connected-component layout** (after salmon's component
  engine). Reads link transcripts only within connected components, which the EM
  never couples; each component's classes are stored together, ordered by their
  heaviest target, and run tiles of up to 64 iterations while in cache;
  single-target classes become a per-transcript constant, and transcripts no read
  references are skipped. The stopping rule stays global and exact (a tile is
  replayed up to the round where it would have stopped). Wall time, EM time and
  peak memory, 16 threads, 0.11.0 to now:

  | | NanoSim NA12878 cDNA | TKSM SQ2 |
  |---|---|---|
  | default model | 13.0 s to 11.1 s; EM 3.3 s to 0.9 s | 20.5 s to 17.1 s; EM 4.9 s to 1.9 s |
  | coverage model | 20.9 s to 13.6 s; EM 9.6 s to 2.6 s; 2.2-2.5 GB to 1.8 GB | 35.0 s to 25.8 s; EM 14.8 s to 5.7 s; 3.79 GB to 3.52 GB |
  | truncation model | 25.4 s to 16.1 s; EM 16.4 s to 6.7 s; 1.96 GB to 1.69 GB | 41.4 s to 28.9 s; EM 28.4 s to 15.6 s; 3.80 GB to 3.21 GB |

  - The default and coverage models' estimates are unchanged up to summation order
    (within 2e-6 reads on these simulations), and the default model's do not
    depend on the thread count.
  - The coverage model no longer reorders the alignment store before the EM (the
    layout keeps that locality), which removes the copy that set its peak memory.
  - The truncation model's frozen-parameter phase runs SQUAREM on each component
    separately, so its estimates move slightly: on the six Panel B simulations,
    Spearman changes by -0.0010 to +0.0001 and expressed-only Spearman, MARD and
    log-Pearson by at most 0.0004. Its layout is built directly from the frozen
    parameters, without the frozen-weight table.
  - Applies to the equivalence-class EM (without `--presence-model`), the per-read
    EM above 4 threads, and the truncation EM (without `--presence-model`).
    `OARFISH_EM_COMPONENTS=0` restores the previous layouts.
- **Faster bootstraps.** Replicates share one component layout and differ only in
  their read counts; truncation-model replicates fit the positional parameters on
  the same 1M-read sample as the point estimate. NanoSim NA12878 cDNA, 16 threads:
  8 replicates 18.5 s to 4.9 s (default model) and 115 s to 29 s (coverage model);
  4 replicates 123 s to 25 s (truncation model). Replicate spread is unchanged.
- **Faster BAM parsing and EM inner loops.** The BAM is read by a BGZF pipeline of
  oarfish's own: each of `-j` workers inflates a block (with libdeflate) and reads its
  records in place, while a second thread filters the alignments and fills the store (and,
  under `--model-coverage`, a third adds the coverage). The equivalence classes are built
  in parallel, the header's seqcol digest is streamed and computed during the parse, and
  the component EMs' and the truncation fit's inner loops do less work per class. Wall
  time, 16 threads, from the component-layout EM above to now (median of 3):

  | | NanoSim NA12878 cDNA | TKSM SQ2 |
  |---|---|---|
  | default model | 11.0 s to 3.6 s | 16.8 s to 5.4 s |
  | coverage model | 13.4 s to 4.4 s | 25.0 s to 7.5 s |
  | truncation model | 15.8 s to 7.5 s | 28.3 s to 14.1 s |

  - At 4 threads, NanoSim NA12878 cDNA with the default model: 21.5 s to 13.0 s.
  - Peak memory is within 1.5% on these runs. On the smaller Panel B simulations it
    rises by up to 0.6 GB (H9 dRNA with the coverage model: 2.26 GB to 2.86 GB).
  - The estimates are unchanged except for summation order: the per-read EM orders each
    block's classes by size (NumReads move by at most 2.3e-8), and the truncation
    model's estimates move by rounding (on the six Panel B simulations, Spearman changes
    by -4e-5 to +7e-5 and MARD by at most 2e-5).
  - libdeflate's C source is vendored by `libdeflate-sys` and built with the `cc` crate,
    so the build still needs only a C compiler.

### Dependencies

- **No open cargo-audit advisories.**
  - arrow2 is replaced by `parquet` 60.1 for the `.infreps.pq` bootstrap file. The schema is
    unchanged, and the file reads the same in pyarrow. This clears RUSTSEC-2025-0038, and,
    through parquet2, RUSTSEC-2026-0403.
  - sendable-swapvec is replaced by a small built-in spill store for the read names that
    `--write-assignment-probs` keeps (RUSTSEC-2026-0041, lz4_flex).
- **No git dependencies.** seqcol_rs and kders come from crates.io (0.4.2 and 0.1.1).
- **Bumps:**
  - noodles: bam 0.98, sam 0.93, bgzf 0.52;
  - bramble-rs 0.2.0, needletail 0.7.4, niffler 3.0.2;
  - sha2 0.11, base64 0.23, tabled 0.22, statrs 0.19, itertools 0.15.

  Outputs and seqcol digests are unchanged.

### Fixed

- **`--short-quant` no longer loses reads.** The EM's updates are multiplicative, so a
  transcript whose short-read abundance was 0 stayed at 0, and the reads aligned only to such
  transcripts were dropped from the estimates (how many depended on the EM path and the thread
  count). Those transcripts now start at a negligible abundance (1e-8 of the mean), so every
  path keeps all reads and the paths agree.
- **Bootstraps no longer panic on input with no aligned reads**; the replicates are all zero, as
  the point estimate is.
- `--convergence-thresh` has its own help text (it repeated `--max-em-iter`'s), and read mode's
  mapping-time log line names the backend in use (it said rammap with `--mapper spearfish` too).

## 0.11.0 - 2026-10-24

### Added

- **`--mapper spearfish`: [spearfish](https://github.com/COMBINE-lab/spearfish) as a read-mode
  backend.** Maps reads to the transcriptome in process with spearfish (crates `spearfish` and
  `long-read-rad`, 0.1), building its index from `--annotated`/`--novel` or opening a `.sfi`
  given as `--index`. The preset follows `--seq-tech`; hits are weighted by their score
  difference from the read's best, and with `--model-truncation` spans are extended to the read
  ends. Transcriptome read mode only.
- **RAD input: `oarfish -a mappings.rad`** quantifies the `bulk_long_read` RAD files that
  `spearfish map --rad` writes, taking the transcripts, score scale, preset and span settings
  from the file.
- **`--model-truncation`: a learned placement model.** Each alignment's score probability is
  multiplied by the probability of where it lands on the transcript under a mixture of
  3'-anchored, 5'-anchored, full-length and uniformly placed reads, whose weights, read-length
  laws and end-gap laws are learned in the EM; the fitted mixture is logged. It adapts to the
  protocol without settings. Compared with the default model, on a direct-RNA spike-in run it
  raises isoform-level Spearman from 0.820 to 0.882; against a spike-in reference padded with
  fabricated isoforms it cuts the mass placed on them from 9.4% to 0.5%; on simulated reads
  Spearman goes from 0.953 to 0.956. Against short-read TPM on 16 matched samples, it agrees
  better than the default model on all 8 direct-RNA samples and on 6 of 8 cDNA samples (all 8
  on expressed transcripts). It works with bootstraps and in genome mode, and conflicts with
  `--model-coverage`. (Supersedes the 3'-anchored model proposed in #75, which mis-ranked
  isoforms whenever reads were not 3'-anchored.)
- **`--presence-model spike-slab|bayes` (experimental).** A per-transcript presence posterior
  that discounts low-evidence transcripts during the EM; reduces false positives at some cost to
  true low-abundance transcripts. Off by default.
- **`meta_info.json` records the oarfish version and the run's start and end
  times**, in the same fields and format as salmon (`oarfish_version`,
  `start_time`, `end_time`; e.g. `"Tue Jan 15 11:10:00 2019"`, local time), for
  bulk and single-cell runs.

### Changed

- **Coverage model (`logistic`): corrected bin arithmetic.** Coverage was filled on bins of rounded
  width but read weights were looked up on a fixed 100-nt grid, so the two disagreed about where
  bins lie; the bin where each read ends was dropped by accident; and the overlap was computed in
  unsigned integers that could underflow. Both now use one exact bin layout (`len / nbins`), with
  the one behaviour the accident encoded made a stated rule: a read's partly covered end bin is
  left out of both its coverage and its weight when it fully covers another bin (counting it
  costs accuracy). Versus 0.10.3's coverage model, Spearman moves by at most 0.001 on simulated
  reads and spike-ins, except for a gain of 0.013 against an over-annotated spike-in reference;
  agreement with short-read TPM improves on 12 of 17 matched samples.
- **`--coverage-model none|logistic|truncation`** selects the positional model in one flag;
  `--model-coverage` and `--model-truncation` remain as shorthands. Contradictory combinations,
  and the truncation model in single-cell mode (previously ignored silently), are errors.
- **`--score-prob-denom` defaults to 3 (was 5).** Sharper weighting of a read's alignments toward
  the best-scoring one. On simulated reads with exact truth it raises Spearman slightly (0.950 to
  0.953 for the default model, 0.954 to 0.956 with `--model-truncation`); against short-read TPM
  it improved 27 of 28 matched samples in a held-out sweep. On spike-ins it raises Spearman but
  slightly worsens MARD. `--score-prob-denom 5` restores the previous behaviour exactly.
- **Tail-neutral alignment scores, for all models.** The score credit a read's poly(A) tail earns
  by aligning to an annotated A-run at a transcript's end is removed, using the aligner's match
  score (from the minimap2 `@PG` line or the rammap preset). On simulated reads with poly(A)
  tails this improves the default model (Spearman 0.948 to 0.950, and 22% fewer false
  positives, with `--score-prob-denom` held at 5); outputs on data without poly(A) tails are
  unchanged.
- **EM engines.** Without a per-read model term, the EM runs over equivalence classes (CSR
  layout); per-read models use an atomic-free CSR layout (`OARFISH_READ_EM=local|atomic` selects
  the alternatives; `--per-read-em` forces the per-read EM). Estimates match the previous
  release up to summation order.
- **Reads are ordered by likeliest transcript before the per-read and truncation EMs**, which
  speeds them up; `--write-assignment-probs` still lists reads in input order.

### Performance

Versus 0.10.3, from transcriptome BAMs of 35M and 75M alignments at 16 threads:

| | time | peak memory |
|---|---|---|
| default model | 15–16× faster | 27–33% lower |
| coverage model | 9–10× faster | 2–22% lower |
| truncation model (new) | about 2× the default model's time | about the coverage model's |

- **BAM parsing about 4.5× faster than 0.10.3.** Parsing is serial (only BGZF decompression is
  threaded), so per-record work is wall time: records are now read undecoded into reused
  buffers, the CIGAR, tags and sequence are read from their bytes only where a filter needs
  them, and records are no longer cloned.
- **Lower peak memory**: the per-alignment record shrinks from 24 to 16 bytes (an unused field
  removed); coverage probabilities exist only under the coverage model; the per-read EM reads the
  alignment store instead of copying it; the truncation EM keeps its frozen weights in store order.
- **EM engines**: the equivalence-class and per-read EMs run without atomics (see above); the
  truncation EM fits its placement parameters on a 1M-read sample when there are more reads,
  then holds them fixed and accelerates the abundance updates (SQUAREM); estimates match a fit on
  all reads (`OARFISH_TRUNC_SAMPLE=0`).
- **Deterministic**: the truncation EM (and the spike-slab presence model) give byte-identical
  output from run to run and across thread counts.
- All parallel work stays within `--threads` (the store reorder and coverage normalization
  previously used every core).
- Genome-BAM and raw-read modes are as fast or faster than 0.10.3 with the same memory (genome
  BAMs about 2.3× faster).

## 0.10.3 - 2026-07-16

Patch release improving the numeric precision of the `--write-assignment-probs`
output.

### Fixed

- **Assignment probabilities in the `.prob`/`.prob.lz4` file could be reported
  as `0.000`.** The values were always printed with 3 decimal places, but the
  default `--display-thresh` (`1e-6`) admits assignments far smaller than
  3-decimal resolution, so a probability that passed the reporting threshold
  could still round to `0.000`. The printed precision now tracks
  `--display-thresh` (`decimals = ceil(-log10(display_thresh))`, clamped to
  `[3, 9]`): the default now prints 6 decimals, and an assignment that passes
  the threshold is never shown as zero. Thresholds of `1e-3` or larger are
  unchanged (still 3 decimals), and the `.prob` file format is otherwise
  identical. See [#71](https://github.com/COMBINE-lab/oarfish/issues/71).

## 0.10.2 - 2026-07-16

Patch release fixing raw-read mapping of ONT direct-RNA data basecalled in RNA
(uracil) mode.

### Fixed

- **ONT direct-RNA reads containing `U` (uracil) mapped at a near-zero rate in
  raw-read mode (`--reads`).** Reads basecalled in RNA mode carry `U`/`u`
  instead of `T`/`t`. The `rammap` seed encoding recognizes only `A/C/G/T` (`U`
  encodes as ambiguous), so `U`-containing k-mers never matched the `T`-based
  transcriptome index and ~99% of reads produced no mapping — even though the
  same reads map normally through minimap2 or the `rammap` CLI, both of which
  normalize `U` internally. oarfish now converts `U`→`T` (and `u`→`t`) once at
  read ingestion, so direct-RNA reads map as expected. BAM-input modes
  (`--alignments`, `--genome-alignments`) were never affected. See
  [#70](https://github.com/COMBINE-lab/oarfish/issues/70).

### Changed

- The raw-read discard table now reports two additional per-read outcomes —
  **read had no mapping** (the mapper returned no alignment record) and **read
  had no valid alignment** (alignments were found but none had a positive best
  score). These reads were previously dropped without being counted, so the
  table did not reconcile to the total number of reads; the new rows make the
  cause of unmapped reads visible (and were what pinpointed #70).

## 0.10.1 - 2026-07-11

Patch release fixing a `cargo install oarfish` build failure. No behavior change.

### Fixed

- **`cargo install oarfish` failed to compile against `bramble-rs` 0.1.5.** A
  fresh install re-resolves dependencies and picked up `bramble-rs` 0.1.5, whose
  `FastaDb::from_seqs` had switched to an internal deterministic-hasher map,
  rejecting the default-`RandomState` `HashMap` oarfish passed (a type mismatch
  in `src/util/aligner.rs`). See [#68](https://github.com/COMBINE-lab/oarfish/issues/68).
  Fixed by requiring `bramble-rs` 0.1.6, which makes `from_seqs` hasher-agnostic
  and exports `SeqMap` + `FastaDb::from_seq_map`; oarfish now builds bramble's map
  type directly and hands it over with no rehash. (`bramble-rs` 0.1.5 was yanked.)

## 0.10.0 - 2026-06-06

This release replaces the read-mode aligner and adds genome-based quantification.
Most users of alignment/BAM input are unaffected; users of raw-read mapping mode
should read the **Breaking changes** below.

### Breaking changes

- **Read-mode mapper swapped from `minimap2-rs` to the pure-Rust `rammap`
  backend.** The `minimap2-rs` backend and the `rammap` Cargo feature that gated
  the two have been removed; `rammap` is now the sole raw-read mapper. Quantification
  from raw reads (`--reads`) may differ slightly from previous `minimap2`-based
  runs. BAM-input modes (`--alignments`, `--genome-alignments`) are **unaffected**
  and still accept output from minimap2, pbmm2, STAR, etc.
- **`--thread-buff-size` removed.** It only configured `minimap2`'s per-thread
  `cap_kalloc` arena. It is replaced by **`--dp-cache-cap-mb`**, which bounds the
  per-thread DP alignment-scratch buffer the rammap mapper retains (default
  128 MB; also settable via the `RAMMAP_DP_CACHE_CAP_MB` environment variable;
  `0` disables the cap).

### Added

- **Genome-based quantification modes** (projection onto annotated transcripts via
  [bramble](https://github.com/COMBINE-lab/bramble)):
  - *Genome read-projection* (`--reads --genome --annotation`): spliced-align raw
    reads to the genome and project their alignments onto the transcripts.
  - *Genome-alignment projection* (`--genome-alignments --annotation`): project an
    existing spliced **genome** BAM (from minimap2/pbmm2/STAR/…) onto the transcripts.
  - Soft-clip **rescue**, on by default, recovers discriminating sequence from
    soft-clipped read ends (`--no-rescue` to disable; `--genome-fasta` to supply the
    rescue reference in genome-BAM mode).
  - `--junctions` (BED12) to hint spliced alignment instead of deriving junctions
    from `--annotation`.
- **`--dp-cache-cap-mb`**: cap (MB) on the mapper's retained per-thread DP scratch
  buffer; bounds peak RSS at high thread counts (default 128 MB; `0` = unbounded).
- **`--score-prob-denom`**: tune the score→probability denominator used to weight a
  read's alignments in the EM (transcriptome mode only).
- Documentation: genome modes, soft-clip rescue, the rammap mapping backend, and a
  GitHub Actions workflow that deploys the docs to GitHub Pages.

### Changed

- Genome-mode peak RSS reduced (~20% at high thread counts) by capping the rammap
  per-thread DP scratch-cache; alignment output is bit-identical.
- A transcriptome `--alignments` BAM produced by an unrecognized aligner now yields
  a clear error message instead of a panic.
- `--single-cell` now explicitly requires `--alignments` and conflicts with
  `--reads`, `--genome`, and `--genome-alignments`. Single-cell quantification has
  always required a transcriptome-aligned BAM; previously combining it with a
  genome mode silently ran bulk projection instead. It now errors clearly.
- **Build no longer requires a C toolchain / clang / bindgen** (the `minimap2-sys`
  dependency is gone), simplifying installation and cross-compilation.

### Notes

- Prebuilt **minimap2 `.mmi` indices remain usable** with `--index`: rammap reads
  both its own RMMI format and legacy `.mmi` files, and oarfish recomputes the
  reference signature from the index when no oarfish digest footer is present.
- The `rammap` dependency is temporarily pinned to the published fork
  `rammap-core-temp` 1.2.0, which carries the DP-cache cap (upstream
  [jwanglab/rammap#9](https://github.com/jwanglab/rammap/pull/9)). It will be
  switched back to `rammap-core` once that change is released upstream.
