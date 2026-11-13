use crate::prog_opts::ProjProbSource;
use crate::util::bgzf_records::RawRecord;
use crate::util::constants::EMPTY_READ_NAME;
use crate::util::name_store::NameStore;
use crate::util::oarfish_types::{
    AlnInfo, AlnRecordLike, InMemoryAlignmentStore, TranscriptInfo, add_coverage_adds,
    prefetch_read,
};
use crate::util::projection::projected_to_records;
use crate::util::read_ends::EndEvidence;
use bramble_rs::g2t::G2TTree;
use bramble_rs::{GenomicAlignment, ProjectionConfig, ProjectionContext, project_group_with};
use noodles_bam as bam;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::header::record::value::map::tag;
use noodles_sam::{Header, alignment::RecordBuf};
use num_format::{Locale, ToFormattedString};
use std::io;
use std::path::Path;
use tracing::{error, info};

const KNOWN_MAPPERS: [&str; 4] = ["minimap2", "pbmm2", "bramble", "rammap"];

pub fn read_and_verify_header<R: io::BufRead>(
    reader: &mut bam::io::Reader<R>,
    aln_file: &Path,
) -> anyhow::Result<Header> {
    // read the bam file header, print out some basic info
    let header = reader.read_header()?;
    info!(
        "read header from BAM file {}, contains {} reference sequences.",
        aln_file.display(),
        header
            .reference_sequences()
            .len()
            .to_formatted_string(&Locale::en)
    );

    // if we have an inner header field (@HD) then check for
    // the sort-order tag ane ensure it is *not* "coordinate".
    if let Some(inner_header) = header.header() {
        let so = tag::Other::try_from(*b"SO")?;
        let so_type_opt = inner_header.other_fields().get(&so);

        // we have an SO flag, ensure it's not "coordinate"
        if let Some(so_type) = so_type_opt {
            if so_type == "coordinate" {
                error!("oarfish is not designed to process coordinate sorted BAM files.");
                anyhow::bail!(
                "You provided a coordinate-sorted BAM, but oarfish does not support processing these.
                    You should provide a BAM file collated by record name (which is the \"natural\" minimap2 order)."
            );
            }
        } else {
            // we had no SO flag
            info!(
                "BAM file had \"inner header\", but that header did not have the @SO field; cannot determine sort order from header alone."
            );
        }
    } else {
        // we had no inner header
        info!(
            "BAM file did not have \"inner header\", so cannot determine the sort order from the header alone."
        );
    }

    let mut matched_prog = None;
    let mut progs = vec![];
    // Check that the alignment was produced by an aligner we have validated
    // (right now, minimap2, pbmm2, bramble, and rammap).
    for (prog, _pmap) in header.programs().roots() {
        if KNOWN_MAPPERS.iter().any(|known| prog == *known) {
            matched_prog = Some(prog);
            break;
        } else {
            progs.push(prog);
        }
    }
    match matched_prog {
        Some(prog) => {
            info!(
                "saw supported mapper {} as a program in the header; proceeding.",
                prog
            );
            Ok(header)
        }
        None => {
            // A previously fatal `assert!` here turned an unrecognized aligner
            // into a panic. Fail gracefully with an actionable message instead.
            anyhow::bail!(
                "Could not find a validated aligner in the BAM @PG header.\n\
                 oarfish's transcriptome alignment mode currently recognizes [{}], but this BAM \
                 listed {:?}.\n\
                 Alignments from other tools may still work (oarfish only requires the per-record \
                 `AS` alignment-score tag), but they have not been validated. If you would like \
                 support for your aligner, please open an issue at \
                 https://github.com/COMBINE-lab/oarfish/issues.",
                KNOWN_MAPPERS.join(", "),
                progs
            )
        }
    }
}

/// Read and verify the header of a *genome*-aligned BAM (genome-BAM projection
/// mode). Like [`read_and_verify_header`] this rejects coordinate-sorted input
/// (projection groups alignments by read name, so the BAM must be name-collated),
/// but it does **not** restrict the aligner: spliced genome alignments can come
/// from minimap2, pbmm2, STAR, and others. The program records are logged for
/// provenance only.
pub fn read_and_verify_genome_header<R: io::BufRead>(
    reader: &mut bam::io::Reader<R>,
    aln_file: &Path,
) -> anyhow::Result<Header> {
    let header = reader.read_header()?;
    info!(
        "read header from genome BAM file {}, contains {} reference sequences.",
        aln_file.display(),
        header
            .reference_sequences()
            .len()
            .to_formatted_string(&Locale::en)
    );

    if let Some(inner_header) = header.header() {
        let so = tag::Other::try_from(*b"SO")?;
        if let Some(so_type) = inner_header.other_fields().get(&so)
            && so_type == "coordinate"
        {
            error!("oarfish cannot project a coordinate-sorted genome BAM.");
            anyhow::bail!(
                "You provided a coordinate-sorted genome BAM, but genome-projection mode requires\n\
                 a BAM collated by read name so that all alignments of a read are contiguous.\n\
                 Please collate it first, e.g. `samtools collate -o collated.bam input.bam`."
            );
        }
    }

    let progs: Vec<_> = header.programs().roots().map(|(prog, _)| prog).collect();
    info!(
        "genome BAM @PG program(s): {:?}; proceeding with projection.",
        progs
    );
    Ok(header)
}

pub enum NextAction {
    SkipUnmapped,
    ProcessSameBarcode,
    NewBarcode,
    EndOfFile,
}

#[inline(always)]
fn is_same_barcode(rec: &RecordBuf, current_barcode: &[u8]) -> anyhow::Result<bool> {
    const CB_TAG: [u8; 2] = *b"CB";
    let same_barcode = match rec.data().get(&CB_TAG) {
        None => anyhow::bail!("could not get CB tag value"),
        Some(v) => match v {
            noodles_sam::alignment::record_buf::data::field::Value::String(x) => {
                x.as_slice() == current_barcode
            }
            _ => anyhow::bail!("CB tag value had unexpected type!"),
        },
    };
    Ok(same_barcode)
}

/// Takes a collection of [RecordBuf]s, a mutable reference to an [InMemoryAlignmentStore]
/// as well as the relevant [TranscriptInfo].
///
/// This function first sorts the `records` by the read name (ensuring that any primary alignment
/// of the read comes first), so that reads are collated by name.  Then, it processes in turn the
/// alignments for each input read, filtering them according to the filters attached to the
/// `store`.  Subsequently, the alignments are summarized in the `store`. If all `records` are
/// processed successfully, [Ok]`()` is returned, otherwise the relevant [anyhow::Error] is
/// returned.
pub fn sort_and_parse_barcode_records(
    records: &mut Vec<RecordBuf>,
    store: &mut InMemoryAlignmentStore,
    txps: &mut [TranscriptInfo],
    records_for_read: &mut Vec<RecordBuf>,
) -> anyhow::Result<()> {
    records_for_read.clear();
    let mut prev_read = String::new();

    // first sort records by read name
    records.sort_unstable_by(|x, y| match x.name().cmp(&y.name()) {
        std::cmp::Ordering::Equal => {
            match (x.flags().is_secondary(), y.flags().is_secondary()) {
                (false, true) => std::cmp::Ordering::Less,
                (true, false) => std::cmp::Ordering::Greater,
                (true, true) => std::cmp::Ordering::Equal,
                // this one shouldn't happen
                (false, false) => std::cmp::Ordering::Equal,
            }
        }
        x => x,
    });

    // Parse the input alignment file, gathering the alignments aggregated
    // by their source read. **Note**: this requires that we have a
    // name-sorted input bam file (currently, aligned against the transcriptome).
    //
    // *NOTE*: this had to be changed from `records` to `record_bufs` or
    // critical information was missing from the records. This happened when
    // moving to the new version of noodles. Track `https://github.com/zaeleus/noodles/issues/230`
    // to see if it's clear why this is the case
    for record in records {
        if let Some(rname) = record.name() {
            let record_copy = record.clone();
            let rstring: String = String::from_utf8_lossy(rname.as_ref()).into_owned();
            // if this is an alignment for the same read, then
            // push it onto our temporary vector.
            if prev_read == rstring {
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            } else {
                // otherwise, record the alignment range for the
                // previous read record.
                if !prev_read.is_empty() {
                    store.add_group(txps, records_for_read);
                    if records_for_read.len() == 1 {
                        store.inc_unique_alignments();
                    }
                    records_for_read.clear();
                }
                // the new "prev_read" name is the current read name
                // so it becomes the first on the new alignment range
                // vector.
                prev_read = rstring;
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            }
        }
    }
    // if we end with a non-empty alignment range vector, then
    // add that group.
    if !records_for_read.is_empty() {
        store.add_group(txps, records_for_read);
        if records_for_read.len() == 1 {
            store.inc_unique_alignments();
        }
        records_for_read.clear();
    }
    Ok(())
}

#[inline(always)]
pub fn parse_alignments_for_barcode<R: io::BufRead>(
    iter: &mut core::iter::Peekable<noodles_bam::io::reader::RecordBufs<R>>,
    current_cb: &[u8],
) -> anyhow::Result<Vec<noodles_sam::alignment::record_buf::RecordBuf>> {
    //records_for_read.clear();
    let mut records_for_barcode =
        Vec::<noodles_sam::alignment::record_buf::RecordBuf>::with_capacity(2048);
    let mut _records_processed = 0_usize;

    // Parse the input alignment file, gathering the alignments aggregated
    // by their source read. **Note**: this requires that we have a
    // name-sorted input bam file (currently, aligned against the transcriptome).
    //
    // *NOTE*: this had to be changed from `records` to `record_bufs` or
    // critical information was missing from the records. This happened when
    // moving to the new version of noodles. Track `https://github.com/zaeleus/noodles/issues/230`
    // to see if it's clear why this is the case

    loop {
        let action = if let Some(result) = iter.peek() {
            if let Ok(record) = result {
                // unmapped reads don't contribute to quantification
                // but we track them.
                if record.flags().is_unmapped() {
                    _records_processed += 1;
                    NextAction::SkipUnmapped
                } else {
                    let same_barcode = is_same_barcode(record, current_cb)?;
                    if !same_barcode {
                        NextAction::NewBarcode
                    } else {
                        NextAction::ProcessSameBarcode
                    }
                }
            } else {
                anyhow::bail!("error parsing record");
            }
        } else {
            NextAction::EndOfFile
        };

        match action {
            NextAction::SkipUnmapped => {
                iter.next().unwrap()?;
            }
            NextAction::ProcessSameBarcode => {
                let record = iter.next().unwrap()?;
                records_for_barcode.push(record.clone());
            }
            NextAction::NewBarcode | NextAction::EndOfFile => {
                break;
            }
        };
    }
    Ok(records_for_barcode)
}

pub fn parse_alignments(
    store: &mut InMemoryAlignmentStore,
    name_vec: &mut Option<NameStore>,
    _header: &Header,
    input: BamRecordStream,
    txps: &mut [TranscriptInfo],
    check_order_thresh: usize,
    quiet: bool,
) -> anyhow::Result<()> {
    //use blart::TreeMap;
    use rustc_hash::FxHashSet;

    let mut read_name_map = FxHashSet::default();
    read_name_map.reserve(check_order_thresh);
    let mut rg_num = 0_usize;

    // we'll need these to keep track of which alignments belong
    // to which reads.
    let mut prev_read: Vec<u8> = Vec::new();
    let mut num_unmapped = 0_u64;

    let pb = if quiet {
        indicatif::ProgressBar::hidden()
    } else {
        indicatif::ProgressBar::new_spinner().with_message("Number of alignments processed")
    };

    pb.set_style(
        indicatif::ProgressStyle::with_template(
            "[{elapsed_precise}] {spinner:4.green/blue} {msg} {human_pos:>12}",
        )
        .unwrap()
        .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
    );
    pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

    // Parse the input alignment file, gathering the alignments aggregated
    // by their source read. **Note**: this requires that we have a
    // name-sorted input bam file (currently, aligned against the transcriptome).
    //
    // The work is pipelined (see `util::bgzf_records`): the workers that inflate the BGZF
    // blocks also read each block's records (in place, as `RawRecord`s, whose fields are read
    // from the inflated bytes on demand) and take from each the few fields the filters and the store
    // use (a `ParsedAln`, each read once), grouping a block's consecutive records of one read
    // name into runs. This thread joins the blocks' runs in input order into reads (a read cut
    // by a block end is joined here), checks the collation, and groups the reads into batches
    // of `PARSE_BATCH_READS` reads; a second thread filters the batches' reads and adds them
    // to the store, in input order (see `add_read_batch`), so the store, the coverage
    // profiles, the discard counts and the read names come out exactly as a serial pass
    // builds them. Only the compact fields cross threads, never the records. (Filtering each
    // batch in parallel on the rayon pool instead was slower than one adder thread.)
    //
    // Under the coverage model the coverage adds (`add_unit_interval`), which would bound the
    // adder, run on a third thread that owns the transcripts: the adder filters with a copy of
    // the transcripts' lengths (all the filters read of them) and sends each batch's kept
    // alignments, in order, to that thread, which makes the adds in the order the adder would
    // have (see `add_coverage_adds`).
    let header = store.aln_header;
    let tail_evidence = store.filter_opts.tail_neutral.is_some();
    let keep_names = name_vec.is_some();
    let defer_coverage = store.filter_opts.model_coverage;
    let lens: Vec<usize> = if defer_coverage {
        txps.iter().map(|t| t.len.get()).collect()
    } else {
        Vec::new()
    };
    let (batch_tx, batch_rx) = std::sync::mpsc::sync_channel::<ReadBatch>(2);
    let (recycle_tx, recycle_rx) = std::sync::mpsc::channel::<ReadBatch>();
    let (adds_tx, adds_rx) = std::sync::mpsc::sync_channel::<Vec<AlnInfo>>(2);
    let (adds_recycle_tx, adds_recycle_rx) = std::sync::mpsc::channel::<Vec<AlnInfo>>();
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let (mut adder_txps, coverage) = if defer_coverage {
            let coverage = scope.spawn(move || {
                for mut adds in adds_rx {
                    add_coverage_adds(txps, &adds);
                    adds.clear();
                    // the adder may have stopped; the list is then dropped
                    let _ = adds_recycle_tx.send(adds);
                }
            });
            let adder_txps = AdderTxps::Deferred {
                lens: &lens,
                adds: Vec::new(),
            };
            (adder_txps, Some(coverage))
        } else {
            (AdderTxps::Direct(txps), None)
        };
        let adder = scope.spawn(move || {
            let mut group = Vec::new();
            let mut coverage_gone = false;
            for mut batch in batch_rx {
                add_read_batch(store, name_vec, &mut adder_txps, &mut batch, &mut group);
                batch.clear();
                // the reader may have stopped early (on an error); the batch is then dropped
                let _ = recycle_tx.send(batch);
                if let AdderTxps::Deferred { adds, .. } = &mut adder_txps {
                    if coverage_gone {
                        adds.clear();
                    } else if !adds.is_empty() {
                        let next = adds_recycle_rx.try_recv().unwrap_or_default();
                        // the coverage thread has gone only if it panicked; joining it says why
                        coverage_gone = adds_tx.send(std::mem::replace(adds, next)).is_err();
                    }
                }
            }
            // dropping `adds_tx` here lets the coverage thread finish the queued lists and stop
        });

        let mut n_records = 0_u64;
        let mut next_tick = 0_u64;
        let mut batch = ReadBatch::default();
        // where the current read's alignments start in `batch.alns`
        let mut read_start = 0_usize;

        let read_result = crate::util::bgzf_records::for_each_block(
            input.file,
            input.first_block_skip,
            input.workers,
            |record: &RawRecord<'_>, out: &mut BlockReads| out.push(record, header, tail_evidence),
            |out: &mut BlockReads| -> anyhow::Result<bool> {
                n_records += out.n_records;
                num_unmapped += out.n_unmapped;
                if n_records >= next_tick {
                    pb.set_position(n_records);
                    next_tick = n_records + 65536;
                }
                let mut name_start = 0;
                let mut alns_start = 0;
                for run in &out.runs {
                    let rname = &out.names[name_start..run.name_end];
                    name_start = run.name_end;
                    // if this is an alignment for the same read, then
                    // push it onto our temporary vector.
                    if prev_read.as_slice() != rname {
                        // otherwise, close the previous read's group (a read none of whose
                        // records has a reference adds nothing, as `add_group` would).
                        if batch.alns.len() > read_start {
                            batch.close_read(keep_names, &prev_read);
                            if batch.len() == PARSE_BATCH_READS {
                                let next = recycle_rx.try_recv().unwrap_or_default();
                                if batch_tx.send(std::mem::replace(&mut batch, next)).is_err() {
                                    // the adder has gone (it panicked); joining it reports why
                                    return Ok(false);
                                }
                            }
                            read_start = batch.alns.len();
                        }
                        // the new "prev_read" name is the current read name
                        // so it becomes the first on the new alignment range
                        // vector.
                        prev_read.clear();
                        prev_read.extend_from_slice(rname);
                        if rg_num < check_order_thresh {
                            if !read_name_map.insert(prev_read.clone()) {
                                error!(
                                    "It appears that the input BAM file is not name-collated. oarfish is not designed to process coordinate sorted BAM files."
                                );
                                anyhow::bail!(
                                    "You appear to have provided a coordinate-sorted BAM, but oarfish does not support processing these.\n\
                                        You should provide a BAM file collated by record name (which is the \"natural\" minimap2 order).\n\
                                        Alignment records for the same read {} were observed twice in a non-contiguous block.",
                                    String::from_utf8_lossy(&prev_read)
                                );
                            }
                            rg_num += 1;
                        }
                    }
                    batch
                        .alns
                        .extend_from_slice(&out.alns[alns_start..run.alns_end]);
                    alns_start = run.alns_end;
                }
                Ok(true)
            },
        );
        let read_result = read_result.map(|()| {
            pb.set_position(n_records);
            // if we end with a non-empty alignment range vector, then
            // add that group.
            if batch.alns.len() > read_start {
                batch.close_read(keep_names, &prev_read);
            }
            if batch.len() > 0 {
                let _ = batch_tx.send(batch);
            }
        });
        // closing the channel lets the adder finish the queued batches and stop
        drop(batch_tx);
        if let Err(panic) = adder.join() {
            std::panic::resume_unwind(panic);
        }
        if let Some(coverage) = coverage
            && let Err(panic) = coverage.join()
        {
            std::panic::resume_unwind(panic);
        }
        read_result
    })?;

    pb.finish_with_message("Finished processing alignments.");

    info!(
        "the alignment file contained {} unmapped read records.",
        num_unmapped.to_formatted_string(&Locale::en)
    );

    Ok(())
}

/// The BAM record stream [`parse_alignments`] reads: a BGZF file positioned at the start of
/// the block holding the first record, how many inflated bytes of that block precede it (the
/// end of the header), and how many workers inflate and read the blocks.
pub struct BamRecordStream {
    pub file: std::fs::File,
    pub first_block_skip: usize,
    pub workers: usize,
}

/// The records of one BGZF block, as a worker of [`parse_alignments`] reads them: the fields
/// of each record with a reference (`alns`), grouped into runs of consecutive records of one
/// read name, with the names of the runs.
#[derive(Default)]
struct BlockReads {
    n_records: u64,
    n_unmapped: u64,
    alns: Vec<ParsedAln>,
    runs: Vec<NameRun>,
    names: Vec<u8>,
    /// where the last run's name starts in `names`
    last_name_start: usize,
}

/// Consecutive records of one read name in a block: the name is `names[..name_end]` after
/// the previous run's name, the alignments `alns[..alns_end]` after the previous run's.
#[derive(Clone, Copy)]
struct NameRun {
    name_end: usize,
    alns_end: usize,
}

impl BlockReads {
    /// Takes one record, as the serial parse loop did: unmapped records are counted and
    /// skipped, records without a name skipped, and a record joins its read's alignments only
    /// when it has a reference.
    #[inline]
    fn push(&mut self, record: &RawRecord<'_>, header: &Header, tail_evidence: bool) {
        self.n_records += 1;
        // unmapped reads don't contribute to quantification
        // but we track them.
        if Flags::from(record.flags()).is_unmapped() {
            self.n_unmapped += 1;
            return;
        }
        let Some(rname) = record.name() else {
            return;
        };
        let same = self
            .runs
            .last()
            .is_some_and(|r| &self.names[self.last_name_start..r.name_end] == rname);
        if !same {
            self.last_name_start = self.names.len();
            self.names.extend_from_slice(rname);
            self.runs.push(NameRun {
                name_end: self.names.len(),
                alns_end: self.alns.len(),
            });
        }
        if record.reference_sequence_id().is_some() {
            self.alns
                .push(ParsedAln::new(record, header, tail_evidence));
            self.runs.last_mut().unwrap().alns_end = self.alns.len();
        }
    }
}

impl crate::util::bgzf_records::BlockOutput for BlockReads {
    fn clear(&mut self) {
        self.n_records = 0;
        self.n_unmapped = 0;
        self.alns.clear();
        self.runs.clear();
        self.names.clear();
        self.last_name_start = 0;
    }
}

/// How many reads [`parse_alignments`] groups into one batch for the adder.
const PARSE_BATCH_READS: usize = 4096;

/// The fields of one BAM alignment record that the filters ([`AlignmentFilters::filter`])
/// and the store read, each taken from the record once. As an [`AlnRecordLike`] it answers
/// every query exactly as the record would.
///
/// [`AlignmentFilters::filter`]: crate::util::oarfish_types::AlignmentFilters::filter
#[derive(Clone, Copy, Debug)]
struct ParsedAln {
    ref_id: Option<usize>,
    /// zero-based start (`None`: the record has no valid position)
    start: Option<u32>,
    span: Option<usize>,
    score: Option<i64>,
    flags: Flags,
    seq_len: usize,
    /// read only when tail-neutral scores are on; otherwise [`EndEvidence::UNKNOWN`]
    ends: EndEvidence,
}

impl ParsedAln {
    /// Takes the fields from the record's bytes as the `noodles_bam::Record` backend of
    /// [`AlnRecordLike`] reads them (a negative reference id or position is no valid one).
    fn new(rec: &RawRecord<'_>, _header: &Header, tail_evidence: bool) -> Self {
        let cigar = rec.cigar();
        Self {
            ref_id: rec
                .reference_sequence_id()
                .and_then(|n| usize::try_from(n).ok()),
            start: rec.alignment_start().and_then(|n| u32::try_from(n).ok()),
            span: raw_aln_span(cigar),
            score: rec
                .int_tag(*b"AS")
                .expect("alignment record should have an AS tag")
                .expect("could not get value"),
            flags: Flags::from(rec.flags()),
            seq_len: rec.sequence_len(),
            ends: if tail_evidence {
                raw_end_evidence(rec, cigar)
            } else {
                EndEvidence::UNKNOWN
            },
        }
    }
}

/// The reference span of packed CIGAR ops (`len << 4 | kind`): M, D, N, = and X (kinds 0, 2,
/// 3, 7 and 8) consume the reference; `None` for a span of 0.
#[inline]
fn raw_aln_span(cigar: &[u8]) -> Option<usize> {
    const CONSUMES_REF: u32 = 1 << 0 | 1 << 2 | 1 << 3 | 1 << 7 | 1 << 8;
    let (ops, _) = cigar.as_chunks::<4>();
    let span: usize = ops
        .iter()
        .map(|b| u32::from_le_bytes(*b))
        .filter(|n| CONSUMES_REF >> (n & 0xf) & 1 == 1)
        .map(|n| (n >> 4) as usize)
        .sum();
    (span > 0).then_some(span)
}

/// What a record says about the read's ends (see [`crate::util::read_ends`]): its terminal
/// clips (kinds 4 and 5) and, when it holds SEQ, the A-run at SEQ's right end, which is the
/// molecule's tail whichever strand the read aligned to (SEQ is in reference orientation).
#[inline]
fn raw_end_evidence(rec: &RawRecord<'_>, ops: &[u8]) -> EndEvidence {
    // the base symbols of the 4-bit codes, as `noodles_bam` decodes them
    const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let clip = |op: Option<&[u8]>| {
        op.map(|b| u32::from_le_bytes(b.try_into().expect("4-byte CIGAR op")))
            .filter(|n| matches!(n & 0xf, 4 | 5))
            .map_or(0, |n| n >> 4)
    };
    let n = rec.sequence_len();
    EndEvidence {
        clip_left: Some(clip(ops.get(..4))),
        clip_right: Some(clip(ops.len().checked_sub(4).map(|i| &ops[i..]))),
        tail: (n > 0).then(|| {
            crate::util::read_ends::terminal_run(
                (0..n).rev().map(|i| BASES[rec.base_code(i) as usize]),
                b'A',
            )
        }),
    }
}

impl AlnRecordLike for ParsedAln {
    fn opt_sequence_len(&self) -> Option<usize> {
        Some(self.seq_len)
    }

    fn is_reverse_complemented(&self) -> bool {
        self.flags.is_reverse_complemented()
    }

    fn is_unmapped(&self) -> bool {
        self.flags.is_unmapped()
    }

    fn ref_id(&self, _header: &Header) -> anyhow::Result<usize> {
        self.ref_id
            .ok_or_else(|| anyhow::anyhow!("alignment record has no valid reference sequence id"))
    }

    fn aln_span(&self) -> Option<usize> {
        self.span
    }

    fn aln_score(&self) -> Option<i64> {
        self.score
    }

    fn aln_start(&self) -> u32 {
        self.start.expect("valid aln start")
    }

    fn aln_end(&self) -> u32 {
        self.aln_start() + self.span.unwrap_or(1) as u32
    }

    fn is_supp(&self) -> bool {
        self.flags.is_supplementary()
    }

    fn end_evidence(&self) -> EndEvidence {
        self.ends
    }

    fn name(&self) -> Option<String> {
        None
    }
}

/// A batch of reads for the adder: the alignments of each read in turn, where each read's
/// end, and (only when read names are kept) each read's name.
#[derive(Default)]
struct ReadBatch {
    alns: Vec<ParsedAln>,
    ends: Vec<usize>,
    names: Vec<String>,
}

impl ReadBatch {
    fn len(&self) -> usize {
        self.ends.len()
    }

    /// Ends the current read, whose alignments are the ones pushed since the last end.
    fn close_read(&mut self, keep_name: bool, name: &[u8]) {
        self.ends.push(self.alns.len());
        if keep_name {
            self.names.push(bstr::BStr::new(name).to_string());
        }
    }

    fn clear(&mut self) {
        self.alns.clear();
        self.ends.clear();
        self.names.clear();
    }
}

/// Filters the reads of `batch` and adds each to `store`, in batch order, exactly as
/// [`InMemoryAlignmentStore::add_group`] would: the discard counts, the coverage profiles, the
/// read names and the count of uniquely aligned reads all come out the same.
fn add_read_batch(
    store: &mut InMemoryAlignmentStore,
    name_vec: &mut Option<NameStore>,
    txps: &mut AdderTxps<'_>,
    batch: &mut ReadBatch,
    group: &mut Vec<ParsedAln>,
) {
    // Under the coverage model the filter reads each alignment's transcript length from the
    // compact `lens` (random reads into it, as the adds' were into the transcripts); prefetch
    // those `PREFETCH_LEN_AHEAD` alignments ahead. A prefetch changes no value.
    const PREFETCH_LEN_AHEAD: usize = 16;
    let n_alns = batch.alns.len();
    let mut len_next = 0_usize;
    let mut start = 0;
    for (i, &end) in batch.ends.iter().enumerate() {
        group.clear();
        group.extend_from_slice(&batch.alns[start..end]);
        start = end;
        let added = match txps {
            AdderTxps::Direct(txps) => store.add_group(txps, group),
            AdderTxps::Deferred { lens, adds } => {
                let len_to = (end + PREFETCH_LEN_AHEAD).min(n_alns);
                for a in &batch.alns[len_next.min(len_to)..len_to] {
                    if let Some(tid) = a.ref_id {
                        prefetch_read(lens.as_ptr().wrapping_add(tid));
                    }
                }
                len_next = len_next.max(len_to);
                store.add_group_deferring_coverage(|tid| lens[tid], group, adds)
            }
        };
        if added {
            // if we are keeping read names, for the purpose of reporting read assignment
            // probabilities
            if let Some(nvec) = name_vec {
                nvec.push(std::mem::take(&mut batch.names[i]))
                    .expect("cannot push name to read name vector");
            }
            if group.len() == 1 {
                store.inc_unique_alignments();
            }
        }
    }
}

/// How the adder of [`parse_alignments`] reaches the transcripts: directly (it filters with
/// their lengths and makes the coverage adds itself), or, under the coverage model, through
/// a copy of their lengths, listing the coverage adds in `adds` for the coverage thread.
enum AdderTxps<'a> {
    Direct(&'a mut [TranscriptInfo]),
    Deferred {
        lens: &'a [usize],
        adds: Vec<AlnInfo>,
    },
}

/// Map a noodles CIGAR [`CigarKind`] to the SAM op code byte expected by
/// bramble's [`GenomicAlignment::cigar`] (`0`=M, `1`=I, `2`=D, `3`=N, `4`=S,
/// `5`=H, `6`=P, `7`==, `8`=X).
#[inline]
fn cigar_kind_to_sam_u8(kind: CigarKind) -> u8 {
    match kind {
        CigarKind::Match => 0,
        CigarKind::Insertion => 1,
        CigarKind::Deletion => 2,
        CigarKind::Skip => 3,
        CigarKind::SoftClip => 4,
        CigarKind::HardClip => 5,
        CigarKind::Pad => 6,
        CigarKind::SequenceMatch => 7,
        CigarKind::SequenceMismatch => 8,
    }
}

/// Read a single-character (`type A`) tag (e.g. `ts`, `XS`) as a `char`.
fn char_tag(rec: &RecordBuf, tag: [u8; 2]) -> Option<char> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&tag) {
        Some(Value::Character(c)) => Some(*c as char),
        _ => None,
    }
}

/// Read an integer tag (e.g. `HI`) as an `i32`, tolerating any integer width.
fn int_tag(rec: &RecordBuf, tag: [u8; 2]) -> Option<i64> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&tag) {
        Some(Value::Int8(v)) => Some(*v as i64),
        Some(Value::UInt8(v)) => Some(*v as i64),
        Some(Value::Int16(v)) => Some(*v as i64),
        Some(Value::UInt16(v)) => Some(*v as i64),
        Some(Value::Int32(v)) => Some(*v as i64),
        Some(Value::UInt32(v)) => Some(*v as i64),
        _ => None,
    }
}

/// Convert a noodles [`RecordBuf`] (one record of a genome-aligned BAM) into a
/// bramble [`GenomicAlignment`]. Returns `None` for unmapped records or those
/// lacking a reference id / start. `query_name` is supplied by the caller (all
/// records in a group share it).
fn record_buf_to_genomic_alignment(rec: &RecordBuf, query_name: &str) -> Option<GenomicAlignment> {
    let flags = rec.flags();
    if flags.is_unmapped() {
        return None;
    }
    let ref_id = rec.reference_sequence_id()? as i32;
    let ref_start = rec.alignment_start()?.get() as i64; // 1-based (SAM POS)

    let cigar: Vec<(u32, u8)> = rec
        .cigar()
        .as_ref()
        .iter()
        .map(|op| (op.len() as u32, cigar_kind_to_sam_u8(op.kind())))
        .collect();

    let seq = rec.sequence();
    let read_len = seq.len();
    let sequence = if read_len == 0 {
        None
    } else {
        Some(seq.as_ref().to_vec())
    };

    Some(GenomicAlignment {
        query_name: query_name.to_string(),
        ref_id,
        ref_start,
        is_reverse: flags.is_reverse_complemented(),
        cigar,
        sequence,
        is_paired: flags.is_segmented(),
        is_first_in_pair: flags.is_first_segment(),
        xs_strand: char_tag(rec, *b"XS"),
        ts_strand: char_tag(rec, *b"ts"),
        hit_index: int_tag(rec, *b"HI").unwrap_or(0) as i32,
        mate_ref_id: rec.mate_reference_sequence_id().map(|x| x as i32),
        mate_ref_start: rec.mate_alignment_start().map(|p| p.get() as i64),
        mate_is_unmapped: flags.is_mate_unmapped(),
        read_len,
    })
}

/// Project and store one read's genomic alignment group.
///
/// Converts the group's records to [`GenomicAlignment`]s, projects them onto the
/// transcriptome with bramble, and adds the (filtered) projected alignments to
/// `store`. Returns `true` if any projected alignment was retained.
#[allow(clippy::too_many_arguments)]
fn project_and_add_group(
    store: &mut InMemoryAlignmentStore,
    txps: &mut [TranscriptInfo],
    g2t: &G2TTree,
    proj_config: &ProjectionConfig,
    beta: f32,
    prob_source: ProjProbSource,
    pctx: &mut ProjectionContext,
    query_name: &str,
    records_for_read: &[RecordBuf],
) -> bool {
    // build GenomicAlignments and a parallel vector of their alignment scores
    // (the BAM `AS` tag, if present) used by the score/combined prob sources.
    let mut alns: Vec<GenomicAlignment> = Vec::with_capacity(records_for_read.len());
    let mut src_scores: Vec<i32> = Vec::with_capacity(records_for_read.len());
    for r in records_for_read {
        if let Some(ga) = record_buf_to_genomic_alignment(r, query_name) {
            alns.push(ga);
            src_scores.push(int_tag(r, *b"AS").unwrap_or(0) as i32);
        }
    }
    if alns.is_empty() {
        return false;
    }
    // the query length is the same across a read's records; take the first with
    // a recorded sequence/length.
    let read_len = alns
        .iter()
        .map(|a| a.read_len)
        .find(|&l| l > 0)
        .unwrap_or(0);

    let projected = project_group_with(&alns, g2t, proj_config, pctx);
    if projected.is_empty() {
        return false;
    }
    let recs = projected_to_records(&projected, &src_scores);
    store.add_projected_group(txps, &recs, read_len, beta, prob_source)
}

/// Parse a name-collated, genome-aligned BAM, projecting each read's alignments
/// onto the transcriptome (via bramble) and accumulating them in `store`.
///
/// This mirrors [`parse_alignments`] (same read-name grouping and collation
/// check), but each group is projected with bramble instead of filtered against
/// the transcriptome directly. `store` must have been created over the
/// *transcriptome* header.
#[allow(clippy::too_many_arguments)]
pub fn parse_genome_alignments<R: io::BufRead>(
    store: &mut InMemoryAlignmentStore,
    name_vec: &mut Option<NameStore>,
    genome_header: &Header,
    g2t: &G2TTree,
    proj_config: &ProjectionConfig,
    beta: f32,
    prob_source: ProjProbSource,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    check_order_thresh: usize,
    quiet: bool,
) -> anyhow::Result<()> {
    use rustc_hash::FxHashSet;

    let mut read_name_map = FxHashSet::default();
    read_name_map.reserve(check_order_thresh);
    let mut rg_num = 0_usize;

    let mut prev_read = String::new();
    let mut num_unmapped = 0_u64;
    let mut records_for_read: Vec<RecordBuf> = vec![];
    // reused across all read groups (avoids per-read allocation of bramble's
    // projection scratch / ksw2 aligner).
    let mut pctx = ProjectionContext::new();

    let pb = if quiet {
        indicatif::ProgressBar::hidden()
    } else {
        indicatif::ProgressBar::new_spinner().with_message("Number of alignments processed")
    };
    pb.set_style(
        indicatif::ProgressStyle::with_template(
            "[{elapsed_precise}] {spinner:4.green/blue} {msg} {human_pos:>12}",
        )
        .unwrap()
        .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
    );
    pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

    let mut add_read_name = |recs: &[RecordBuf]| {
        if let Some(nvec) = name_vec {
            let read_name = recs
                .first()
                .and_then(|r| r.name())
                .map(|n| n.to_string())
                .unwrap_or_else(|| EMPTY_READ_NAME.to_string());
            nvec.push(read_name)
                .expect("cannot push name to read name vector");
        }
    };

    for result in reader.record_bufs(genome_header) {
        let record = result?;
        pb.inc(1);

        if record.flags().is_unmapped() {
            num_unmapped += 1;
            continue;
        }

        let Some(rname) = record.name() else { continue };
        let rstring: String = String::from_utf8_lossy(rname.as_ref()).into_owned();

        if prev_read == rstring {
            if record.reference_sequence_id().is_some() {
                records_for_read.push(record.clone());
            }
        } else {
            if !prev_read.is_empty()
                && project_and_add_group(
                    store,
                    txps,
                    g2t,
                    proj_config,
                    beta,
                    prob_source,
                    &mut pctx,
                    &prev_read,
                    &records_for_read,
                )
            {
                add_read_name(&records_for_read);
                if records_for_read.len() == 1 {
                    store.inc_unique_alignments();
                }
            }
            records_for_read.clear();

            prev_read = rstring;
            if rg_num < check_order_thresh {
                if !read_name_map.insert(prev_read.clone()) {
                    error!("It appears that the input genome BAM file is not name-collated.");
                    anyhow::bail!(
                        "You appear to have provided a coordinate-sorted genome BAM, but genome-projection\n\
                         mode requires a BAM collated by read name. Alignment records for read {} were\n\
                         observed twice in a non-contiguous block. Try `samtools collate`.",
                        prev_read
                    );
                }
                rg_num += 1;
            }
            if record.reference_sequence_id().is_some() {
                records_for_read.push(record.clone());
            }
        }
    }

    if !records_for_read.is_empty()
        && project_and_add_group(
            store,
            txps,
            g2t,
            proj_config,
            beta,
            prob_source,
            &mut pctx,
            &prev_read,
            &records_for_read,
        )
    {
        add_read_name(&records_for_read);
        if records_for_read.len() == 1 {
            store.inc_unique_alignments();
        }
    }

    pb.finish_with_message("Finished processing alignments.");
    info!(
        "the genome alignment file contained {} unmapped read records.",
        num_unmapped.to_formatted_string(&Locale::en)
    );

    Ok(())
}
