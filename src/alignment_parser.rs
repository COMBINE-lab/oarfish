use crate::prog_opts::ProjProbSource;
use crate::util::constants::EMPTY_READ_NAME;
use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use crate::util::projection::projected_to_records;
use bramble_rs::g2t::G2TTree;
use bramble_rs::{GenomicAlignment, ProjectionConfig, ProjectionContext, project_group_with};
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::header::record::value::map::tag;
use noodles_sam::{Header, alignment::RecordBuf};
use num_format::{Locale, ToFormattedString};
use std::io;
use std::path::Path;
use swapvec::SwapVec;
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
        let so = tag::Other::try_from([b'S', b'O'])?;
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
        let so = tag::Other::try_from([b'S', b'O'])?;
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
    const CB_TAG: [u8; 2] = [b'C', b'B'];
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

    // Parse the input alignemnt file, gathering the alignments aggregated
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

    // Parse the input alignemnt file, gathering the alignments aggregated
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

pub fn parse_alignments<R: io::BufRead>(
    store: &mut InMemoryAlignmentStore,
    name_vec: &mut Option<SwapVec<String>>,
    header: &Header,
    reader: &mut bam::io::Reader<R>,
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
    let mut prev_read = String::new();
    let mut num_unmapped = 0_u64;
    let mut records_for_read = vec![];

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

    // Adds the read name for the read corresponding to the provided alignment group
    // `recs`, **if** we are keeping read names for the purpose of reporting read
    // assignment probabilities.
    let mut add_read_name = |recs: &Vec<noodles_sam::alignment::record_buf::RecordBuf>| {
        if let Some(nvec) = name_vec {
            let first_aln = recs.first().expect("alignment group should be non-empty");
            let read_name = first_aln
                .name()
                .unwrap_or(bstr::BStr::new(EMPTY_READ_NAME))
                .to_string();
            nvec.push(read_name)
                .expect("cannot push name to read name vector");
        }
    };

    // Parse the input alignemnt file, gathering the alignments aggregated
    // by their source read. **Note**: this requires that we have a
    // name-sorted input bam file (currently, aligned against the transcriptome).
    //
    // *NOTE*: this had to be changed from `records` to `record_bufs` or
    // critical information was missing from the records. This happened when
    // moving to the new version of noodles. Track `https://github.com/zaeleus/noodles/issues/230`
    // to see if it's clear why this is the case
    for result in reader.record_bufs(header) {
        let record = result?;
        pb.inc(1);

        // unmapped reads don't contribute to quantification
        // but we track them.
        if record.flags().is_unmapped() {
            num_unmapped += 1;
            continue;
        }
        let record_copy = record.clone();
        if let Some(rname) = record.name() {
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
                    if store.add_group(txps, &mut records_for_read) {
                        add_read_name(&records_for_read);
                        if records_for_read.len() == 1 {
                            store.inc_unique_alignments();
                        }
                    }
                    records_for_read.clear();
                }
                // the new "prev_read" name is the current read name
                // so it becomes the first on the new alignment range
                // vector.
                prev_read = rstring;
                if rg_num < check_order_thresh {
                    if !read_name_map.insert(prev_read.clone()) {
                        error!(
                            "It appears that the input BAM file is not name-collated. oarfish is not designed to process coordinate sorted BAM files."
                        );
                        anyhow::bail!(
                            "You appear to have provided a coordinate-sorted BAM, but oarfish does not support processing these.\n\
                                    You should provide a BAM file collated by record name (which is the \"natural\" minimap2 order).\n\
                                    Alignment records for the same read {} were observed twice in a non-contiguous block.",
                            &prev_read
                        );
                    }
                    rg_num += 1;
                }
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            }
        }
    }
    // if we end with a non-empty alignment range vector, then
    // add that group.
    if !records_for_read.is_empty() {
        // if we are using read names and we added the group here
        if store.add_group(txps, &mut records_for_read) {
            add_read_name(&records_for_read);
            if records_for_read.len() == 1 {
                store.inc_unique_alignments();
            }
        }
        records_for_read.clear();
    }

    pb.finish_with_message("Finished processing alignments.");

    info!(
        "the alignment file contained {} unmapped read records.",
        num_unmapped.to_formatted_string(&Locale::en)
    );

    Ok(())
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
        xs_strand: char_tag(rec, [b'X', b'S']),
        ts_strand: char_tag(rec, [b't', b's']),
        hit_index: int_tag(rec, [b'H', b'I']).unwrap_or(0) as i32,
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
            src_scores.push(int_tag(r, [b'A', b'S']).unwrap_or(0) as i32);
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
    name_vec: &mut Option<SwapVec<String>>,
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
                        &prev_read
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
