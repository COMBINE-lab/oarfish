use crate::util::constants::EMPTY_READ_NAME;
use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use noodles_bam as bam;
use noodles_sam::header::record::value::map::tag;
use noodles_sam::{alignment::RecordBuf, Header};
use num_format::{Locale, ToFormattedString};
use std::io;
use std::path::Path;
use swapvec::SwapVec;
use tracing::{error, info};

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

    // sort order tag
    let so = tag::Other::try_from([b'S', b'O'])?;
    let so_type = header
        .header()
        .expect("has inner header")
        .other_fields()
        .get(&so)
        .expect("BAM file should have @SO field");

    if so_type == "coordinate" {
        error!("oarfish is not designed to process coordinate sorted BAM files.");
        anyhow::bail!(
        "You provided a coordinate-sorted BAM, but oarfish does not support processing these.
         You should provide a BAM file collated by record name (which is the \"natural\" minimap2 order)."
        );
    }

    let mut saw_minimap2 = false;
    let mut progs = vec![];
    // explicitly check that alignment was done with a supported
    // aligner (right now, just minimap2).
    for (prog, _pmap) in header.programs().roots() {
        if prog == "minimap2" {
            saw_minimap2 = true;
            break;
        } else {
            progs.push(prog);
        }
    }
    assert!(
        saw_minimap2,
        "Currently, only minimap2 is supported as an aligner. The bam file listed {:?}.",
        progs
    );
    info!("saw minimap2 as a program in the header; proceeding.");
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
                        error!("It appears that the input BAM file is not name-collated. oarfish is not designed to process coordinate sorted BAM files.");
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
