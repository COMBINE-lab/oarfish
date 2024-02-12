use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use noodles_bam as bam;
use noodles_sam::Header;
use num_format::{Locale, ToFormattedString};
use std::io;
use tracing::{info, trace};

pub fn parse_alignments<R: io::BufRead>(
    store: &mut InMemoryAlignmentStore,
    header: &Header,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
) -> anyhow::Result<()> {
    // we'll need these to keep track of which alignments belong
    // to which reads.
    let mut prev_read = String::new();
    let mut num_unmapped = 0_u64;
    let mut records_for_read = vec![];

    // Parse the input alignemnt file, gathering the alignments aggregated
    // by their source read. **Note**: this requires that we have a
    // name-sorted input bam file (currently, aligned against the transcriptome).
    //
    // *NOTE*: this had to be changed from `records` to `record_bufs` or
    // critical information was missing from the records. This happened when
    // moving to the new version of noodles. Track `https://github.com/zaeleus/noodles/issues/230`
    // to see if it's clear why this is the case
    for (i, result) in reader.record_bufs(header).enumerate() {
        let record = result?;

        if i % 100_000 == 1 {
            trace!("processed {i} alignment records");
        }
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
                    store.add_group(txps, &mut records_for_read);
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
        store.add_group(txps, &mut records_for_read);
        records_for_read.clear();
    }

    info!(
        "the alignment file contained {} unmapped read records.",
        num_unmapped.to_formatted_string(&Locale::en)
    );

    Ok(())
}
