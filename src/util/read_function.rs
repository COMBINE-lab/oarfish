use crate::util::oarfish_types::ShortReadRecord;
use anyhow::bail;
use csv::ReaderBuilder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use tracing::warn;

/// Read the short read quantification from the file `short_read_path`
pub fn read_short_quant_vec(
    short_read_path: &str,
    txps_name: &[String],
) -> anyhow::Result<Vec<f64>> {
    // try to open the short read file
    let file = File::open(short_read_path)?;

    // read the data as a tab-delimited TSV file.
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(file);

    // deserialize CSV records into a Vec of ShortReadRecords
    let records: HashMap<String, ShortReadRecord> = rdr
        .deserialize()
        .collect::<Result<Vec<ShortReadRecord>, csv::Error>>()
        .unwrap_or_else(|err| {
            eprintln!("Failed to deserialize CSV records: {}", err);
            std::process::exit(1);
        })
        .into_iter()
        .map(|rec| (rec.name.clone(), rec))
        .collect();

    // txps_name are the transcript names in the BAM header. We expect
    // that all elements that we see in the short read file should
    // exist in txps_name. If they don't then this is an error.
    {
        // create a scope here to free the HashSet as soon as we are done checking
        let txps_name_set: HashSet<&str> = txps_name.iter().map(|x| x.as_str()).collect();
        if !records
            .iter()
            .all(|(k, _v)| txps_name_set.contains(k.as_str()))
        {
            bail!(
                "There were transcripts in the short read quantification file that didn't appear in the BAM header; cannot proceed."
            );
        }
    }

    // project the short read records to a vector in the same order of
    // txps_name, filling in records missing in the short read HashMap
    // with default values
    let mut num_missing = 0;
    let ordered_rec: Vec<f64> = txps_name
        .iter()
        .map(|name| {
            records.get(name).map_or_else(
                // we are missing this record, return 0 abundance
                || {
                    num_missing += 1;
                    0_f64
                },
                // return the number of reads from the record
                |rec| rec.num_reads,
            )
        })
        .collect();

    if num_missing > 0 {
        warn!(
            "There were {} transcripts appearing in the BAM header but missing from the short read quatifications; they have been assumed to have 0 abunance.",
            num_missing
        );
    }

    Ok(ordered_rec)
}
