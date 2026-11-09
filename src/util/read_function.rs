use crate::util::oarfish_types::ShortReadRecord;
use anyhow::bail;
use csv::ReaderBuilder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use tracing::{info, warn};

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
    let mut ordered_rec: Vec<f64> = txps_name
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

    let floored = floor_zero_abundances(&mut ordered_rec);
    if floored > 0 {
        info!(
            "{} transcripts have no short-read abundance; they start the EM at a negligible abundance, so that reads aligned only to them are kept.",
            floored
        );
    }

    Ok(ordered_rec)
}

/// The EM's updates are multiplicative, so a transcript that starts at 0 stays at 0, and the
/// reads aligned only to such transcripts would be lost. This raises every abundance that is not
/// positive to `SHORT_QUANT_FLOOR` times the mean (or to 1 if none is positive): such a
/// transcript still gets the reads only it explains, and otherwise stays negligible and is pruned
/// at the end of the EM. Returns how many abundances it raised.
pub(crate) fn floor_zero_abundances(v: &mut [f64]) -> usize {
    let total: f64 = v.iter().filter(|x| x.is_finite() && **x > 0.0).sum();
    let mean = total / v.len().max(1) as f64;
    let floor = if mean > 0.0 {
        mean * SHORT_QUANT_FLOOR
    } else {
        1.0
    };
    let mut n = 0;
    for x in v.iter_mut() {
        if x.is_nan() || *x <= 0.0 {
            *x = floor;
            n += 1;
        }
    }
    n
}

/// The starting abundance of a transcript the short reads gave none, relative to the mean.
const SHORT_QUANT_FLOOR: f64 = 1e-8;

#[cfg(test)]
mod tests {
    use super::floor_zero_abundances;

    #[test]
    fn zero_abundances_are_floored() {
        let mut v = vec![0.0, 3.0, 3.0, f64::NAN];
        assert_eq!(floor_zero_abundances(&mut v), 2);
        assert_eq!(&v[1..3], &[3.0, 3.0]);
        assert_eq!(v[0], 1.5 * 1e-8);
        assert_eq!(v[3], 1.5 * 1e-8);

        let mut v = vec![0.0; 3];
        assert_eq!(floor_zero_abundances(&mut v), 3);
        assert_eq!(v, [1.0; 3]);

        let mut v: Vec<f64> = vec![];
        assert_eq!(floor_zero_abundances(&mut v), 0);
    }
}
