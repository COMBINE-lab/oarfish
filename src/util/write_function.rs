use crate::prog_opts::ReadAssignmentProbOut;
use crate::util::oarfish_types::EMInfo;
use crate::util::parquet_utils;
use itertools::izip;

use crate::util::name_store::NameStore;
use arrow_array::ArrayRef;
use either::Either;
use lz4::EncoderBuilder;
use path_tools::WithAdditionalExtension;

use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::{
    fs,
    fs::File,
    fs::OpenOptions,
    fs::create_dir_all,
    io::{self, BufWriter, Write},
};

/// When this run started, as recorded by [`record_start_time`].
static START_TIME: OnceLock<String> = OnceLock::new();

/// The local time in the format salmon uses in `meta_info.json` (that of C's `asctime`, e.g.
/// "Tue Jan 15 11:10:00 2019").
fn meta_info_time_now() -> String {
    chrono::Local::now()
        .format("%a %b %e %H:%M:%S %Y")
        .to_string()
}

/// Records the start time written to `meta_info.json`. Called once, when oarfish starts.
pub fn record_start_time() {
    START_TIME.get_or_init(meta_info_time_now);
}

/// Adds the fields salmon records in `meta_info.json` about the run itself: the version that
/// produced the output, and when the run started and ended (the end being when the metadata is
/// written).
fn add_run_info(info: &mut serde_json::Value) {
    if let Some(m) = info.as_object_mut() {
        let end_time = meta_info_time_now();
        let start_time = START_TIME
            .get()
            .cloned()
            .unwrap_or_else(|| end_time.clone());
        m.insert(
            "oarfish_version".to_string(),
            env!("CARGO_PKG_VERSION").into(),
        );
        m.insert("start_time".to_string(), start_time.into());
        m.insert("end_time".to_string(), end_time.into());
    }
}

pub fn write_single_cell_output(
    output: &PathBuf,
    mut info: serde_json::Value,
    header: &noodles_sam::header::Header,
    counts: &sprs::TriMatI<f32, u32>,
) -> io::Result<()> {
    // if there is a parent directory
    if let Some(p) = output.parent() {
        // unless this was a relative path with one component,
        // which we should treat as the file prefix, then grab
        // the non-empty parent and create it.
        if p != Path::new("") {
            create_dir_all(p)?;
        }
    }

    {
        let info_path = output.with_additional_extension(".meta_info.json");
        let write = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(info_path)
            .expect("Couldn't create output file");

        add_run_info(&mut info);
        serde_json::ser::to_writer_pretty(write, &info)?;
    }

    let out_path = output.with_additional_extension(".count.mtx");
    sprs::io::write_matrix_market(out_path, counts)?;

    let out_path = output.with_additional_extension(".features.txt");
    File::create(&out_path)?;
    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    for (rseq, _rmap) in header.reference_sequences().iter() {
        writeln!(writer, "{}", rseq).expect("Couldn't write to output file.");
    }
    Ok(())
}

//this part is taken from dev branch
pub fn write_output(
    output: &PathBuf,
    mut info: serde_json::Value,
    header: &noodles_sam::header::Header,
    counts: &[f64],
    aux_counts: &[crate::util::aux_counts::CountInfo],
) -> io::Result<()> {
    // if there is a parent directory
    if let Some(p) = output.parent() {
        // unless this was a relative path with one component,
        // which we should treat as the file prefix, then grab
        // the non-empty parent and create it.
        if p != Path::new("") {
            create_dir_all(p)?;
        }
    }

    {
        let info_path = output.with_additional_extension(".meta_info.json");
        let write = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(info_path)
            .expect("Couldn't create output file");

        add_run_info(&mut info);
        serde_json::ser::to_writer_pretty(write, &info)?;
    }

    let out_path = output.with_additional_extension(".quant");
    File::create(&out_path)?;

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "tname\tlen\tnum_reads").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.

    for (i, (rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        writeln!(writer, "{}\t{}\t{}", rseq, rmap.length(), counts[i])
            .expect("Couldn't write to output file.");
    }

    // write the auxiliary count info
    let out_path = output.with_additional_extension(".ambig_info.tsv");
    File::create(&out_path)?;

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "unique_reads\tambig_reads\ttotal_reads")
        .expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.

    for (i, (_rseq, _rmap)) in header.reference_sequences().iter().enumerate() {
        let total = aux_counts[i].total_count;
        let unique = aux_counts[i].unique_count;
        let ambig = total.saturating_sub(unique);
        writeln!(writer, "{}\t{}\t{}", unique, ambig, total)
            .expect("Couldn't write to output file.");
    }

    Ok(())
}

#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn write_out_cdf(
    output: &String,
    prob: &str,
    rate: &str,
    bins: &u32,
    alpha: f64,
    beta: f64,
    emi: &EMInfo,
    txps_name: &[String],
) -> io::Result<()> {
    let output_directory = format!("{}/{}/CDFOutput", output, bins);
    fs::create_dir_all(output_directory.clone())?;

    let out_path: String = if prob == "entropy" {
        format!(
            "{}/{}_{}_{}_{}_cdf.tsv",
            output_directory, prob, rate, alpha, beta
        )
    } else {
        format!("{}/{}_{}_cdf.tsv", output_directory, prob, rate)
    };

    File::create(out_path.clone())?;

    let write_cdf = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer_cdf = BufWriter::new(write_cdf);

    writeln!(writer_cdf, "Txps_Name\tCDF_Values").expect("Couldn't write to output file.");
    for (i, txp) in txps_name.iter().enumerate() {
        let cdf_values: String = emi.txp_info[i]
            .coverage_prob
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<String>>()
            .join("\t");

        writeln!(writer_cdf, "{}\t{}", *txp, cdf_values,).expect("Couldn't write to output file.");
    }

    Ok(())
}

/// Writes the bootstrap replicates `columns` (name, values) to `<output_path>.infreps.pq`.
pub(crate) fn write_infrep_file(
    output_path: &Path,
    columns: Vec<(String, ArrayRef)>,
) -> anyhow::Result<()> {
    let output_path = output_path
        .to_path_buf()
        .with_additional_extension(".infreps.pq");
    parquet_utils::write_columns(&output_path, columns)
}

/// Number of decimal places to use when printing assignment probabilities to
/// the `.prob` file. Chosen so that an assignment that passes `--display-thresh`
/// never rounds to `0.000`: with the fixed 3-decimal format an entry could pass
/// the (default `1e-6`) threshold yet print as zero (COMBINE-lab/oarfish#71).
/// Precision tracks the threshold — `decimals = ceil(-log10(display_thresh))` —
/// clamped to `[3, 9]` so it stays at the historical 3 decimals for thresholds
/// >= 1e-3 and never grows without bound (e.g. the `none` sentinel).
fn prob_display_decimals(display_thresh: f64) -> usize {
    if display_thresh > 0.0 && display_thresh.is_finite() {
        (-display_thresh.log10()).ceil().clamp(3.0, 9.0) as usize
    } else {
        9
    }
}

/// The probability that a read comes from each transcript it aligns to, from the weight
/// (abundance times alignment probability) of each of its alignments.
///
/// A read can have several alignments to one transcript; their weights are added, so that
/// each transcript is listed once, with the probability of the transcript and not of one
/// alignment to it. Transcripts are in the order of their first alignment. Those below
/// `display_thresh` are left out, and the rest are renormalized to sum to 1.
fn read_assignment_probs(
    weights: impl Iterator<Item = (usize, f64)>,
    display_thresh: f64,
    txps: &mut Vec<usize>,
    txp_probs: &mut Vec<f64>,
) {
    txps.clear();
    txp_probs.clear();
    for (target_id, w) in weights {
        // a read has few alignments, so a scan is cheaper than a map
        match txps.iter().position(|&t| t == target_id) {
            Some(i) => txp_probs[i] += w,
            None => {
                txps.push(target_id);
                txp_probs.push(w);
            }
        }
    }
    let denom: f64 = txp_probs.iter().sum();

    let mut kept = 0;
    let mut kept_sum = 0.0_f64;
    for i in 0..txps.len() {
        let nprob = (txp_probs[i] / denom).clamp(0.0, 1.0);
        if nprob >= display_thresh {
            txps[kept] = txps[i];
            txp_probs[kept] = nprob;
            kept_sum += nprob;
            kept += 1;
        }
    }
    txps.truncate(kept);
    txp_probs.truncate(kept);
    for p in txp_probs.iter_mut() {
        *p /= kept_sum;
    }
}

pub fn write_out_prob(
    output: &PathBuf,
    emi: &EMInfo,
    counts: &[f64],
    names_vec: NameStore,
    txps_name: &[String],
    display_thresh: f64,
) -> anyhow::Result<()> {
    if let Some(p) = output.parent() {
        // unless this was a relative path with one component,
        // which we should treat as the file prefix, then grab
        // the non-empty parent and create it.
        if p != Path::new("") {
            create_dir_all(p)?;
        }
    }

    let compressed = matches!(
        emi.eq_map.filter_opts.write_assignment_probs_type,
        Some(ReadAssignmentProbOut::Compressed)
    );

    let extension = if compressed { ".prob.lz4" } else { ".prob" };
    let out_path = output.with_additional_extension(extension);
    File::create(&out_path)?;

    let write_prob = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");

    let mut writer_prob = if compressed {
        Either::Right(EncoderBuilder::new().level(4).build(write_prob)?)
    } else {
        Either::Left(BufWriter::with_capacity(1024 * 1024, write_prob))
    };

    writeln!(writer_prob, "{}\t{}", txps_name.len(), emi.eq_map.len())
        .expect("couldn't write to prob output file");
    for tname in txps_name {
        writeln!(writer_prob, "{}", tname).expect("couldn't write to prob output file");
    }

    let model_coverage = emi.eq_map.filter_opts.model_coverage;
    //let names_vec = emi.eq_map.take_read_names_vec()?;

    let names_iter = names_vec.into_iter();

    let mut txps = Vec::<usize>::new();
    let mut txp_probs = Vec::<f64>::new();

    // print enough decimals that a probability passing `display_thresh` is not
    // shown as 0.000 (see `prob_display_decimals` / oarfish#71).
    let prob_decimals = prob_display_decimals(display_thresh);

    for ((alns, probs, coverage_probs), name) in izip!(emi.eq_map.iter_input_order(), names_iter) {
        let rn = name.expect("could not extract read name from file");
        let read = rn.trim_end_matches('\0');

        write!(writer_prob, "{}\t", read).expect("couldn't write to prob output file");

        let weights = izip!(alns, probs, coverage_probs).map(|(a, p, cp)| {
            let target_id = a.ref_id as usize;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            (target_id, counts[target_id] * (*p as f64) * cov_prob)
        });
        read_assignment_probs(weights, display_thresh, &mut txps, &mut txp_probs);

        let txp_ids = txps
            .iter()
            .map(|x| format!("{}", x))
            .collect::<Vec<String>>()
            .join("\t");
        let prob_vals = txp_probs
            .iter()
            .map(|x| format!("{:.*}", prob_decimals, x))
            .collect::<Vec<String>>()
            .join("\t");
        writeln!(writer_prob, "{}\t{}\t{}", txps.len(), txp_ids, prob_vals)
            .expect("couldn't write to prob output file");
    }

    if let Either::Right(lz4) = writer_prob {
        let (_output, result) = lz4.finish();
        result?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{add_run_info, prob_display_decimals, read_assignment_probs};

    fn assign(weights: &[(usize, f64)], thresh: f64) -> (Vec<usize>, Vec<f64>) {
        let (mut t, mut p) = (vec![99], vec![9.0]);
        read_assignment_probs(weights.iter().copied(), thresh, &mut t, &mut p);
        (t, p)
    }

    #[test]
    fn alignments_to_one_transcript_are_summed() {
        // two alignments to transcript 0: one entry, not two of 0.5
        assert_eq!(assign(&[(0, 1.0), (0, 1.0)], 1e-6), (vec![0], vec![1.0]));
        // transcripts stay in the order of their first alignment
        let (t, p) = assign(&[(3, 1.0), (1, 2.0), (3, 1.0)], 1e-6);
        assert_eq!(t, [3, 1]);
        assert_eq!(p, [0.5, 0.5]);
    }

    #[test]
    fn the_threshold_applies_to_the_transcript_total() {
        // each alignment to transcript 5 is below 0.1, their sum is not
        let (t, p) = assign(&[(2, 0.85), (5, 0.075), (5, 0.075)], 0.1);
        assert_eq!(t, [2, 5]);
        assert!((p[0] - 0.85).abs() < 1e-12 && (p[1] - 0.15).abs() < 1e-12);
        // a transcript below the threshold is left out and the rest sum to 1
        let (t, p) = assign(&[(2, 0.95), (5, 0.05)], 0.1);
        assert_eq!((t, p), (vec![2], vec![1.0]));
    }

    #[test]
    fn run_info_matches_salmon_fields() {
        let mut info = serde_json::json!({ "prob_model": "no_coverage" });
        add_run_info(&mut info);
        assert_eq!(info["prob_model"], "no_coverage");
        assert_eq!(info["oarfish_version"], env!("CARGO_PKG_VERSION"));
        for k in ["start_time", "end_time"] {
            // C's asctime layout, e.g. "Tue Jan 15 11:10:00 2019" or "Fri Nov  6 09:05:00 2026".
            let t = info[k].as_str().expect("a string");
            let f: Vec<&str> = t.split_whitespace().collect();
            assert_eq!(t.len(), 24, "{t}");
            assert_eq!(f.len(), 5, "{t}");
            assert_eq!(f[3].len(), 8, "{t}");
        }
    }

    #[test]
    fn prob_decimals_track_threshold() {
        // >= 1e-3 keeps the historical 3-decimal output
        assert_eq!(prob_display_decimals(1e-2), 3);
        assert_eq!(prob_display_decimals(1e-3), 3);
        assert_eq!(prob_display_decimals(0.5), 3);
        // below 1e-3 the precision grows to match the threshold
        assert_eq!(prob_display_decimals(1e-6), 6); // the default
        assert_eq!(prob_display_decimals(1e-4), 4);
        // capped at 9 (e.g. the `none` sentinel = f64::MIN_POSITIVE)
        assert_eq!(prob_display_decimals(1e-12), 9);
        assert_eq!(prob_display_decimals(f64::MIN_POSITIVE), 9);
        // degenerate inputs fall back to the max precision
        assert_eq!(prob_display_decimals(0.0), 9);
    }

    #[test]
    fn passing_prob_never_prints_as_zero() {
        // an entry at the (default) threshold must not render as all-zeros
        let thresh = 1e-6;
        let d = prob_display_decimals(thresh);
        let s = format!("{:.*}", d, thresh);
        assert!(
            s.chars().any(|c| c.is_ascii_digit() && c != '0'),
            "value at threshold printed as zero: {s}"
        );
    }
}
