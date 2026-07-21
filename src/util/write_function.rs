use crate::prog_opts::ReadAssignmentProbOut;
use crate::util::oarfish_types::EMInfo;
use crate::util::parquet_utils;
use itertools::izip;

use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
};
use either::Either;
use lz4::EncoderBuilder;
use path_tools::WithAdditionalExtension;
use swapvec::SwapVec;

use std::path::{Path, PathBuf};
use std::{
    fs,
    fs::File,
    fs::OpenOptions,
    fs::create_dir_all,
    io::{self, BufWriter, Write},
};

pub fn write_single_cell_output(
    output: &PathBuf,
    info: serde_json::Value,
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
    info: serde_json::Value,
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

pub(crate) fn write_infrep_file(
    output_path: &Path,
    fields: Vec<Field>,
    chunk: Chunk<Box<dyn Array>>,
) -> anyhow::Result<()> {
    let output_path = output_path
        .to_path_buf()
        .with_additional_extension(".infreps.pq");
    let schema = Schema::from(fields);
    parquet_utils::write_chunk_to_file(output_path.to_str().unwrap(), schema, chunk)
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

pub fn write_out_prob(
    output: &PathBuf,
    emi: &EMInfo,
    counts: &[f64],
    names_vec: SwapVec<String>,
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

    for ((alns, probs, coverage_probs), name) in izip!(emi.eq_map.iter(), names_iter) {
        let mut denom = 0.0_f64;

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            denom += counts[target_id] * prob * cov_prob;
        }

        let rn = name.expect("could not extract read name from file");
        let read = rn.trim_end_matches('\0');

        write!(writer_prob, "{}\t", read).expect("couldn't write to prob output file");

        txps.clear();
        txp_probs.clear();

        let mut denom2 = 0.0_f64;

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            let nprob = ((counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);

            if nprob >= display_thresh {
                txps.push(target_id);
                txp_probs.push(nprob);
                denom2 += nprob;
            }
        }

        for p in txp_probs.iter_mut() {
            *p /= denom2;
        }

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

/// Write raw candidate-relative signals for coverage-model development.
///
/// The candidate arrays preserve equivalence classes while deterministic
/// sampling bounds output volume. Probabilities are the alignment-score and
/// final coverage terms before abundance is applied by the EM.
pub fn write_coverage_signals(
    output: &PathBuf,
    emi: &EMInfo,
    names_vec: SwapVec<String>,
    txps_name: &[String],
    sample_rate: usize,
) -> anyhow::Result<()> {
    if let Some(parent) = output.parent() {
        if parent != Path::new("") {
            create_dir_all(parent)?;
        }
    }
    let path = output.with_additional_extension(".coverage_signals.tsv");
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, file);
    writeln!(
        writer,
        "read\tread_index\teq_size\ttranscripts\tlengths\tstarts\tends\tstrands\tscore_probs\tcoverage_probs"
    )?;
    let rate = sample_rate.max(1);
    for (read_index, ((alns, score, coverage), name)) in
        emi.eq_map.iter().zip(names_vec.into_iter()).enumerate()
    {
        if read_index % rate != 0 {
            continue;
        }
        let read = name.expect("could not extract read name");
        let join = |values: Vec<String>| values.join(",");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            read.trim_end_matches('\0'),
            read_index,
            alns.len(),
            join(
                alns.iter()
                    .map(|a| txps_name[a.ref_id as usize].clone())
                    .collect()
            ),
            join(
                alns.iter()
                    .map(|a| emi.txp_info[a.ref_id as usize].len.get().to_string())
                    .collect()
            ),
            join(alns.iter().map(|a| a.start.to_string()).collect()),
            join(alns.iter().map(|a| a.end.to_string()).collect()),
            join(alns.iter().map(|a| format!("{:?}", a.strand)).collect()),
            join(score.iter().map(|p| format!("{p:.8}")).collect()),
            join(coverage.iter().map(|p| format!("{p:.12}")).collect()),
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::prob_display_decimals;

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
