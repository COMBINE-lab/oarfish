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

use std::path::{Path, PathBuf};
use std::{
    fs,
    fs::create_dir_all,
    fs::File,
    fs::OpenOptions,
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

pub fn write_out_prob(output: &PathBuf, emi: &EMInfo, txps_name: &[String]) -> io::Result<()> {
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

    let mut txps = Vec::<usize>::new();
    let mut txp_probs = Vec::<f64>::new();

    for (alns, probs, coverage_probs, name) in emi.eq_map.iter_with_names() {
        let mut denom = 0.0_f64;

        for (_a, p, cp) in izip!(alns, probs, coverage_probs) {
            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            denom += prob * cov_prob;
        }

        let read = if let Some(rn) = name {
            rn
        } else {
            "no_read_name_available"
        };
        write!(writer_prob, "{}\t", read).expect("couldn't write to prob output file");

        txps.clear();
        txp_probs.clear();

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            txps.push(target_id);
            txp_probs.push(((prob * cov_prob) / denom).clamp(0.0, 1.0));
        }

        let txp_ids = txps
            .iter()
            .map(|x| format!("{}", x))
            .collect::<Vec<String>>()
            .join("\t");
        let prob_vals = txp_probs
            .iter()
            .map(|x| format!("{:.3}", x))
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
