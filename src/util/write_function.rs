use crate::util::oarfish_types::EMInfo;
use crate::util::parquet_utils;

use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
};
use path_tools::WithAdditionalExtension;

use std::path::{Path, PathBuf};
use std::{
    fs,
    fs::create_dir_all,
    fs::File,
    fs::OpenOptions,
    io::{self, BufWriter, Write},
};

//this part is taken from dev branch
pub fn write_output(
    output: &PathBuf,
    info: serde_json::Value,
    header: &noodles_sam::header::Header,
    counts: &[f64],
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
