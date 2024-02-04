use crate::util::oarfish_types::{TranscriptInfo, EMInfo};
use std::{
    fs,
    fs::OpenOptions,
    fs::File,
    io::{self, BufWriter, Write},
};
use noodles_sam as sam;

//this part is taken from dev branch
pub fn write_out_count(
    output: &String,
    prob_flag: &bool, 
    bins: &u32, 
    header: &noodles_sam::header::Header,
    counts: &[f64],
) -> io::Result<()> {

    let prob = if *prob_flag {"binomial"} else {"NoCoverage"};

    let out_path: String;
    let output_directory = format!("{}/{}/QuantOutput", output, bins);
    fs::create_dir_all(output_directory.clone())?;
    out_path = format!("{}/{}_quant.tsv", output_directory, prob);
    File::create(out_path.clone())?;

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
        writeln!(
            writer,
            "{}\t{}\t{}",
            rseq,
            rmap.length(),
            counts[i]
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}


pub fn write_out_cdf(
    output: &String,
    prob: &str, 
    rate: &str, 
    bins: &u32, 
    alpha: f64, 
    beta: f64,
    emi: &EMInfo,
    txps_name: &Vec<String>,
) -> io::Result<()> {

    let out_path: String;
    let output_directory = format!("{}/{}/CDFOutput", output, bins);
    fs::create_dir_all(output_directory.clone())?;
    if prob == "entropy" {
        out_path = format!("{}/{}_{}_{}_{}_cdf.tsv", output_directory, prob, rate, alpha, beta);
    } else {
        out_path = format!("{}/{}_{}_cdf.tsv", output_directory, prob, rate);
    }
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
        let cdf_values: String = emi.txp_info[i].coverage_prob.iter()
            .map(|value| value.to_string())
            .collect::<Vec<String>>()
            .join("\t");

        writeln!(
            writer_cdf,
            "{}\t{}",
            *txp,
            cdf_values,
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}
