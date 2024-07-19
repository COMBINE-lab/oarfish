use crate::util::oarfish_types::EMInfo;
use path_tools::WithAdditionalExtension;
use serde_json::json;
use std::path::{Path, PathBuf};
use std::{
    fs,
    fs::create_dir_all,
    fs::File,
    fs::OpenOptions,
    io::{self, BufWriter, Write},
};
use itertools::izip;
use std::collections::{HashSet, HashMap};

//this part is taken from dev branch
pub fn write_output(
    output: &PathBuf,
    prob_flag: &bool,
    bins: &u32,
    em_info: &EMInfo,
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

    let prob = if *prob_flag {
        "binomial"
    } else {
        "no_coverage"
    };

    {
        let info = json!({
            "prob_model" : prob,
            "num_bins" : bins,
            "filter_options" : em_info.eq_map.filter_opts,
            "discard_table" : em_info.eq_map.discard_table
        });

        let info_path = output.join(format!("{}.meta_info.json", prob));
        let write = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(info_path)
            .expect("Couldn't create output file");

        serde_json::ser::to_writer_pretty(write, &info)?;
    }

    let out_path = output.join(format!("{}.quant", prob));
    File::create(&out_path)?;

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "tname\tlen\tnum_reads\tentropy\tentropy_max\tentropy_ratio\t1-entropy_ratio\tread_count").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.

    for (i, (rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        writeln!(writer, 
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                    rseq, 
                    rmap.length(), 
                    counts[i], 
                    em_info.txp_info[i].entropy, 
                    em_info.txp_info[i].entropy_max, 
                    em_info.txp_info[i].entropy_ratio, 
                    em_info.txp_info[i].entropy_ratio_1,
                    em_info.txp_info[i].num_read)

            .expect("Couldn't write to output file.");
    }

    Ok(())
}



pub fn write_cdf(
    output: &PathBuf,
    em_info: &EMInfo,
    header: &noodles_sam::header::Header,
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

    let out_path = output.join("binomial.cdf");
    File::create(&out_path)?;

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "Txps_Name\tCDF_Values").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (i, (rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        let cdf_values: String = em_info.txp_info[i].coverage_prob.iter()
            .map(|value| value.to_string())
            .collect::<Vec<String>>()
            .join("\t");

        writeln!(
            writer,
            "{}\t{}",
            rseq,
            cdf_values,
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}


pub fn write_info(
    output: &PathBuf,
    prob_flag: &bool,
    em_info: &EMInfo,
    txps_name: &Vec<String>,
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

    let prob = if *prob_flag {
        "binomial"
    } else {
        "no_coverage"
    };


    let out_path = output.join(format!("{}.info", prob));
    File::create(&out_path)?;

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    let out_path_kde = output.join(format!("{}.kde", prob));
    File::create(&out_path_kde)?;

    let write_kde = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path_kde)
        .expect("Couldn't create output file");
    let mut writer_kde = BufWriter::new(write_kde);


    writeln!(writer, "Read\ttname\ttlen\tread_start\tread_end\tAS\tAS_prob\tcoverage_prob\tkde_prob\tfinal_prob\tentropy\tentropy_max\tentropy_ratio\t1-entropy_ratio\tTIN").expect("Couldn't write to output file.");
    writeln!(writer_kde, "Read\ttname\ttlen\tread_len\tkde_prob").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.

    for (alns, probs, cov_probs, as_vec, read_vec, kde_probs) in em_info.eq_map.iter() {

        let mut prob_vec = vec![];
        for (a, p, cov_p, as_val, rname, kdep) in izip!(alns, probs, cov_probs, as_vec, read_vec, kde_probs) {
            prob_vec.push((*p) as f64 * (if prob == "binomial" { *cov_p } else { 1.0 }) * (*kdep));
        }
        let summation: f64 = prob_vec.iter().sum();
        for (a, p, cov_p, as_val, rname, p_final, kde_p) in izip!(alns, probs, cov_probs, as_vec, read_vec, prob_vec, kde_probs) {
            let target_id = a.ref_id as usize;
            let txp_n = txps_name[target_id].clone();
            let rstring = rname;
            let read_start = a.start as u32;
            let read_end = a.end as u32;
            let as_val = as_val;
            let as_prob = *p;
            let txp_len = em_info.txp_info[target_id].len.get();
            let p_final_value = p_final/summation;

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                rstring,
                txp_n,
                txp_len,
                read_start,
                read_end,
                as_val,
                as_prob,
                *cov_p,
                *kde_p,
                p_final_value,
                em_info.txp_info[target_id].entropy,
                em_info.txp_info[target_id].entropy_max,
                em_info.txp_info[target_id].entropy_ratio,
                em_info.txp_info[target_id].entropy_ratio_1,
                em_info.txp_info[target_id].tin,
            )
            .expect("Couldn't write to output file.");

            writeln!(
                writer_kde,
                "{}\t{}\t{}\t{}\t{}",
                rstring,
                txp_n,
                txp_len,
                read_end - read_start,
                *kde_p,
            )
            .expect("Couldn't write to output file.");
        }
    }

    Ok(())
}


pub fn write_EquivalenceClass(
    output: &PathBuf,
    prob_flag: &bool,
    equivalence_classes: &HashMap<Vec<usize>, (usize, Vec<String>)>,
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

    let prob = if *prob_flag {
        "binomial"
    } else {
        "no_coverage"
    };


    let out_path = output.join(format!("{}.EquivalenceClasses", prob));
    File::create(&out_path)?;

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(out_path)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    //writeln!(writer, "Txps_Name\tCDF_Values").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (key, (count, names)) in equivalence_classes.iter() {
        let key_values: String = {
            let key_length = key.len().to_string();
            let values = key.iter().map(|value| value.to_string()).collect::<Vec<String>>().join("\t");
            format!("{}\t{}", key_length, values)
        };

        let read_names: String = names.iter()
            .map(|value| value.to_string())
            .collect::<Vec<String>>()
            .join("\t");
        

        writeln!(
            writer,
            "{}\t{}\t{}",
            key_values,
            count,
            read_names,
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}
