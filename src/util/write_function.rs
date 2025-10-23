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

//pub fn write_out_prob(
//    output: &PathBuf,
//    emi: &EMInfo,
//    counts: &[f64],
//    names_vec: &SwapVec<String>,
//    txps_name: &[String],
//) -> anyhow::Result<()> {
//    if let Some(p) = output.parent() {
//        // unless this was a relative path with one component,
//        // which we should treat as the file prefix, then grab
//        // the non-empty parent and create it.
//        if p != Path::new("") {
//            create_dir_all(p)?;
//        }
//    }
//
//    let compressed = matches!(
//        emi.eq_map.filter_opts.write_assignment_probs_type,
//        Some(ReadAssignmentProbOut::Compressed)
//    );
//
//    let extension = if compressed { ".prob.lz4" } else { ".prob" };
//    let out_path = output.with_additional_extension(extension);
//    File::create(&out_path)?;
//
//    let write_prob = OpenOptions::new()
//        .write(true)
//        .create(true)
//        .truncate(true)
//        .open(out_path)
//        .expect("Couldn't create output file");
//
//    let mut writer_prob = if compressed {
//        Either::Right(EncoderBuilder::new().level(4).build(write_prob)?)
//    } else {
//        Either::Left(BufWriter::with_capacity(1024 * 1024, write_prob))
//    };
//
//    writeln!(writer_prob, "{}\t{}", txps_name.len(), emi.eq_map.len())
//        .expect("couldn't write to prob output file");
//    for tname in txps_name {
//        writeln!(writer_prob, "{}", tname).expect("couldn't write to prob output file");
//    }
//
//    let model_coverage = emi.eq_map.filter_opts.model_coverage;
//    //let names_vec = emi.eq_map.take_read_names_vec()?;
//
//    let names_iter = names_vec.into_iter();
//
//    let mut txps = Vec::<usize>::new();
//    let mut txp_probs = Vec::<f64>::new();
//
//    for ((alns, probs, coverage_probs), name) in izip!(emi.eq_map.iter(), names_iter) {
//        let mut denom = 0.0_f64;
//
//        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
//            let target_id = a.ref_id as usize;
//            let prob = *p as f64;
//            let cov_prob = if model_coverage { *cp } else { 1.0 };
//            denom += counts[target_id] * prob * cov_prob;
//        }
//
//        let rn = name.expect("could not extract read name from file");
//        let read = rn.trim_end_matches('\0');
//
//        write!(writer_prob, "{}\t", read).expect("couldn't write to prob output file");
//
//        txps.clear();
//        txp_probs.clear();
//
//        const DISPLAY_THRESH: f64 = 0.001;
//        let mut denom2 = 0.0_f64;
//
//        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
//            let target_id = a.ref_id as usize;
//            let prob = *p as f64;
//            let cov_prob = if model_coverage { *cp } else { 1.0 };
//            let nprob = ((counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
//            if nprob >= DISPLAY_THRESH {
//                txps.push(target_id);
//                txp_probs.push(nprob);
//                denom2 += nprob;
//            }
//        }
//
//        for p in txp_probs.iter_mut() {
//            *p /= denom2;
//        }
//
//        let txp_ids = txps
//            .iter()
//            .map(|x| format!("{}", x))
//            .collect::<Vec<String>>()
//            .join("\t");
//        let prob_vals = txp_probs
//            .iter()
//            .map(|x| format!("{:.3}", x))
//            .collect::<Vec<String>>()
//            .join("\t");
//        writeln!(writer_prob, "{}\t{}\t{}", txps.len(), txp_ids, prob_vals)
//            .expect("couldn't write to prob output file");
//    }
//
//    if let Either::Right(lz4) = writer_prob {
//        let (_output, result) = lz4.finish();
//        result?;
//    }
//
//    Ok(())
//}
//
//
//
//
//
//
//pub fn write_out_assignment(
//    output: &PathBuf,
//    emi: &EMInfo,
//    counts: &[usize],
//    names_vec: &SwapVec<String>,
//    txps_name: &[String],
//) -> anyhow::Result<()> {
//    if let Some(p) = output.parent() {
//        // unless this was a relative path with one component,
//        // which we should treat as the file prefix, then grab
//        // the non-empty parent and create it.
//        if p != Path::new("") {
//            create_dir_all(p)?;
//        }
//    }
//
//    let compressed = matches!(
//        emi.eq_map.filter_opts.write_assignment_probs_type,
//        Some(ReadAssignmentProbOut::Compressed)
//    );
//
//    let extension = if compressed { ".hard_assignment.lz4" } else { ".hard_assignment" };
//    let out_path = output.with_additional_extension(extension);
//    File::create(&out_path)?;
//
//    let write_prob = OpenOptions::new()
//        .write(true)
//        .create(true)
//        .truncate(true)
//        .open(out_path)
//        .expect("Couldn't create output file");
//
//    let mut writer_prob = if compressed {
//        Either::Right(EncoderBuilder::new().level(4).build(write_prob)?)
//    } else {
//        Either::Left(BufWriter::with_capacity(1024 * 1024, write_prob))
//    };
//
//
//    let names_iter = names_vec.into_iter();
//    let count_iter = counts.into_iter();
//
//    for ((alns, probs, coverage_probs), name, tx_id) in izip!(emi.eq_map.iter(), names_iter, count_iter) {
//
//        let rn = name.expect("could not extract read name from file");
//        let read = rn.trim_end_matches('\0');
//        let tx_name = &txps_name[*tx_id];
//
//        write!(writer_prob, "{}\t{}\n", read, tx_name).expect("couldn't write to prob output file");
//
//    }
//
//
//    if let Either::Right(lz4) = writer_prob {
//        let (_output, result) = lz4.finish();
//        result?;
//    }
//
//    Ok(())
//}




pub fn write_out_prob_assignment(
    output: &PathBuf,
    emi: &EMInfo,
    counts: &[f64],
    assignment_tid: &[usize],
    read_candidate: &Vec<Vec<usize>>,
    names_vec: SwapVec<String>,
    txps_name: &[String],
    prob_flag: bool,
    assign_flag:bool,
) -> anyhow::Result<()> {

    let names_iter: Vec<String> = names_vec
        .into_iter()
        .map(|opt| opt.unwrap_or_default()) // each element is Option<String>
        .collect();

    if prob_flag{
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


        let mut txps = Vec::<usize>::new();
        let mut txp_probs = Vec::<f64>::new();

        for ((alns, probs, coverage_probs), name) in izip!(emi.eq_map.iter(), &names_iter) {
            let mut denom = 0.0_f64;

            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                denom += counts[target_id] * prob * cov_prob;
            }

            //let rn = name.expect("could not extract read name from file");
            let read = name.trim_end_matches('\0');

            write!(writer_prob, "{}\t", read).expect("couldn't write to prob output file");

            txps.clear();
            txp_probs.clear();

            const DISPLAY_THRESH: f64 = 0.001;
            let mut denom2 = 0.0_f64;

            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                let nprob = ((counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
                if nprob >= DISPLAY_THRESH {
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
    }

    if assign_flag{
        if let Some(p) = output.parent() {
        // unless this was a relative path with one component,
        // which we should treat as the file prefix, then grab
        // the non-empty parent and create it.
        if p != Path::new("") {
            create_dir_all(p)?;
            }
        }

        let compressed = matches!(
            emi.eq_map.filter_opts.collapsed_gibbs_sampler_type,
            Some(ReadAssignmentProbOut::Compressed)
        );

        let extension = if compressed { ".hard_assignment.lz4" } else { ".hard_assignment" };
        let out_path = output.with_additional_extension(extension);
        File::create(&out_path)?;

        let write_assign = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(out_path)
            .expect("Couldn't create output file");

        let mut write_assign = if compressed {
            Either::Right(EncoderBuilder::new().level(4).build(write_assign)?)
        } else {
            Either::Left(BufWriter::with_capacity(1024 * 1024, write_assign))
        };


        writeln!(write_assign, "{}\t{}", txps_name.len(), emi.eq_map.len())
            .expect("couldn't write to prob output file");
        for tname in txps_name {
            writeln!(write_assign, "{}", tname).expect("couldn't write to prob output file");
        }


        let assignment_tid_iter = assignment_tid.into_iter();

        for (name, tx_id) in izip!(&names_iter, assignment_tid_iter) {

            //let rn = name.expect("could not extract read name from file");
            let read = name.trim_end_matches('\0');

            write!(write_assign, "{}\t{}\n", read, tx_id).expect("couldn't write to prob output file");

        }

        if let Either::Right(lz4) = write_assign {
            let (_output, result) = lz4.finish();
            result?;
        }



        //write the all gibbs iterations
        let extension_gibbs = if compressed { ".Gibbs_iterations.lz4" } else { ".Gibbs_iterations" };
        let out_path_gibbs = output.with_additional_extension(extension_gibbs);
        File::create(&out_path_gibbs)?;

        let write_gibbs = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(out_path_gibbs)
            .expect("Couldn't create output file");

        let mut write_gibbs = if compressed {
            Either::Right(EncoderBuilder::new().level(4).build(write_gibbs)?)
        } else {
            Either::Left(BufWriter::with_capacity(1024 * 1024, write_gibbs))
        };


        writeln!(write_gibbs, "{}\t{}", txps_name.len(), emi.eq_map.len())
            .expect("couldn't write to prob output file");
        for tname in txps_name {
            writeln!(write_gibbs, "{}", tname).expect("couldn't write to prob output file");
        }


        for gibbs_iteration in read_candidate.iter(){
            let gibbs_join = gibbs_iteration
            .iter()
            .map(|x| format!("{}", x))
            .collect::<Vec<String>>()
            .join("\t");
            write!(write_gibbs, "{}\n", gibbs_join).expect("couldn't write to prob output file");
        }

        if let Either::Right(lz4) = write_gibbs {
            let (_output, result) = lz4.finish();
            result?;
        }

    }

    Ok(())
}






pub fn write_out_prob_assignment_greedy(
    output: &PathBuf,
    emi: &EMInfo,
    counts: &[f64],
    assignment_tid: &[usize],
    names_vec: SwapVec<String>,
    txps_name: &[String],
    prob_flag: bool,
    assign_flag:bool,
) -> anyhow::Result<()> {

    let names_iter: Vec<String> = names_vec
        .into_iter()
        .map(|opt| opt.unwrap_or_default()) // each element is Option<String>
        .collect();

    if prob_flag{
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


        let mut txps = Vec::<usize>::new();
        let mut txp_probs = Vec::<f64>::new();

        for ((alns, probs, coverage_probs), name) in izip!(emi.eq_map.iter(), &names_iter) {
            let mut denom = 0.0_f64;

            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                denom += counts[target_id] * prob * cov_prob;
            }

            //let rn = name.expect("could not extract read name from file");
            let read = name.trim_end_matches('\0');

            write!(writer_prob, "{}\t", read).expect("couldn't write to prob output file");

            txps.clear();
            txp_probs.clear();

            const DISPLAY_THRESH: f64 = 0.001;
            let mut denom2 = 0.0_f64;

            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                let nprob = ((counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
                if nprob >= DISPLAY_THRESH {
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
    }

    if assign_flag{
        if let Some(p) = output.parent() {
        // unless this was a relative path with one component,
        // which we should treat as the file prefix, then grab
        // the non-empty parent and create it.
        if p != Path::new("") {
            create_dir_all(p)?;
            }
        }

        let compressed = matches!(
            emi.eq_map.filter_opts.collapsed_gibbs_sampler_type,
            Some(ReadAssignmentProbOut::Compressed)
        );

        let extension = if compressed { ".hard_assignment.lz4" } else { ".hard_assignment" };
        let out_path = output.with_additional_extension(extension);
        File::create(&out_path)?;

        let write_assign = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(out_path)
            .expect("Couldn't create output file");

        let mut write_assign = if compressed {
            Either::Right(EncoderBuilder::new().level(4).build(write_assign)?)
        } else {
            Either::Left(BufWriter::with_capacity(1024 * 1024, write_assign))
        };


        writeln!(write_assign, "{}\t{}", txps_name.len(), emi.eq_map.len())
            .expect("couldn't write to prob output file");
        for tname in txps_name {
            writeln!(write_assign, "{}", tname).expect("couldn't write to prob output file");
        }


        let assignment_tid_iter = assignment_tid.into_iter();

        for (name, tx_id) in izip!(&names_iter, assignment_tid_iter) {

            //let rn = name.expect("could not extract read name from file");
            let read = name.trim_end_matches('\0');

            write!(write_assign, "{}\t{}\n", read, tx_id).expect("couldn't write to prob output file");

        }

        if let Either::Right(lz4) = write_assign {
            let (_output, result) = lz4.finish();
            result?;
        }

    }

    Ok(())
}