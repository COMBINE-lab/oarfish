extern crate libc;

mod variables;
mod common_functions;
mod nocoverage_prob;
mod uniform_prob;
mod poisson_prob;
mod lander_waterman_prob;
mod multinomial_prob;
mod binomial_continuous;
mod entropy_distribution;
//mod kde_prob;

use clap::Parser;
use statrs::function::beta;

use core::panic;
use std::{
    fs::File,
    io::{self, BufReader},
    vec,
};

use noodles_bam as bam;

use std::time::Instant;
use rand::Rng;




/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser, required = true)]
    alignments: String,
    #[clap(short, long, value_parser, required = true)]
    output_path: String,
    #[clap(short, long, value_parser, default_value_t = String::from("yes"))]
    coverage: String,
    #[clap(short, long, value_parser, default_value_t = String::from("multinomial"))]
    prob: String,
    #[clap(short, long, value_parser, default_value_t = String::from("dr"))]
    rate: String,
    #[clap(short, long, value_parser, default_value_t = 1)]
    bins: u32,
    #[clap(short, long, value_parser, default_value_t = 1)]
    threads: usize,
    #[clap(short, long, value_parser, default_value_t = String::from("none"))]
    short_quant: String,
    #[clap(short, long, value_parser, default_value_t = 0.3)]
    alpha: f64,
    #[clap(short, long, value_parser, default_value_t = 20.0)]
    beta: f64,

    // Maximum allowable distance of the right-most end of an alignment from the 3' transcript end
    #[clap(short, long, value_parser, default_value_t = u32::MAX as i64)]
    three_prime_clip: i64,
    // Maximum allowable distance of the left-most end of an alignment from the 5' transcript end
    #[clap(short, long, value_parser, default_value_t = u32::MAX)]
    five_prime_clip: u32,
    // Fraction of the best possible alignment score that a secondary alignment must have for
    // consideration
    #[clap(short, long, value_parser, default_value_t = 0.95)]
    score_threshold: f32,
    // Fraction of a query that must be mapped within an alignemnt to consider the alignemnt
    // valid
    #[clap(short, long, value_parser, default_value_t = 0.5)]
    min_aligned_fraction: f32,
    // Minimum number of nucleotides in the aligned portion of a read
    #[clap(short = 'l', long, value_parser, default_value_t = 50)]
    min_aligned_len: u32,
    // Allow both forward-strand and reverse-complement alignments
    #[clap(short = 'n', long, value_parser, default_value_t = true)]
    allow_negative_strand: bool,
}



#[inline]
fn m_step(
    eq_map: &variables::InMemoryAlignmentStore,
    tinfo: &[variables::TranscriptInfo],
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) -> (usize, f64) {
    let mut dropped: usize=0;
    let mut likelihood: f64 =0.0;
    let count_summation: f64 = prev_count.iter().sum();

    for (alns, probs, coverage_probs) in eq_map.iter() {
        let mut denom = 0.0_f64;
        let mut prob_summation = 0.0_f64;
        for (a, (p, cp)) in alns.iter().zip(probs.iter().zip(coverage_probs.iter())) {
            let target_id: usize = a.reference_sequence_id().unwrap();
            let prob = *p as f64;
            //let cov_prob = if tinfo[target_id].coverage < 0.5 { 1e-5 } else { 1.0 };//powf(2.0) as f64;
            //let cov_prob = tinfo[target_id].coverage;//powf(2.0) as f64;
            let cov_prob = *cp;
            //let cov_prob = 1.0;
            denom += prev_count[target_id] * prob * cov_prob;
            prob_summation += prob * cov_prob;
        }


        if denom > 1e-30_f64 {
            likelihood += (denom / (prob_summation * count_summation)).ln();
            for (a, (p,cp)) in alns.iter().zip(probs.iter().zip(coverage_probs.iter())) {
                let target_id = a.reference_sequence_id().unwrap();
                let prob = *p as f64;
                //let cov_prob = if tinfo[target_id].coverage < 0.5 { 1e-5 } else { 1.0 };//powf(2.0) as f64;
                //let cov_prob = tinfo[target_id].coverage;
                let cov_prob = *cp;
                //let cov_prob = 1.0;
                curr_counts[target_id] += (prev_count[target_id] * prob * cov_prob) / denom;
            }
        } else {
            dropped +=1;
        }
    }
    (dropped, likelihood)
    //eprintln!("dropped:{:?}", dropped);
}

//fn em(em_info: &variables::EMInfo, short_read_path: String, txps_name: &Vec<String>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>) {
//    let eq_map = em_info.eq_map;
//    let tinfo = em_info.txp_info;
//    let max_iter = em_info.max_iter;
//    let total_weight: f64 = eq_map.num_aligned_reads() as f64;
//    let mut num_dropped_reads: usize = 0;
//
//    let mut current_vec: Vec<Vec<f64>> = Vec::new();
//    let mut initial_vec: Vec<Vec<f64>> = Vec::new();
//    let mut likelihood_vec: Vec<f64> = Vec::new();
//
//    for num_initial in 0..30 {
//        // init
//        let mut prev_counts: Vec<f64> = Vec::new();
//        if num_initial == 0 {
//            prev_counts = common_functions::short_quant_vec(short_read_path.clone(), txps_name);
//
//        } else if num_initial == 1 {
//            prev_counts = common_functions::short_quant_vec(short_read_path.clone(), txps_name);
//            // Replace zero values with 1e-8 using map
//            let modified_values: Vec<f64> = prev_counts
//            .iter_mut()
//            .map(| value| if *value == 0.0 { 1e-8 } else { *value })
//            .collect();
//            // Update the original vector with modified values
//            prev_counts.copy_from_slice(&modified_values);
//
//        } else if num_initial == 2 {
//            prev_counts = common_functions::short_quant_vec(short_read_path.clone(), txps_name);
//            // Replace zero values with 1e-8 using map
//            let modified_values: Vec<f64> = prev_counts
//            .iter_mut()
//            .map(| value| if *value == 0.0 { 1e-20 } else { *value })
//            .collect();
//            // Update the original vector with modified values
//            prev_counts.copy_from_slice(&modified_values);
//
//        } else if num_initial == 3 {
//            prev_counts = common_functions::short_quant_vec(short_read_path.clone(), txps_name);
//            // Replace zero values with 1e-8 using map
//            let modified_values: Vec<f64> = prev_counts
//            .iter_mut()
//            .map(| value| if *value == 0.0 { 1e-30 } else { *value })
//            .collect();
//            // Update the original vector with modified values
//            prev_counts.copy_from_slice(&modified_values);
//
//        } else if num_initial == 4 {
//            let avg = total_weight / (tinfo.len() as f64);
//            prev_counts = vec![avg; tinfo.len()];
//
//        } else {
//            let mut rng = rand::thread_rng();
//            for _ in 0..tinfo.len() {
//                let mut value: f64 = rng.gen_range(0.1..1.0);
//                while value.abs() < 1e-8 {
//                    value = rng.gen_range(0.1..1.0);
//                }
//                prev_counts.push(value);
//            }
//            let sum: f64 = prev_counts.iter().sum();
//            let scale_factor = total_weight / sum;
//            prev_counts.iter_mut().for_each(|value| *value *= scale_factor);
//        }
//
//        initial_vec.push(prev_counts.clone());
//
//        let mut curr_counts = vec![0.0f64; tinfo.len()];
//
//        let mut rel_diff = 0.0_f64;
//        let mut niter = 0_u32;
//
//        while niter < max_iter {
//            m_step(eq_map, tinfo, &mut prev_counts, &mut curr_counts);
//
//            //std::mem::swap(&)
//            for i in 0..curr_counts.len() {
//                if prev_counts[i] > 1e-8 {
//                    let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
//                    rel_diff = if rel_diff > rd { rel_diff } else { rd };
//                }
//            }
//
//            std::mem::swap(&mut prev_counts, &mut curr_counts);
//            curr_counts.fill(0.0_f64);
//
//            if (rel_diff < 1e-3) && (niter > 50) {
//                break;
//            }
//            niter += 1;
//            if niter % 10 == 0 {
//                eprintln!("iteration {}; rel diff {}", niter, rel_diff);
//            }
//            rel_diff = 0.0_f64;
//        }
//
//        prev_counts.iter_mut().for_each(|x| {
//            if *x < 1e-8 {
//                *x = 0.0
//            }
//        });
//        let mut likelihood = 0.0_f64;
//        (num_dropped_reads, likelihood) = m_step(eq_map, tinfo, &mut prev_counts, &mut curr_counts);
//        current_vec.push(curr_counts);
//        likelihood_vec.push(likelihood);
//    }
//
//    (current_vec, initial_vec, likelihood_vec)
//}


fn em(em_info: &variables::EMInfo, short_read_path: String, txps_name: &Vec<String>) -> (Vec<f64>, usize, f64) {
    let eq_map = em_info.eq_map;
    let tinfo = em_info.txp_info;
    let mut max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;
    let mut num_dropped_reads: usize = 0;

    // init
    let mut prev_counts: Vec<f64>;
    if short_read_path == "none" {
        let avg = total_weight / (tinfo.len() as f64);
        prev_counts = vec![avg; tinfo.len()];
    } else {
        prev_counts = common_functions::short_quant_vec(short_read_path, txps_name);

        //// Replace zero values with 1e-8 using map
        //let modified_values: Vec<f64> = prev_counts
        //.iter_mut()
        //.map(| value| if *value == 0.0 { 1e-8 } else { *value })
        //.collect();
//
        //// Update the original vector with modified values
        //prev_counts.copy_from_slice(&modified_values);

        //max_iter = 1;
    }

    let mut curr_counts = vec![0.0f64; tinfo.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(eq_map, tinfo, &mut prev_counts, &mut curr_counts);

        //std::mem::swap(&)
        for i in 0..curr_counts.len() {
            if prev_counts[i] > 1e-8 {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if (rel_diff < 1e-3) && (niter > 50) {
            break;
        }
        niter += 1;
        if niter % 10 == 0 {
            eprintln!("iteration {}; rel diff {}", niter, rel_diff);
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| {
        if *x < 1e-8 {
            *x = 0.0
        }
    });
    let mut likelihood = 0.0_f64;
    (num_dropped_reads, likelihood) = m_step(eq_map, tinfo, &mut prev_counts, &mut curr_counts);

    (curr_counts, num_dropped_reads, likelihood)
}


fn normalize_read_probs(eq_map: &mut variables::InMemoryAlignmentStore, txp_info: &Vec<variables::TranscriptInfo>, coverage: bool, prob: &str, rate: &str) {
    //eprintln!("in the normalized function");
    let mut normalize_probs_temp: Vec<f64> = vec![];
    let mut kde_c_tinfo: Vec<usize> = vec![0 ; txp_info.len()];
    let mut normalized_coverage_prob: Vec<f64> = vec![];

    //eprintln!("before for loop");
    for (alns, _as_probs, _coverage_prob) in eq_map.iter() {
        //eprintln!("before if");
        let mut normalized_prob_section: Vec<f64> = vec![];
        if coverage {
            for a in alns.iter() {
                let target_id = a.reference_sequence_id().unwrap();
                let mut cov_prob: f64;
                if prob == "kde" || (prob == "kde_c" && rate == "1D") {
                    let start_aln = a.alignment_start().unwrap().get() as usize;
                    let end_aln = a.alignment_end().unwrap().get() as usize;
                    let coverage_prob: &Vec<f32> = &txp_info[target_id].coverage_prob;
                    //cov_prob = (tinfo[target_id].coverage_prob[end_aln] - tinfo[target_id].coverage_prob[start_aln]) as f64;
                    cov_prob = (get_coverag_prob(coverage_prob, end_aln) - get_coverag_prob(coverage_prob, start_aln)) as f64;
                
                } else if prob == "multinomial" || prob == "LW" || prob == "DisPoisson" || prob == "ConPoisson" || prob == "binomial" || prob == "entropy" {
                    //eprintln!("it is here");
                    let start_aln = a.alignment_start().unwrap().get() as usize;
                    let end_aln = a.alignment_end().unwrap().get() as usize;
                    //eprintln!("1: {:?}, 2: {:?}", tinfo[target_id].coverage_prob[end_aln], tinfo[target_id].coverage_prob[start_aln]);
                    cov_prob = (txp_info[target_id].coverage_prob[end_aln] - txp_info[target_id].coverage_prob[start_aln]) as f64;
                    //cov_prob = cov_prob * (end_aln - start_aln) as f64 / txp_info[target_id].len.get() as f64;
                    if cov_prob.is_nan() || cov_prob.is_infinite() {
                        eprintln!("cov_prob: {:?}", cov_prob);
                        eprintln!("end: {:?}", txp_info[target_id].coverage_prob[end_aln]);
                        eprintln!("start: {:?}", txp_info[target_id].coverage_prob[end_aln]);
                    }
                    
                } else if prob == "kde_c" && rate == "2D" {
                    let coverage_prob: &Vec<f32> = &txp_info[target_id].coverage_prob;
                    cov_prob = coverage_prob[kde_c_tinfo[target_id]] as f64;
                    kde_c_tinfo[target_id] += 1;
                    
                } else {
                    let start_aln = a.alignment_start().unwrap().get() as usize;
                    let end_aln = a.alignment_end().unwrap().get() as usize;
                    cov_prob = txp_info[target_id].coverage;
                    //cov_prob = cov_prob * (end_aln - start_aln) as f64 / txp_info[target_id].len.get() as f64;
                }
                //} else {
                //    let start_aln = a.alignment_start().unwrap().get() as usize;
                //    let end_aln = a.alignment_end().unwrap().get() as usize;
                //    cov_prob = (end_aln - start_aln) as f64 / txp_info[target_id].len.get() as f64;
                //}
                
                if cov_prob.is_nan() || cov_prob.is_infinite() {

                    panic!("Error: Invalid result. normalize_read_probs function.");
                }
                normalize_probs_temp.push(cov_prob);
            }
            let sum_normalize_probs_temp: f64 = if normalize_probs_temp.iter().sum::<f64>() > 0.0 {normalize_probs_temp.iter().sum()} else {1.0};
            normalized_prob_section = normalize_probs_temp.iter().map(|&prob| prob/sum_normalize_probs_temp).collect();

        } else {
            normalized_prob_section = vec![1.0 as f64; alns.len()];
        }
        normalized_coverage_prob.extend(normalized_prob_section);
        normalize_probs_temp.clear();
    }
    eq_map.coverage_probabilities = normalized_coverage_prob;
    //eprintln!("after for loop");
}


fn get_coverag_prob(coverage_prob: &Vec<f32>, position: usize) -> f32 {

    let mut value = std::f32::NAN;

    for index in (0..=position).rev() {
        if let Some(&element) = coverage_prob.get(index) {
            if !element.is_nan() {
                value = element;
                break;  // Exit the loop once the first valid value is found
            }
        }
    }

    value
}


fn num_discarded_reads_fun(store: &variables::InMemoryAlignmentStore, t: &Vec<variables::TranscriptInfo>, num_bins: &u32) -> usize {
    //eprintln!("in the num_discarded_reads_fun");
    let mut num_discarded_reads: usize = 0;


    for (alns, _probs, _coverage_prob) in store.iter() {
        //eprintln!("in the for loop");
        let mut num_discarded_alignment: usize = 0;

        for a in alns.iter() {
            let target_id: usize = a.reference_sequence_id().unwrap();

            let transcript_len = t[target_id].len.get() as u32; //transcript length
            let bin_size = transcript_len / *num_bins;
            let bins: Vec<std::ops::Range<u32>> = (0..*num_bins).map(|i| i * bin_size..(i + 1) * bin_size).collect();

            let dr_threshold = (bins[0].end - bins[0].start)/2;

            
            let mut discarded_read_flag: bool = true;
            let read_start = a.alignment_start().unwrap().get() as u32;
            let read_end = a.alignment_end().unwrap().get() as u32;
        
            for (i, bin) in bins.iter().enumerate(){
                if read_start <= bin.start && read_end >= bin.end{
                    discarded_read_flag = false;
                } 
                else if read_start >= bin.start && read_end >= bin.start && read_start <=bin.end && read_end <=bin.end  {
                    if (read_end - read_start) >= dr_threshold {
                        discarded_read_flag = false;
                    }
                }
                else if read_start >= bin.start && read_start <bin.end  {
                    if (bin.end - read_start) >= dr_threshold {
                        discarded_read_flag = false;
                    }
                }
                else if read_end > bin.start && read_end <=bin.end  {
                    if (read_end - bin.start) >= dr_threshold {
                        discarded_read_flag = false;
                    }
                }
            }

            if discarded_read_flag {
                num_discarded_alignment = num_discarded_alignment + 1;
            }
        }

        if num_discarded_alignment == alns.len(){

            num_discarded_reads = num_discarded_reads + 1;
        }
    }


    num_discarded_reads

}



fn main() -> io::Result<()> {

    let args = Args::parse();

    let filter_opts = variables::AlignmentFilters::builder()
        .five_prime_clip(args.five_prime_clip)
        .three_prime_clip(args.three_prime_clip)
        .score_threshold(args.score_threshold)
        .min_aligned_fraction(args.min_aligned_fraction)
        .min_aligned_len(args.min_aligned_len)
        .allow_rc(args.allow_negative_strand)
        .build();

    //check if the "no_coverage", "prob", "bin" argumnets exist or not and if they do not exist provide the intial values for them
    let mut coverage = true; 
    if args.coverage == "no"{
        coverage = false;
    }
    //the prob can have one of these values [uniform, DisPoisson, ConPoisson, multinomial]
    let prob = args.prob.as_str(); 
    //the rate could be one of these values [shr, shr_tlen, dr]
    let rate = args.rate.as_str(); 
    let bins: u32 = args.bins; 
    

    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(bam::Reader::new)?;

    let header = reader.read_header()?;

    for (prog, _pmap) in header.programs().iter() {
        eprintln!("program: {}", prog);
    }

    let mut txps: Vec<variables::TranscriptInfo> = Vec::with_capacity(header.reference_sequences().len());

    // loop over the transcripts in the header and fill in the relevant
    // information here.
    let mut txps_name: Vec<String> = Vec::new();
    for (rseq, rmap) in header.reference_sequences().iter() {
        txps_name.push(rseq.to_string());
        txps.push(variables::TranscriptInfo::with_len(rmap.length()));
    }

    //let mut rmap = HashMap<usize, ::new();
    //
    let mut prev_read = String::new();
    let mut num_mapped = 0_u64;
    let mut records_for_read = vec![];
    let mut store: variables::InMemoryAlignmentStore = variables::InMemoryAlignmentStore::new(filter_opts);

    for result in reader.records(&header) {
        //eprintln!("in the first for loop");
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        let record_copy = record.clone();
        if let Some(rname) = record.read_name() {
            let rstring: String =
                <noodles_sam::record::read_name::ReadName as AsRef<str>>::as_ref(rname).to_owned();
            //eprintln!("rstring: {:?}", rstring);
            // if this is an alignment for the same read, then 
            // push it onto our temporary vector.
            if prev_read == rstring {
                if let (Some(ref_id), false) = (
                    record.reference_sequence_id(),
                    record.flags().is_supplementary(),
                ) {
                    records_for_read.push(record_copy);
                    //txps[ref_id].ranges.push(
                    //    (record.alignment_start().unwrap().get() as u32)
                    //        ..(record.alignment_end().unwrap().get() as u32),
                    //);
                }
            } else {
                if !prev_read.is_empty() {
                    //println!("the previous read had {} mappings", records_for_read.len());
                    store.add_group(&mut txps,&mut records_for_read);
                    records_for_read.clear();
                    num_mapped += 1;
                }
                prev_read = rstring;
                if let (Some(ref_id), false) = (
                    record.reference_sequence_id(),
                    record.flags().is_supplementary(),
                ) {
                    records_for_read.push(record_copy);
                    //txps[ref_id].ranges.push(
                    //    (record.alignment_start().unwrap().get() as u32)
                    //        ..(record.alignment_end().unwrap().get() as u32),
                    //);
                }
            }
        }
    }
    if !records_for_read.is_empty() {
        store.add_group(&mut txps,&mut records_for_read);
        records_for_read.clear();
        num_mapped += 1;
    }

    eprintln!("Number of mapped reads : {}", num_mapped);
    eprintln!("normalizing alignment scores");
    store.normalize_scores();
    eprintln!("Total number of alignment records : {}", store.total_len());
    eprintln!("number of aligned reads : {}", store.num_aligned_reads());
    eprintln!("computing coverages");

    //defining the global threading
    //rayon::ThreadPoolBuilder::new()
    //.num_threads(args.threads)
    //.build_global()
    //.unwrap();
    
    let num_alignments: Vec<usize>;
    let num_discarded_alignments: Vec<usize>;
    let num_discarded_reads_decision_rule: usize;

    if !coverage {
        num_alignments = nocoverage_prob::no_coverage_prob(&mut txps, args.threads);
        num_discarded_alignments = vec![0; txps.len()];
        num_discarded_reads_decision_rule = 0;

    }else{
        match prob {
            "uniform" => {
                num_alignments = uniform_prob::uniform_prob(&mut txps, args.threads);
                num_discarded_alignments = vec![0; txps.len()];
                num_discarded_reads_decision_rule = 0;
            }
            "DisPoisson" => {
                (num_alignments, num_discarded_alignments) = poisson_prob::discrete_poisson_prob(&mut txps, rate, &bins, args.threads);
                num_discarded_reads_decision_rule = num_discarded_reads_fun(&store, &txps, &bins);
            }
            "ConPoisson" => {
                (num_alignments, num_discarded_alignments) = poisson_prob::continuous_poisson_prob(&mut txps, rate, &bins, args.threads);
                num_discarded_reads_decision_rule = 0;
            }
            "multinomial" => {
                (num_alignments, num_discarded_alignments) = multinomial_prob::multinomial_prob(&mut txps, rate, &bins, args.threads);
                num_discarded_reads_decision_rule = num_discarded_reads_fun(&store, &txps, &bins);
                //num_discarded_reads_decision_rule = 0;
            }
            "LW" => {
                (num_alignments, num_discarded_alignments) = lander_waterman_prob::lander_waterman_prob(&mut txps, rate, &bins, args.threads);
                num_discarded_reads_decision_rule = if rate == "dr" {num_discarded_reads_fun(&store, &txps, &bins)} else {0};
            }
            "binomial" => {
                (num_alignments, num_discarded_alignments) = binomial_continuous::binomial_continuous_prob(&mut txps, rate, &bins, args.threads);
                num_discarded_reads_decision_rule = 0;
            }
            "entropy" => {
                (num_alignments, num_discarded_alignments) = entropy_distribution::entropy_prob(&mut txps, rate, &bins, args.threads, args.alpha, args.beta); //, &args.output_path, &txps_name);
                num_discarded_reads_decision_rule = 0;
            }
            //"kde" => {
            //    let start_time = Instant::now();
            //    num_alignments = kde_prob::kde_prob(&mut txps, rate, args.threads);
            //    let end_time = Instant::now();
            //    let elapsed_time = end_time.duration_since(start_time);
            //    eprintln!("Time taken for kde: {:?}", elapsed_time);
//
            //    num_discarded_alignments = vec![0; txps.len()]; 
            //    num_discarded_reads_decision_rule = 0;
            //}
            //"kde_c" => {
            //    let start_time = Instant::now();
            //    num_alignments = kde_prob::kde_c_prob(&mut txps, rate, args.threads);
            //    let end_time = Instant::now();
            //    let elapsed_time = end_time.duration_since(start_time);
            //    eprintln!("Time taken for kde: {:?}", elapsed_time);
//
            //    num_discarded_alignments = vec![0; txps.len()]; 
            //    num_discarded_reads_decision_rule = 0;
            //}
            _ => {
                panic!("{:?} probability function is not defined here!", prob);
            }
        }

    }

    //["multinomial", "LW"]
    eprintln!("done");
    eprintln!("normalizing read probabilities");

    normalize_read_probs(&mut store, &mut txps, coverage, &prob, &rate);
    eprintln!("done");

    let emi = variables::EMInfo {
        eq_map: &store,
        txp_info: &txps,
        max_iter: 1000,
    };

    eprintln!("computing estimated counts in the EM");
    let (counts, num_discarded_reads_em, likelihood)  = em(&emi, args.short_quant, &txps_name);
    //let (current_vec, initial_vec, likelihood_vec) = em(&emi, args.short_quant, &txps_name);
    eprintln!("done");

    //write the stat output
    eprintln!("write output");
    let num_reads = store.num_aligned_reads();
    //common_functions::write_out_probs(&args.output_path, prob, rate, &bins, args.alpha, args.beta, &emi, &txps_name).expect("Failed to write probs output");
    //common_functions::write_out_cdf(&args.output_path, prob, rate, &bins, args.alpha, args.beta, &emi, &txps_name).expect("Failed to write probs output");
    //common_functions::write_out_stat(&args.output_path, prob, rate, &bins, args.alpha, args.beta, &header, &num_alignments, &num_discarded_alignments, &num_reads, &num_discarded_reads_decision_rule, &num_discarded_reads_em, &counts).expect("Failed to write stat output");
    common_functions::write_out_count(&args.output_path, prob, rate, &bins, args.alpha, args.beta, &header, &txps, &counts).expect("Failed to write count output");
    //common_functions::write_out_initial_count_likelihood(&args.output_path,  &bins, &header, &current_vec, &initial_vec, &likelihood_vec).expect("Failed to write count output");
    //eprintln!("likelihood is: {:?}", likelihood);

    Ok(())
}


