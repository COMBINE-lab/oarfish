extern crate libc;

mod variables;
mod common_functions;
mod nocoverage_prob;
mod uniform_prob;
mod poisson_prob;
mod lander_waterman_prob;
mod multinomial_prob;
mod kde_prob;

use clap::Parser;


use std::{
    fs::File,
    io::{self, BufReader},
    vec,
};

use noodles_bam as bam;

use std::time::Instant;




/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser, required = true)]
    alignments: String,
    #[clap(short, long, value_parser, required = true)]
    stat_output: String,
    #[clap(short, long, value_parser, required = true)]
    out_cov_prob: String,
    #[clap(short, long, value_parser, default_value_t = String::from("yes"))]
    coverage: String,
    #[clap(short, long, value_parser, default_value_t = String::from("multinomial"))]
    prob: String,
    #[clap(short, long, value_parser, default_value_t = String::from("dr"))]
    rate: String,
    #[clap(short, long, value_parser, default_value_t = 1)]
    bins: u32,
}



#[inline]
fn m_step(
    eq_map: &variables::InMemoryAlignmentStore,
    tinfo: &[variables::TranscriptInfo],
    read_coverage_prob: &Vec<Vec<f64>>,
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) -> usize {
    let mut dropped: usize=0;
    for ((alns, probs), coverage_probs) in eq_map.iter().zip(read_coverage_prob.iter()) {
        let mut denom = 0.0_f64;
        for (a, (p, cp)) in alns.iter().zip(probs.iter().zip(coverage_probs.iter())) {
            let target_id: usize = a.reference_sequence_id().unwrap();
            let prob = *p as f64;
            //let cov_prob = if tinfo[target_id].coverage < 0.5 { 1e-5 } else { 1.0 };//powf(2.0) as f64;
            //let cov_prob = tinfo[target_id].coverage;//powf(2.0) as f64;
            let cov_prob = *cp;
            //let cov_prob = 1.0;
            denom += prev_count[target_id] * prob * cov_prob;
        }

        if denom > 1e-30_f64 {
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
    dropped
    //eprintln!("dropped:{:?}", dropped);
}

fn em(em_info: &variables::EMInfo, read_coverage_prob: &Vec<Vec<f64>>) -> (Vec<f64>, usize) {
    let eq_map = em_info.eq_map;
    let tinfo = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;
    let mut num_dropped_reads: usize = 0;

    // init
    let avg = total_weight / (tinfo.len() as f64);
    let mut prev_counts = vec![avg; tinfo.len()];
    let mut curr_counts = vec![0.0f64; tinfo.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(eq_map, tinfo, read_coverage_prob, &mut prev_counts, &mut curr_counts);

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
    num_dropped_reads = m_step(eq_map, tinfo, read_coverage_prob, &mut prev_counts, &mut curr_counts);

    (curr_counts, num_dropped_reads)
}




fn normalize_read_probs(em_info: &variables::EMInfo, coverage: bool, prob: &str, rate: &str) -> Vec<Vec<f64>>{

    let eq_map: &variables::InMemoryAlignmentStore = em_info.eq_map;
    let tinfo: &Vec<variables::TranscriptInfo> = em_info.txp_info;
    let mut normalize_probs_vec: Vec<Vec<f64>> = Vec::new();
    let mut normalize_probs_temp: Vec<f64> = vec![];
    let mut kde_c_tinfo: Vec<usize> = vec![0 ; tinfo.len()];

    for (alns, _probs) in eq_map.iter() {
        let mut normalized_prob_section: Vec<f64> = vec![];
        
        if coverage {
            for a in alns.iter() {
                let target_id = a.reference_sequence_id().unwrap();
                let cov_prob: f64;
                if prob == "kde" || (prob == "kde_c" && rate == "1D") {
                    let start_aln = a.alignment_start().unwrap().get() as usize;
                    let end_aln = a.alignment_end().unwrap().get() as usize;
                    let coverage_prob: &Vec<f32> = &tinfo[target_id].coverage_prob;
                    //cov_prob = (tinfo[target_id].coverage_prob[end_aln] - tinfo[target_id].coverage_prob[start_aln]) as f64;
                    cov_prob = (get_coverag_prob(coverage_prob, end_aln) - get_coverag_prob(coverage_prob, start_aln)) as f64;

                } else if prob == "multinomial" || prob == "LW" {
                    //eprintln!("it is here");
                    let start_aln = a.alignment_start().unwrap().get() as usize;
                    let end_aln = a.alignment_end().unwrap().get() as usize;
                    //eprintln!("1: {:?}, 2: {:?}", tinfo[target_id].coverage_prob[end_aln], tinfo[target_id].coverage_prob[start_aln]);
                    cov_prob = (tinfo[target_id].coverage_prob[end_aln] - tinfo[target_id].coverage_prob[start_aln]) as f64;
                    
                } else if prob == "kde_c" && rate == "2D" {
                    let coverage_prob: &Vec<f32> = &tinfo[target_id].coverage_prob;
                    cov_prob = coverage_prob[kde_c_tinfo[target_id]] as f64;
                    kde_c_tinfo[target_id] += 1;
                    
                }else {
                    cov_prob = tinfo[target_id].coverage;
                } 
                
                if cov_prob.is_nan() || cov_prob.is_infinite() {
                    //eprintln!("prob is: {:?}", prob);
                    //eprintln!("cov_prob is: {:?}", cov_prob);
                    panic!("Error: Invalid result. normalize_read_probs function.");
                }
                normalize_probs_temp.push(cov_prob);
            }  
            let sum_normalize_probs_temp: f64 = if normalize_probs_temp.iter().sum::<f64>() > 0.0 {normalize_probs_temp.iter().sum()} else {1.0};    
            normalized_prob_section = normalize_probs_temp.iter().map(|&prob| prob/sum_normalize_probs_temp).collect();
            //normalized_prob_section = normalize_probs_temp.clone();
        } else {
            normalized_prob_section = vec![1.0 as f64; alns.len()];
        }
        
        normalize_probs_vec.push(normalized_prob_section);
        normalize_probs_temp.clear();
    }

    //panic!("end of normalize function");
    normalize_probs_vec
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

    let mut num_discarded_reads: usize = 0;


    for (alns, _probs) in store.iter() {
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
    for (_rseq, rmap) in header.reference_sequences().iter() {
        // println!("ref: {}, rmap : {:?}", rseq, rmap.length());
        txps.push(variables::TranscriptInfo::with_len(rmap.length()));
    }

    //let mut rmap = HashMap<usize, ::new();
    //
    let mut prev_read = String::new();
    let mut num_mapped = 0_u64;
    let mut records_for_read = vec![];
    let mut store: variables::InMemoryAlignmentStore = variables::InMemoryAlignmentStore::new();

    for result in reader.records(&header) {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        let record_copy = record.clone();
        if let Some(rname) = record.read_name() {
            let rstring: String =
                <noodles_sam::record::read_name::ReadName as AsRef<str>>::as_ref(rname).to_owned();
            // if this is an alignment for the same read, then 
            // push it onto our temporary vector.
            if prev_read == rstring {
                if let (Some(ref_id), false) = (
                    record.reference_sequence_id(),
                    record.flags().is_supplementary(),
                ) {
                    records_for_read.push(record_copy);
                    txps[ref_id].ranges.push(
                        (record.alignment_start().unwrap().get() as u32)
                            ..(record.alignment_end().unwrap().get() as u32),
                    );
                }
            } else {
                if !prev_read.is_empty() {
                    //println!("the previous read had {} mappings", records_for_read.len());
                    store.add_group(&records_for_read);
                    records_for_read.clear();
                    num_mapped += 1;
                }
                prev_read = rstring;
                if let (Some(ref_id), false) = (
                    record.reference_sequence_id(),
                    record.flags().is_supplementary(),
                ) {
                    records_for_read.push(record_copy);
                    txps[ref_id].ranges.push(
                        (record.alignment_start().unwrap().get() as u32)
                            ..(record.alignment_end().unwrap().get() as u32),
                    );
                }
            }
        }
    }
    if !records_for_read.is_empty() {
        store.add_group(&records_for_read);
        records_for_read.clear();
        num_mapped += 1;
    }

    eprintln!("Number of mapped reads : {}", num_mapped);
    eprintln!("normalizing alignment scores");
    store.normalize_scores();
    eprintln!("Total number of alignment records : {}", store.total_len());
    eprintln!("number of aligned reads : {}", store.num_aligned_reads());

    eprintln!("computing coverages");
    
    let num_alignments: Vec<usize>;
    let num_discarded_alignments: Vec<usize>;
    let num_discarded_reads_decision_rule: usize;

    if !coverage {
        (num_alignments, num_discarded_alignments) = nocoverage_prob::no_coverage_prob(&mut txps);
        num_discarded_reads_decision_rule = 0;

    }else{
        match prob {
            "uniform" => {
                (num_alignments, num_discarded_alignments) = uniform_prob::uniform_prob(&mut txps);
                num_discarded_reads_decision_rule = 0;
            }
            "DisPoisson" => {
                (num_alignments, num_discarded_alignments) = poisson_prob::discrete_poisson_prob(&mut txps, rate, &bins);
                num_discarded_reads_decision_rule = num_discarded_reads_fun(&store, &txps, &bins);
            }
            "ConPoisson" => {
                (num_alignments, num_discarded_alignments) = poisson_prob::continuous_poisson_prob(&mut txps, rate, &bins);
                num_discarded_reads_decision_rule = 0;
            }
            "multinomial" => {
                (num_alignments, num_discarded_alignments) = multinomial_prob::multinomial_prob(&mut txps, rate, &bins);
                num_discarded_reads_decision_rule = num_discarded_reads_fun(&store, &txps, &bins);
            }
            "LW" => {
                (num_alignments, num_discarded_alignments) = lander_waterman_prob::lander_waterman_prob(&mut txps, rate, &bins);
                num_discarded_reads_decision_rule = if rate == "dr" {num_discarded_reads_fun(&store, &txps, &bins)} else {0};
            }
            "kde" => {
                let start_time = Instant::now();
                num_alignments = kde_prob::kde_prob(&mut txps, rate);
                let end_time = Instant::now();
                let elapsed_time = end_time.duration_since(start_time);
                eprintln!("Time taken for kde: {:?}", elapsed_time);

                num_discarded_alignments = vec![0; txps.len()]; 
                num_discarded_reads_decision_rule = 0;
            }
            "kde_c" => {
                let start_time = Instant::now();
                num_alignments = kde_prob::kde_c_prob(&mut txps, rate);
                let end_time = Instant::now();
                let elapsed_time = end_time.duration_since(start_time);
                eprintln!("Time taken for kde: {:?}", elapsed_time);

                num_discarded_alignments = vec![0; txps.len()]; 
                num_discarded_reads_decision_rule = 0;
            }
            _ => {
                panic!("{:?} probability function is not defined here!", prob);
            }
        }

    }

    //["multinomial", "LW"]
    eprintln!("done");
     
    let emi = variables::EMInfo {
        eq_map: &store,
        txp_info: &txps,
        max_iter: 1000,
    };


    let read_coverage_probs: Vec<Vec<f64>> = normalize_read_probs(&emi, coverage, &prob, &rate);
    
    common_functions::write_read_coverage(args.out_cov_prob, &read_coverage_probs).expect("Failed to write output");
    let (counts, num_discarded_reads_em)  = em(&emi, &read_coverage_probs);

    //write the stat output
    let num_reads = store.num_aligned_reads();
    common_functions::write_output(args.stat_output, &header, &num_alignments, &num_discarded_alignments, &num_reads, &num_discarded_reads_decision_rule, &num_discarded_reads_em, &counts).expect("Failed to write output");

    println!("tname\tcoverage\tlen\tnum_reads"); 
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (i, (_rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        println!("{}\t{}\t{}\t{}", _rseq, txps[i].coverage, rmap.length(), counts[i]);
    }

    Ok(())
}


