use clap::Parser;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader},
    num::NonZeroUsize,
};

use bio_types::annot::loc::Loc;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::Strand;
use noodles_bam as bam;
use noodles_gtf as gtf;
use noodles_gtf::record::Strand as NoodlesStrand;
use noodles_sam as sam;
use sam::record::data::field::tag;

use bio_types::annot::contig::Contig;
use coitrees::{COITree, IntervalNode};
use nested_intervals::IntervalSet;
use special::Gamma;

//use statrs::distribution::{Multinomial, Discrete};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser)]
    alignments: String,
    coverage: String,
    prob: String,
    rate: String,
    bins: u32,
}

#[derive(Debug, PartialEq)]
struct TranscriptInfo {
    len: NonZeroUsize,
    ranges: Vec<std::ops::Range<u32>>,
    coverage: f64,
}

impl TranscriptInfo {
    fn new() -> Self {
        Self {
            len: NonZeroUsize::new(0).unwrap(),
            ranges: Vec::new(),
            coverage: 0.0,
        }
    }
    fn with_len(len: NonZeroUsize) -> Self {
        Self {
            len,
            ranges: Vec::new(),
            coverage: 0.0,
        }
    }
}

#[derive(Debug)]
struct InMemoryAlignmentStore {
    alignments: Vec<sam::alignment::record::Record>,
    probabilities: Vec<f32>,
    // holds the boundaries between records for different reads
    boundaries: Vec<usize>,
}

struct InMemoryAlignmentStoreIter<'a> {
    store: &'a InMemoryAlignmentStore,
    idx: usize,
}

impl<'a> Iterator for InMemoryAlignmentStoreIter<'a> {
    type Item = (&'a [sam::alignment::record::Record], &'a [f32]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx + 1 >= self.store.boundaries.len() {
            None
        } else {
            let start = self.store.boundaries[self.idx];
            let end = self.store.boundaries[self.idx + 1];
            self.idx += 1;
            Some((
                &self.store.alignments[start..end],
                &self.store.probabilities[start..end],
            ))
        }
    }
}

impl InMemoryAlignmentStore {
    fn new() -> Self {
        InMemoryAlignmentStore {
            alignments: vec![],
            probabilities: vec![],
            boundaries: vec![0],
        }
    }

    fn iter(&self) -> InMemoryAlignmentStoreIter {
        InMemoryAlignmentStoreIter {
            store: &self,
            idx: 0,
        }
    }

    fn add_group(&mut self, ag: &Vec<sam::alignment::record::Record>) {
        self.alignments.extend_from_slice(&ag);
        self.boundaries.push(self.alignments.len());
    }

    fn total_len(&self) -> usize {
        self.alignments.len()
    }

    fn num_aligned_reads(&self) -> usize {
        if self.boundaries.len() > 0 {
            self.boundaries.len() - 1
        } else {
            0
        }
    }

    fn normalize_scores(&mut self) {
        self.probabilities = vec![0.0_f32; self.alignments.len()];
        for w in self.boundaries.windows(2) {
            let s: usize = w[0];
            let e: usize = w[1];
            if e - s > 1 {
                let mut max_score = 0_i32;
                let mut scores = Vec::<i32>::with_capacity(e - s);
                for a in &self.alignments[s..e] {
                    let score_value = a
                        .data()
                        .get(&tag::ALIGNMENT_SCORE)
                        .expect("could not get value");
                    let score = score_value.as_int().unwrap() as i32;
                    scores.push(score);
                    if score > max_score {
                        max_score = score;
                    }
                }
                for (i, score) in scores.iter().enumerate() {
                    let f = ((*score as f32) - (max_score as f32)) / 10.0_f32;
                    self.probabilities[s + i] = f.exp();
                }
            } else {
                self.probabilities[s] = 1.0
            }
        }
    }
}

/// Holds the info relevant for running the EM algorithm
struct EMInfo<'eqm, 'tinfo> {
    eq_map: &'eqm InMemoryAlignmentStore,
    txp_info: &'tinfo Vec<TranscriptInfo>,
    max_iter: u32,
}

#[inline]
fn m_step(
    eq_map: &InMemoryAlignmentStore,
    tinfo: &[TranscriptInfo],
    read_coverage_prob: &Vec<Vec<f64>>,
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) {
    for ((alns, probs), coverage_probs) in eq_map.iter().zip(read_coverage_prob.iter()) {
        let mut denom = 0.0_f64;
        for (a, (p, cp)) in alns.iter().zip(probs.iter().zip(coverage_probs.iter())) {
            let target_id = a.reference_sequence_id().unwrap();
            let prob = *p as f64;
            //let cov_prob = if tinfo[target_id].coverage < 0.5 { 1e-5 } else { 1.0 };//powf(2.0) as f64;
            //let cov_prob = tinfo[target_id].coverage;//powf(2.0) as f64;
            let cov_prob = *cp;
            //let cov_prob = 1.0;
            denom += prev_count[target_id] * prob * cov_prob;
        }

        if denom > 1e-8 {
            for (a, (p,cp)) in alns.iter().zip(probs.iter().zip(coverage_probs.iter())) {
                let target_id = a.reference_sequence_id().unwrap();
                let prob = *p as f64;
                //let cov_prob = if tinfo[target_id].coverage < 0.5 { 1e-5 } else { 1.0 };//powf(2.0) as f64;
                //let cov_prob = tinfo[target_id].coverage;
                let cov_prob = *cp;
                //let cov_prob = 1.0;
                curr_counts[target_id] += (prev_count[target_id] * prob * cov_prob) / denom;
            }
        }
    }
}

fn em(em_info: &EMInfo, read_coverage_prob: &Vec<Vec<f64>>) -> Vec<f64> {
    let eq_map = em_info.eq_map;
    let tinfo = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

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
    m_step(eq_map, tinfo, read_coverage_prob, &mut prev_counts, &mut curr_counts);

    curr_counts
}


fn bin_transcript_decision_rule(t: &TranscriptInfo, num_bins: &u32) -> (Vec<u32>, Vec<u32>) {

    let transcript_len = t.len.get() as u32; //transcript length
    let bin_size = transcript_len / *num_bins;
    let bins: Vec<std::ops::Range<u32>> = (0..*num_bins)
        .map(|i| i * bin_size..(i + 1) * bin_size)
        .collect();

    let mut bin_counts: Vec<u32> = vec![0; bins.len()];
    let bin_lengths: Vec<u32> = vec![bins[0].end -bins[0].start; bins.len()];

    let dr_threshold = (bins[0].end - bins[0].start)/2;

    for read in t.ranges.iter(){
        for (i, bin) in bins.iter().enumerate(){
            if read.start <= bin.start && read.end >= bin.end{
                bin_counts[i] = bin_counts[i] + 1;
            } 
            else if read.start >= bin.start && read.end >= bin.start && read.start <=bin.end && read.end <=bin.end  {
                if (read.end - read.start) >= dr_threshold {
                    bin_counts[i] = bin_counts[i] + 1;
                }
            }
            else if read.start >= bin.start && read.start <bin.end  {
                if (bin.end - read.start) >= dr_threshold {
                    bin_counts[i] = bin_counts[i] + 1;
                }
            }
            else if read.end > bin.start && read.end <=bin.end  {
                if (read.end - bin.start) >= dr_threshold {
                    bin_counts[i] = bin_counts[i] + 1;
                }
            }
            
        }
    }

 (bin_counts, bin_lengths)

}


fn bin_transcript_normalize_counts(t: &TranscriptInfo, num_bins: &u32) -> (Vec<u32>, Vec<u32>) {

    let transcript_len = t.len.get() as u32; //transcript length
    let bin_size = transcript_len / *num_bins;
    let bins: Vec<std::ops::Range<u32>> = (0..*num_bins)
        .map(|i| i * bin_size..(i + 1) * bin_size)
        .collect();

    let mut bin_counts: Vec<u32> = vec![0; bins.len()];
    let bin_lengths: Vec<u32> = vec![bins[0].end -bins[0].start; bins.len()];

    for read in t.ranges.iter(){
        for (i, bin) in bins.iter().enumerate(){
            if read.start <= bin.start && read.end >= bin.end{
                bin_counts[i] = bin_counts[i] + 1;
            } 
            else if read.start >= bin.start && read.end >= bin.start && read.start <=bin.end && read.end <=bin.end  {
                    bin_counts[i] = bin_counts[i] + ((read.end - read.start) / (bin.end - bin.start));
            }
            else if read.start >= bin.start && read.start <bin.end  {
                    bin_counts[i] = bin_counts[i] + ((bin.end - read.start) / (bin.end - bin.start));
            }
            else if read.end > bin.start && read.end <=bin.end  {
                    bin_counts[i] = bin_counts[i] + ((read.end - bin.start) / (bin.end - bin.start));
                }
        }
            
    }

 (bin_counts, bin_lengths)

}

//https://en.wikipedia.org/wiki/DNA_sequencing_theory#Early_uses_derived_from_elementary_probability_theory
fn lander_waterman(n: u32, l:u32, g:u32) -> f64{

    1.0 - (-((n * l) as f64/g as f64)).exp()

}

//https://en.wikipedia.org/wiki/Poisson_distribution#Computational_methods
fn continuous_poisson_probability(k: u32, t: f64, r: f64) -> f64{

    if t < 0.0 {
        panic!("Error: Invalid input. Interval length must be a positive value.");
    }

    if r < 0.0 {
        panic!("Error: Invalid input. Rate must be a positive value.");
    }

    if r == 0.0 || t == 0.0 {
        return 0.0;
    }

    let rate = r * t;
    let log_numerator = (k as f64 * rate.ln()) - rate;
    let log_denominator = ((k as f64 + 1.0).gamma()).ln();

    let result = (log_numerator - log_denominator).exp();

    if result.is_nan() || result.is_infinite() {
        
        panic!("Error: Invalid result. Poisson Probability value is NaN or infinity.");
    }
    
    result

}

fn poisson_probability(k: u32, t: f64, r: f64) -> f64{

    if t < 0.0 {
        panic!("Error: Invalid input. Interval length must be a positive value.");
    }

    if r < 0.0 {
        panic!("Error: Invalid input. Rate must be a positive value.");
    }

    if r == 0.0 || t == 0.0 {
        return 0.0;
    }

    let log_numerator1 = (r * t).ln() * k as f64;
    let log_numerator2 = (-r * t).exp().ln();
    let log_denominator = factorial_ln(k) as f64;
    
    let result = (log_numerator1 + log_numerator2 - log_denominator).exp();
    
    if result.is_nan() || result.is_infinite() {
        
        panic!("Error: Invalid result. Poisson Probability value is NaN or infinity.");
    }
    
    result
    
}

fn factorial_ln(n: u32) -> f64 {
    if n == 0 {
        0.0
    } else {
        (1..=n).map(|x| (x as f64).ln()).sum()
    }
}

fn multinomial_probability(interval_count: Vec<u32>, interval_length: Vec<u32>, distinct_rate: f64) -> f64{

    let mut interval_counts = interval_count;
    let interval_lengths = interval_length;
    if interval_counts.iter().sum::<u32>() == 0 {
        return 0.0;
    }

    if distinct_rate == 0.0 {
        return 0.0;
    }

    
    let mut probabilities: Vec<f64> = interval_counts
                                  .iter()
                                  .zip(interval_lengths.iter())
                                  .map(|(&count, &length)| {
                                  if count ==0 || length == 0 {
                                      0.0 
                                  } else {
                                      (count as f64) / (length as f64 * distinct_rate)
                                  }
                                  }).collect();

    //obtain the indices of the probabilities that have the zero values
    let mut zero_value_indices: Vec<usize> = Vec::new();
    for (i, &prob) in probabilities.iter().enumerate(){
        if prob == 0.0{
            zero_value_indices.push(i);
        }
    }
    //remove the zero_valu_indices from the probabilities and interval_counts
    for &i in zero_value_indices.iter().rev(){
        probabilities.remove(i);
        interval_counts.remove(i);
    }

    //check if the sum of probabilities are equal to 1
    let probabilities_sum: f64 = probabilities.iter().sum();
    if probabilities_sum < 0.9999 || probabilities_sum > 1.0001 {
        panic!("Invalid probabilities: Sum is not equal to 1.0");
    }

    //check if the length of the probabilities is equal to the length of interval_counts
    if probabilities.len() != interval_counts.len(){
        panic!("the length of probabilities and interval_counts are not equal!");
    }

    //let multinomial = Multinomial::new(probabilities.as_slice(), interval_counts.iter().sum::<u32>() as u64).expect("Invalid probabilities");
    //let prob_dr_multinomial: f64 = multinomial.ln_pmf(&interval_counts.iter().map(|&count| count as u64).collect::<Vec<u64>>());
    //let result = prob_dr_multinomial.exp();

    let log_numerator1: f64 = factorial_ln(interval_counts.iter().sum::<u32>()); 
    let log_denominator: f64 =  interval_counts.iter().map(|&count| factorial_ln(count)).sum();//factorial_ln(k) as f64;
    let log_numerator2: f64 = probabilities.iter().zip(interval_counts.iter()).map(|(&prob, &count)| prob.ln() * (count as f64)).sum();
    
    let mut result = (log_numerator1 - log_denominator + log_numerator2 ).exp();

    if result.is_nan() || result.is_infinite() {

        panic!("Error: Invalid result. Poisson Probability value is NaN or infinity.");
    }
    
    if (result > 1.0) && (result < 1.0001){

        result = 1.0_f64;
    }
    
    result    
    
}


fn normalize_read_probs(em_info: &EMInfo, coverage: bool) -> Vec<Vec<f64>>{

    let eq_map = em_info.eq_map;
    let tinfo = em_info.txp_info;
    let mut normalize_probs_vec: Vec<Vec<f64>> = Vec::new();
    let mut normalize_probs_temp: Vec<f64> = vec![];

    for (alns, _probs) in eq_map.iter() {
        let mut normalized_prob_section: Vec<f64> = vec![];
        
        if coverage {
            for a in alns.iter() {

                let target_id = a.reference_sequence_id().unwrap();
                let cov_prob = tinfo[target_id].coverage;
                normalize_probs_temp.push(cov_prob);
            }  
            let sum_normalize_probs_temp: f64 = normalize_probs_temp.iter().sum();    
            normalized_prob_section = normalize_probs_temp.iter().map(|&prob| prob/sum_normalize_probs_temp).collect();
        } else {
            normalized_prob_section = vec![1.0 as f64; alns.len()];
        }
        
        normalize_probs_vec.push(normalized_prob_section);
        normalize_probs_temp.clear();
    }

    //panic!("end of normalize function");
    normalize_probs_vec
}


fn uniform_prob(t: &TranscriptInfo) -> f64 {

    let len = t.len.get() as u32; //transcript length
    let interval_set = IntervalSet::new(&t.ranges).expect("couldn't build interval set");
    let mut interval_set = interval_set.merge_connected();
    let covered = interval_set.covered_units();
    let prob_uniform = (covered as f64) / (len as f64);

    prob_uniform
}


fn discrete_poisson_prob(t: &TranscriptInfo, rate: &str, bins: &u32) -> f64{

    if *bins != 0 {
        let bin_counts: Vec<u32>;
        let bin_lengths: Vec<u32>;
        (bin_counts, bin_lengths) = bin_transcript_decision_rule(t, bins); //binning the transcript length and obtain the counts and length vectors
        match rate {
            "shr" | "shr_tlen" => {
                let shared_rate: f64 = (bin_counts.iter().sum::<u32>() as f64) / (bin_lengths.iter().sum::<u32>() as f64);
                let prob_shr: f64 = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, shared_rate)).product();
                return prob_shr;

            }
            "dr" => {
                let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: f64 = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, distinct_rate)).product();
                return prob_dr;
            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }
    }else {
        //not binning the transcript length
        let len = t.len.get() as u32; //transcript length

        //obtain the start and end of reads aligned to these transcripts
        let mut start_end_ranges: Vec<u32> = t.ranges.iter().map(|range| vec![range.start, range.end]).flatten().collect();
        start_end_ranges.push(0); // push the first position of the transcript
        start_end_ranges.push(len); // push the last position of the transcript
        start_end_ranges.sort(); // Sort the vector in ascending order
        start_end_ranges.dedup(); // Remove consecutive duplicates
        //convert the sorted vector of starts and ends into a vector of consecutive ranges
        let distinct_interval: Vec<std::ops::Range<u32>> = start_end_ranges.windows(2).map(|window| window[0]..window[1]).collect(); 
        let interval_length : Vec<u32> = start_end_ranges.windows(2).map(|window| window[1] - window[0]).collect(); 
        //obtain the number of reads aligned in each distinct intervals
        let mut interval_counts: Vec<u32> = Vec::new();
        for interval in distinct_interval
        {
            interval_counts.push(t.ranges.iter().filter(|range| range.start <= interval.start && range.end >= interval.end).count() as u32);
        }

        match rate {
            "shr" => {
                let shared_rate: f64 = (interval_counts.iter().sum::<u32>() as f64) / (interval_length.iter().sum::<u32>() as f64);
                let prob_shr: f64 = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, shared_rate)).product();
                return prob_shr;

            }
            "shr_tlen" => {
                let shared_rate_tlen: f64 = (interval_counts.iter().sum::<u32>() as f64) / (len as f64);
                let prob_shr_tlen: f64 = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, shared_rate_tlen)).product();
                return prob_shr_tlen;
            }
            "dr" => {
                let distinct_rate: f64 =  interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: f64 = poisson_probability(interval_counts.iter().sum(), len as f64, distinct_rate);
                return prob_dr;
            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }
    }
    


}


fn continuous_poisson_prob(t: &TranscriptInfo, rate: &str, bins: &u32) -> f64{

    if *bins != 0 {
        let bin_counts: Vec<u32>;
        let bin_lengths: Vec<u32>;
        (bin_counts, bin_lengths) = bin_transcript_normalize_counts(t, bins); //binning the transcript length and obtain the counts and length vectors
        match rate {
            "shr" | "shr_tlen" => {
                let shared_rate: f64 = (bin_counts.iter().sum::<u32>() as f64) / (bin_lengths.iter().sum::<u32>() as f64);
                let prob_shr: f64 = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, shared_rate)).product();
                return prob_shr;

            }
            "dr" => {
                let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: f64 = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, distinct_rate)).product();
                return prob_dr;
            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }
    }else {
        //not binning the transcript length
        let len = t.len.get() as u32; //transcript length

        //obtain the start and end of reads aligned to these transcripts
        let mut start_end_ranges: Vec<u32> = t.ranges.iter().map(|range| vec![range.start, range.end]).flatten().collect();
        start_end_ranges.push(0); // push the first position of the transcript
        start_end_ranges.push(len); // push the last position of the transcript
        start_end_ranges.sort(); // Sort the vector in ascending order
        start_end_ranges.dedup(); // Remove consecutive duplicates
        //convert the sorted vector of starts and ends into a vector of consecutive ranges
        let distinct_interval: Vec<std::ops::Range<u32>> = start_end_ranges.windows(2).map(|window| window[0]..window[1]).collect(); 
        let interval_length : Vec<u32> = start_end_ranges.windows(2).map(|window| window[1] - window[0]).collect(); 
        //obtain the number of reads aligned in each distinct intervals
        let mut interval_counts: Vec<u32> = Vec::new();
        for interval in distinct_interval
        {
            interval_counts.push(t.ranges.iter().filter(|range| range.start <= interval.start && range.end >= interval.end).count() as u32);
        }

        match rate {
            "shr" => {
                let shared_rate: f64 = (interval_counts.iter().sum::<u32>() as f64) / (interval_length.iter().sum::<u32>() as f64);
                let prob_shr: f64 = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, shared_rate)).product();
                return prob_shr;

            }
            "shr_tlen" => {
                let shared_rate_tlen: f64 = (interval_counts.iter().sum::<u32>() as f64) / (len as f64);
                let prob_shr_tlen: f64 = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, shared_rate_tlen)).product();
                return prob_shr_tlen;
            }
            "dr" => {
                let distinct_rate: f64 =  interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: f64 = continuous_poisson_probability(interval_counts.iter().sum(), len as f64, distinct_rate);
                return prob_dr;
            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }
    }
    


}


fn multinomial_prob(t: &TranscriptInfo, rate: &str, bins: &u32) -> f64{

    if *bins != 0 {
        let bin_counts: Vec<u32>;
        let bin_lengths: Vec<u32>;
        (bin_counts, bin_lengths) = bin_transcript_decision_rule(t, bins); //binning the transcript length and obtain the counts and length vectors
        match rate {
            "shr" | "shr_tlen" => {
                panic!("Shared Rate cannot be used for Multinomial model");
            }
            "dr" => {
                let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: f64 = multinomial_probability(bin_counts, bin_lengths, distinct_rate);
                return prob_dr;
            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }
    }else {
        //not binning the transcript length
        let len = t.len.get() as u32; //transcript length

        //obtain the start and end of reads aligned to these transcripts
        let mut start_end_ranges: Vec<u32> = t.ranges.iter().map(|range| vec![range.start, range.end]).flatten().collect();
        start_end_ranges.push(0); // push the first position of the transcript
        start_end_ranges.push(len); // push the last position of the transcript
        start_end_ranges.sort(); // Sort the vector in ascending order
        start_end_ranges.dedup(); // Remove consecutive duplicates
        //convert the sorted vector of starts and ends into a vector of consecutive ranges
        let distinct_interval: Vec<std::ops::Range<u32>> = start_end_ranges.windows(2).map(|window| window[0]..window[1]).collect(); 
        let interval_length : Vec<u32> = start_end_ranges.windows(2).map(|window| window[1] - window[0]).collect(); 
        //obtain the number of reads aligned in each distinct intervals
        let mut interval_counts: Vec<u32> = Vec::new();
        for interval in distinct_interval
        {
            interval_counts.push(t.ranges.iter().filter(|range| range.start <= interval.start && range.end >= interval.end).count() as u32);
        }

        match rate {
            "shr" | "shr_tlen" => {
                panic!("Shared Rate cannot be used for Multinomial model");
            }
            "dr" => {
                let distinct_rate: f64 =  interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: f64 = multinomial_probability(interval_counts, interval_length, distinct_rate);
                return prob_dr;
            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }
    }
    


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
    //the rate could be one of these values [shr, shr_tlen, disr]
    let rate = args.rate.as_str(); 
    let bins: u32 = args.bins; 


    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(bam::Reader::new)?;

    let header = reader.read_header()?;

    for (prog, _pmap) in header.programs().iter() {
        eprintln!("program: {}", prog);
    }

    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(header.reference_sequences().len());

    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (_rseq, rmap) in header.reference_sequences().iter() {
        // println!("ref: {}, rmap : {:?}", rseq, rmap.length());
        txps.push(TranscriptInfo::with_len(rmap.length()));
    }

    //let mut rmap = HashMap<usize, ::new();
    //
    let mut prev_read = String::new();
    let mut num_mapped = 0_u64;
    let mut records_for_read = vec![];
    let mut store = InMemoryAlignmentStore::new();

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

    eprintln!("computing coverages");
    for t in txps.iter_mut() {

        if !coverage {
            t.coverage = 1.0 as f64;

        }else{
            match prob {
                "DisPoisson" => {
                    t.coverage = discrete_poisson_prob(&t, rate, &bins);
                }
                "ConPoisson" => {
                    t.coverage = continuous_poisson_prob(&t, rate, &bins);
                }
                "uniform" => {
                    t.coverage = uniform_prob(&t);
                }
                "multinomial" => {
                    t.coverage = multinomial_prob(&t, rate, &bins);
                }
                _ => {
                    panic!("{:?} probability function is not defined here!", prob);
                }
            }

        }
        

    }
    eprintln!("done");

    eprintln!("Number of mapped reads : {}", num_mapped);
    eprintln!("normalizing alignment scores");
    store.normalize_scores();
    eprintln!("Total number of alignment records : {}", store.total_len());
    eprintln!("number of aligned reads : {}", store.num_aligned_reads());

    let emi = EMInfo {
        eq_map: &store,
        txp_info: &txps,
        max_iter: 1000,
    };
    let read_coverage_probs: Vec<Vec<f64>> = normalize_read_probs(&emi, coverage);
    let counts = em(&emi, &read_coverage_probs);

    println!("tname\tcoverage\tlen\tnum_reads"); 
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (i, (_rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        println!("{}\t{}\t{}\t{}", _rseq, txps[i].coverage, rmap.length(), counts[i]);
    }

    Ok(())
}








//
// ignore anything below this line for now 
//

#[allow(unused)]
fn main_old() -> io::Result<()> {
    let args = Args::parse();

    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(gtf::Reader::new)?;
    let mut evec = Vec::new();
    let mut tvec = Vec::new();
    let mut tmap = HashMap::new();

    for result in reader.records() {
        let record = result?;
        match record.ty() {
            "exon" => {
                let s: isize = (usize::from(record.start()) as isize) - 1;
                let e: isize = usize::from(record.end()) as isize;
                let l: usize = (e - s).try_into().unwrap();
                let mut t = String::new();
                for e in record.attributes().iter() {
                    if e.key() == "transcript_id" {
                        t = e.value().to_owned();
                    }
                }

                let ni = tmap.len();
                let tid = *tmap.entry(t.clone()).or_insert(ni);

                // if this is what we just inserted
                if ni == tid {
                    tvec.push(Spliced::new(0, 1, 1, Strand::Forward));
                }

                let strand = match record.strand().unwrap() {
                    NoodlesStrand::Forward => Strand::Forward,
                    NoodlesStrand::Reverse => Strand::Reverse,
                };
                let c = Contig::new(tid, s, l, strand);
                evec.push(c);
            }
            "transcript" => {
                let mut t = String::new();
                for e in record.attributes().iter() {
                    if e.key() == "transcript_id" {
                        t = e.value().to_owned();
                    }
                }
                let ni = tmap.len();
                let tid = *tmap.entry(t.clone()).or_insert(ni);

                // if this is what we just inserted
                if ni == tid {
                    tvec.push(Spliced::new(0, 1, 1, Strand::Forward));
                }
            }
            _ => {}
        }
    }

    let mut txp_to_exon = HashMap::new();

    let mut l = 0;
    let mut max_len = 0;
    for (i, e) in evec.iter().enumerate() {
        let mut v = txp_to_exon.entry(e.refid()).or_insert(vec![]);
        v.push(i);
        l = v.len();
        if l > max_len {
            max_len = l;
        }
    }

    let mut txp_features: HashMap<usize, _> = HashMap::new();

    for (k, v) in txp_to_exon.iter_mut() {
        let strand = evec[v[0]].strand();
        v.sort_unstable_by_key(|x| evec[*x as usize].start());
        let s = evec[v[0]].start();
        let starts: Vec<usize> = v
            .iter()
            .map(|e| (evec[*e as usize].start() - s) as usize)
            .collect();
        let lens: Vec<usize> = v.iter().map(|e| evec[*e as usize].length()).collect();
        println!("lens = {:?}, starts = {:?}", lens, starts);
        txp_features.insert(
            **k,
            Spliced::with_lengths_starts(k, s, &lens, &starts, strand).unwrap(),
        );
    }

    let interval_vec: Vec<IntervalNode<usize, usize>> = evec
        .iter()
        .enumerate()
        .map(|(i, e)| {
            IntervalNode::new(
                e.start() as i32,
                (e.start() + e.length() as isize) as i32,
                i,
            )
        })
        .collect();
    let ct = COITree::new(interval_vec);

    println!("parsed {} exons", evec.len());
    println!("parsed {} transcripts", tvec.len());
    println!("max exon transcript had {} exons", max_len);
    Ok(())
}
