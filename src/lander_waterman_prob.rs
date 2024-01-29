use crate::variables::TranscriptInfo;
use crate::common_functions::bin_transcript_normalize_counts;
use crate::common_functions::bin_transcript_decision_rule;
use rayon::prelude::*;
use std::sync::{RwLock, Arc};


//https://en.wikipedia.org/wiki/DNA_sequencing_theory#Early_uses_derived_from_elementary_probability_theory
pub fn lander_waterman(n: f64, coverage: f64) -> f64{

    let result = 1.0 - (-(n * coverage)).exp();
    //result = if result < 10e-8_f64 {10e8_f64} else {1.0 / result};

    if result.is_nan() || result.is_infinite() {
        //eprintln!("result is: {:?}", result);
        panic!("Error: Invalid result. Poisson Probability value is NaN or infinity.");
    }
    result
}

pub fn lander_waterman_prob(txps: &mut Vec<TranscriptInfo>, rate: &str, bins: &u32, threads: usize) -> (Vec<usize>, Vec<usize>) {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let mut num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));
    let mut num_discarded_reads = Arc::new(RwLock::new(vec![0; txps.len()]));

    txps.par_iter_mut().enumerate().for_each(|(i, t)| {
        
        let mut temp_prob: f64 = 0.0;

        if *bins != 0 {
            let mut num_discarded_read_temp: usize = 0;
            let bin_coverage: Vec<f64>;
            
            //==============================================================================================
            ////fill out the number of reads and discarded reads
            {
                let mut w1 =num_reads.write().unwrap();
                (*w1)[i] = t.ranges.len();
            }
            //==============================================================================================
            match rate {
                "shr" | "shr_tlen" => {
                    let bin_counts: Vec<f32>;
                    let bin_lengths: Vec<f32>;
                    (bin_counts, bin_lengths, num_discarded_read_temp, bin_coverage) = bin_transcript_normalize_counts(t, bins); //binning the transcript length and obtain the counts and length vectors
                    {
                        let mut w2 =num_discarded_reads.write().unwrap();
                        (*w2)[i] = num_discarded_read_temp;
                    }
                    let prob_shr: Vec<f64> = bin_counts.iter().zip(bin_coverage.iter()).map(|(&count, &coverage)| lander_waterman(count as f64, coverage)).collect();

                    let tlen = t.len.get(); //transcript length
                    let bin_length = bin_lengths[0];
                    
                    // Compute the sum of probabilities
                    let sum: f64 = prob_shr.iter().sum();
                    // Normalize the probabilities by dividing each element by the sum
                    let normalized_prob: Vec<f64> = prob_shr.iter().map(|&prob| if prob == 0.0 || sum == 0.0 {0.0} else {(prob.ln() - (bin_length as f64).ln() - sum.ln()).exp()}).collect();
                    
                    let mut prob_vec = vec![0.0; tlen + 1];
                    let mut bin_start = 0;
                    for (i, bin_length) in bin_lengths.iter().enumerate() {
                        let bin_end = bin_start + (*bin_length).floor() as usize;
                    
                        let start_index = bin_start;
                        let end_index = if i != (bin_lengths.len() - 1) {bin_end} else {tlen + 1};
                    
                        prob_vec[start_index as usize..end_index as usize]
                            .iter_mut()
                            .for_each(|v| *v = normalized_prob[i as usize]);
                    
                        bin_start = bin_end;
                    }
                    //temp_prob = prob_shr.iter().filter(|&x| *x != 0.0).fold(1.0, |acc, x| acc * x);
                    //t.coverage = temp_prob;
                    let cdf: Vec<f64> = prob_vec.iter().scan(0.0, |acc, &prob| {
                        *acc += prob;
                        Some(*acc)
                    }).collect();
                    t.coverage_prob = cdf.iter().map(|&x| x as f32).collect();

                }
                "dr" => {
                    let bin_counts: Vec<u32>;
                    let bin_lengths: Vec<f32>;
                    (bin_counts, bin_lengths, num_discarded_read_temp, bin_coverage) = bin_transcript_decision_rule(t, bins); //binning the transcript length and obtain the counts and length vectors
                    {
                        let mut w2 =num_discarded_reads.write().unwrap();
                        (*w2)[i] = num_discarded_read_temp;
                    }
                    let prob_dr: Vec<f64> = bin_counts.iter().zip(bin_coverage.iter()).map(|(&count, &coverage)| lander_waterman(count as f64, coverage)).collect();

                    let tlen = t.len.get(); //transcript length
                    let bin_length = bin_lengths[0];

                    // Compute the sum of probabilities
                    let sum: f64 = prob_dr.iter().sum();
                    // Normalize the probabilities by dividing each element by the sum
                    let normalized_prob: Vec<f64> = prob_dr.iter().map(|&prob| if prob == 0.0 || sum == 0.0 {0.0} else {(prob.ln() - (bin_length as f64).ln() - sum.ln()).exp()}).collect();

                    let mut prob_vec = vec![0.0; tlen + 1];
                    let mut bin_start = 0;
                    for (i, bin_length) in bin_lengths.iter().enumerate() {
                        let bin_end = bin_start + (*bin_length).floor() as usize;
                    
                        let start_index = bin_start;
                        let end_index = if i != (bin_lengths.len() - 1) {bin_end} else {tlen + 1};
                    
                        prob_vec[start_index as usize..end_index as usize]
                            .iter_mut()
                            .for_each(|v| *v = normalized_prob[i as usize]);
                    
                        bin_start = bin_end;
                    }

                    //temp_prob = prob_dr.iter().filter(|&x| *x != 0.0).fold(1.0, |acc, x| acc * x);
                    //t.coverage = temp_prob;
                    let cdf: Vec<f64> = prob_vec.iter().scan(0.0, |acc, &prob| {
                        *acc += prob;
                        Some(*acc)
                    }).collect();
                    t.coverage_prob = cdf.iter().map(|&x| x as f32).collect();
                }
                _ => {
                    panic!("{:?} rate is not defined in the program", rate);
                }
            }
        } else {
            panic!("the number of bins should be at least one for lander waterman probability!");
        }
    });

    let mut nr1: Vec<usize> = vec![];
    {
        nr1 = num_reads.read().unwrap().to_vec().clone();
    }
    let mut nr2: Vec<usize> = vec![];
    {
        nr2 = num_discarded_reads.read().unwrap().to_vec().clone();
    }
    (nr1, nr2) 

}