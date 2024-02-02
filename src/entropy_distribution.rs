use crate::variables::TranscriptInfo;
use crate::common_functions::{bin_transcript_decision_rule, bin_transcript_normalize_counts, all_points_counts};
use rayon::prelude::*;
use std::sync::{RwLock, Arc};
use std::f64;

fn entropy_function(probability: &Vec<f64>) -> f64 {
    let sum: f64 = probability.iter().sum();
    let normalized_prob: Vec<f64> = probability.iter().map(|&p| p / sum).collect();
    let entropy: f64  = normalized_prob.iter().map(|&p| if p < 10e-8 {10e-8 * (10e-8 as f64).log2()} else {p * (p).log2()}).sum();

    1.0_f64 + (entropy / (probability.len() as f64).log2())
}

pub fn multinomial_binomial_entropy(interval_count: Vec<f32>, interval_length: Vec<f32>, distinct_rate: f64, tlen: usize, num_reads: usize, alpha: f64, beta: f64) -> Vec<f64>{

    let interval_counts = interval_count;
    let interval_lengths = interval_length;

    if interval_counts.iter().sum::<f32>() == 0.0 {
        return vec![0.0; tlen + 1];
    }

    if distinct_rate == 0.0 {
        return vec![0.0; tlen + 1];
    }
    //eprintln!("counts: {:?}", interval_counts);
    //eprintln!("length: {:?}", interval_lengths);
    let mut probabilities: Vec<f64> = interval_counts
                                  .iter()
                                  .zip(interval_lengths.iter())
                                  .map(|(&count, &length)| {
                                  if count == 0.0 || length == 0.0 {
                                      0.0 
                                  } else {
                                      //eprintln!("count: {:?}, length: {:?}, rate:{:?}", count, length, distinct_rate);
                                      (count as f64) / (length as f64 * distinct_rate)
                                  }
                                  }).collect();


    let num_bins = interval_lengths.len() as u32;
    // Compute the sum of probabilities
    let sum: f64 = probabilities.iter().sum();
    // Normalize the probabilities by dividing each element by the sum
    let normalized_prob: Vec<f64> = probabilities.iter().zip(interval_lengths.iter()).map(|(&prob, &bin_length)| prob / (bin_length as f64 * sum)).collect();


    let mut prob_vec = vec![0.0; tlen + 1];
    let mut bin_start = 0;
    for i in 0..num_bins {
        //let bin_start = i * interval_lengths[i as usize];
        let start_index = bin_start;
        let bin_end = (i + 1) as f32 * interval_lengths[i as usize];
    
        //let start_index = if i == 0 { bin_start } else { bin_start + 1 };
        let end_index = if i + 1 == num_bins { (tlen + 1) as u32 } else { bin_end.floor() as u32 };
    
        prob_vec[start_index as usize..end_index as usize]
            .iter_mut()
            .for_each(|v| *v = normalized_prob[i as usize]);

        bin_start = end_index;
    }

    let sum_prob_vec: f64 = prob_vec.iter().sum();
    let p1: Vec<f64> = prob_vec.iter().map(|&prob| prob / (sum_prob_vec)).collect();
    let uniform_prob = vec![1.0 / (p1.len() as f64); p1.len()];
    let number_reads = num_reads; //t.ranges.len();
    let coefficient1 = ((-alpha) * (number_reads - 1) as f64).exp();
    let coefficient2 = ((-beta) * (1.0 - entropy_function(&p1))).exp();


    let final_prob: Vec<f64> = uniform_prob.iter().zip(p1.iter()).map(|(&uni_prob, &prob1)| (coefficient1 * uni_prob) + ((1.0 - coefficient1) * coefficient2 * prob1)).collect();

    // Compute the sum of probabilities
    let sum_final_prob: f64 = final_prob.iter().sum();
    // Normalize the probabilities by dividing each element by the sum
    let normalized_final_prob: Vec<f64> = final_prob.iter().map(|prob| prob / sum_final_prob).collect();

    let cdf: Vec<f64> = normalized_final_prob.iter().scan(0.0, |acc, &prob| {
    *acc += prob;
    Some(*acc)
    }).collect();

    
    cdf
    
}


pub fn one_bin_entropy(interval_count: Vec<f32>, tlen: usize, num_reads: usize, alpha: f64, beta: f64) -> Vec<f64>{

    let interval_counts = interval_count;

    if interval_counts.iter().sum::<f32>() == 0.0 {
        return vec![0.0; tlen + 1];
    }


    // Compute the sum of probabilities
    let sum: f64 = interval_counts.iter().map(|&x| x as f64).sum();
    // Normalize the probabilities by dividing each element by the sum
    let p1: Vec<f64> = interval_counts.iter().map(|&prob| prob as f64 / sum).collect();


    let uniform_prob = vec![1.0 / (p1.len() as f64); p1.len()];
    let number_reads = num_reads; //t.ranges.len();
    let coefficient1 = ((-alpha) * (number_reads - 1) as f64).exp();
    let coefficient2 = ((-beta) * (1.0 - entropy_function(&p1))).exp();


    let final_prob: Vec<f64> = uniform_prob.iter().zip(p1.iter()).map(|(&uni_prob, &prob1)| (coefficient1 * uni_prob) + ((1.0 - coefficient1) * coefficient2 * prob1)).collect();

    // Compute the sum of probabilities
    let sum_final_prob: f64 = final_prob.iter().sum();
    // Normalize the probabilities by dividing each element by the sum
    let normalized_final_prob: Vec<f64> = final_prob.iter().map(|prob| prob / sum_final_prob).collect();

    let cdf: Vec<f64> = normalized_final_prob.iter().scan(0.0, |acc, &prob| {
    *acc += prob;
    Some(*acc)
    }).collect();

    
    cdf
    
}


pub fn entropy_prob(txps: &mut Vec<TranscriptInfo>, rate: &str, bins: &u32, threads: usize, alpha: f64, beta: f64) -> (Vec<usize>, Vec<usize>) {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));
    let num_discarded_reads = Arc::new(RwLock::new(vec![0; txps.len()]));

    txps.par_iter_mut().enumerate().for_each(|(i, t)| {
        
        let mut temp_prob: Vec<f32>;
        let number_reads = t.ranges.len();
        let tlen = t.len.get() as usize;

        let bin_lengths: Vec<f32>;
        let mut num_discarded_read_temp: usize = 0;
        let bin_coverage: Vec<f64>;

        match rate {
            "binomial" => {
                let bin_counts: Vec<f32>;
                (bin_counts, bin_lengths, num_discarded_read_temp, bin_coverage) = bin_transcript_normalize_counts(t, bins); //binning the transcript length and obtain the counts and length vectors
                //==============================================================================================
                ////fill out the number of reads and discarded reads
                {
                    let mut w1 =num_reads.write().unwrap();
                    (*w1)[i] = t.ranges.len();
                }
                {
                    let mut w2 =num_discarded_reads.write().unwrap();
                    (*w2)[i] = num_discarded_read_temp;
                }
                //==============================================================================================
                let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let prob_dr: Vec<f64> = multinomial_binomial_entropy(bin_counts, bin_lengths, distinct_rate, tlen, number_reads, alpha, beta);
                temp_prob = prob_dr.iter().map(|&x| x as f32).collect();

            }
            "multinomial" => {
                let bin_counts: Vec<u32>;
                (bin_counts, bin_lengths, num_discarded_read_temp, bin_coverage) = bin_transcript_decision_rule(t, bins); //binning the transcript length and obtain the counts and length vectors
                //==============================================================================================
                ////fill out the number of reads and discarded reads
                {
                    let mut w1 =num_reads.write().unwrap();
                    (*w1)[i] = t.ranges.len();
                }
                {
                    let mut w2 =num_discarded_reads.write().unwrap();
                    (*w2)[i] = num_discarded_read_temp;
                }
                //==============================================================================================
                let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                let f32_bin_count = bin_counts.iter().map(|&x| x as f32).collect();
                let prob_dr: Vec<f64> = multinomial_binomial_entropy(f32_bin_count, bin_lengths, distinct_rate, tlen, number_reads, alpha, beta);
                temp_prob = prob_dr.iter().map(|&x| x as f32).collect();

            }
            "1bin" => {
                let bin_counts: Vec<u32>;
                (bin_counts, bin_lengths, num_discarded_read_temp) = all_points_counts(t);
                //==============================================================================================
                ////fill out the number of reads and discarded reads
                {
                    let mut w1 =num_reads.write().unwrap();
                    (*w1)[i] = t.ranges.len();
                }
                {
                    let mut w2 =num_discarded_reads.write().unwrap();
                    (*w2)[i] = num_discarded_read_temp;
                }
                //==============================================================================================
                let f32_bin_count = bin_counts.iter().map(|&x| x as f32).collect();
                let prob_dr: Vec<f64> = one_bin_entropy(f32_bin_count, tlen, number_reads, alpha, beta);
                temp_prob = prob_dr.iter().map(|&x| x as f32).collect();

            }
            _ => {
                panic!("{:?} rate is not defined in the program", rate);
            }
        }

        t.coverage_prob = temp_prob;
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