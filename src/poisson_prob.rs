use crate::variables::TranscriptInfo;
use crate::common_functions::{bin_transcript_normalize_counts, bin_transcript_decision_rule, factorial_ln};
use special::Gamma;
use rayon::prelude::*;
use std::sync::{RwLock, Arc};


//https://en.wikipedia.org/wiki/Poisson_distribution#Computational_methods
pub fn continuous_poisson_probability(k: f32, t: f64, r: f64) -> f64{

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

pub fn poisson_probability(k: u32, t: f64, r: f64) -> f64{

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

pub fn discrete_poisson_prob(txps: &mut Vec<TranscriptInfo>, rate: &str, bins: &u32, threads: usize) -> (Vec<usize>, Vec<usize>) {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let mut num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));
    let mut num_discarded_reads = Arc::new(RwLock::new(vec![0; txps.len()]));

    txps.par_iter_mut().enumerate().for_each(|(i, t)| {
        
        let mut temp_prob: Vec<f64>;
        let mut start_elemtns: Vec<f32> = vec![];

        if *bins != 0 {
            let bin_counts: Vec<u32>;
            let bin_lengths: Vec<f32>;
            let mut num_discarded_read_temp: usize = 0;
            let bin_coverage: Vec<f64>;
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
            start_elemtns = bin_lengths.clone();
            //==============================================================================================
            match rate {
                "shr" | "shr_tlen" => {
                    let shared_rate: f64 = (bin_counts.iter().sum::<u32>() as f64) / (bin_lengths.iter().sum::<f32>() as f64);
                    let prob_shr: Vec<f64> = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, shared_rate)).collect();
                    temp_prob = prob_shr;

                }
                "dr" => {
                    let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                    let prob_dr: Vec<f64> = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, distinct_rate)).collect();
                    temp_prob = prob_dr;
                }
                _ => {
                    panic!("{:?} rate is not defined in the program", rate);
                }
            }
        }else {

            //fill out the number of reads and discarded reads
            {
                let mut w1 =num_reads.write().unwrap();
                (*w1)[i] = t.ranges.len();
            }

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
            //==============================================================================================
            start_elemtns = interval_length.iter().map(|&len| len as f32).collect::<Vec<_>>().clone();
            //==============================================================================================
            match rate {
                "shr" => {
                    let shared_rate: f64 = (interval_counts.iter().sum::<u32>() as f64) / (interval_length.iter().sum::<u32>() as f64);
                    let prob_shr = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, shared_rate)).collect();
                    temp_prob = prob_shr;

                }
                "shr_tlen" => {
                    let shared_rate_tlen: f64 = (interval_counts.iter().sum::<u32>() as f64) / (len as f64);
                    let prob_shr_tlen = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, shared_rate_tlen)).collect();
                    temp_prob = prob_shr_tlen;
                }
                "dr" => {
                    let distinct_rate: f64 =  interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                    let prob_dr= interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| poisson_probability(count, length as f64, distinct_rate)).collect();
                    temp_prob = prob_dr;
                }
                _ => {
                    panic!("{:?} rate is not defined in the program", rate);
                }
            }
        }

        let tlen = t.len.get();
        let mut prob_vec = vec![0.0; tlen + 1];
        let mut bin_start = 0;
        for (i, bin_length) in start_elemtns.iter().enumerate() {
            let bin_end = bin_start + (*bin_length).floor() as usize;
        
            let start_index = bin_start;
            let end_index = if i != (start_elemtns.len() - 1) {bin_end} else {tlen + 1};

            prob_vec[start_index as usize..end_index as usize]
                .iter_mut()
                .for_each(|v| *v = temp_prob[i as usize]);

            bin_start = bin_end;
        }
        let cdf: Vec<f64> = prob_vec.iter().scan(0.0, |acc, &prob| {
            *acc += prob;
            Some(*acc)
        }).collect();
        t.coverage_prob = cdf.iter().map(|&x| x as f32).collect();
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


pub fn continuous_poisson_prob(txps: &mut Vec<TranscriptInfo>, rate: &str, bins: &u32, threads: usize) -> (Vec<usize>, Vec<usize>) {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let mut num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));
    let mut num_discarded_reads = Arc::new(RwLock::new(vec![0; txps.len()]));

    txps.par_iter_mut().enumerate().for_each(|(i, t)| {
        
        let mut temp_prob: Vec<f64>;
        let mut start_elemtns: Vec<f32> = vec![];

        if *bins != 0 {
            let bin_counts: Vec<f32>;
            let bin_lengths: Vec<f32>;
            let mut num_discarded_read_temp: usize = 0;
            let bin_coverage: Vec<f64>;
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
            start_elemtns = bin_lengths.clone();
            //==============================================================================================
            match rate {
                "shr" | "shr_tlen" => {
                    let shared_rate: f64 = (bin_counts.iter().sum::<f32>() as f64) / (bin_lengths.iter().sum::<f32>() as f64);
                    //let prob_shr: f64 = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, shared_rate)).product();
                    let prob_shr: Vec<f64> = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, shared_rate)).collect();
                    temp_prob = prob_shr;

                }
                "dr" => {
                    let distinct_rate: f64 =  bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                    //let prob_dr: f64 = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, distinct_rate)).product();
                    let prob_dr: Vec<f64> = bin_counts.iter().zip(bin_lengths.iter()).map(|(&count, &length)| continuous_poisson_probability(count, length as f64, distinct_rate)).collect();
                    temp_prob = prob_dr;
                }
                _ => {
                    panic!("{:?} rate is not defined in the program", rate);
                }
            }
        }else {

            //fill out the number of reads and discarded reads
            {
                let mut w1 =num_reads.write().unwrap();
                (*w1)[i] = t.ranges.len();
            }

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
            //==============================================================================================
            start_elemtns = interval_length.iter().map(|&len| len as f32).collect::<Vec<_>>().clone();
            //==============================================================================================
            match rate {
                "shr" => {
                    let shared_rate: f64 = (interval_counts.iter().sum::<u32>() as f64) / (interval_length.iter().sum::<u32>() as f64);
                    let prob_shr = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| continuous_poisson_probability(count as f32, length as f64, shared_rate)).collect();
                    temp_prob = prob_shr;

                }
                "shr_tlen" => {
                    let shared_rate_tlen: f64 = (interval_counts.iter().sum::<u32>() as f64) / (len as f64);
                    let prob_shr_tlen = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| continuous_poisson_probability(count as f32, length as f64, shared_rate_tlen)).collect();
                    temp_prob = prob_shr_tlen;
                }
                "dr" => {
                    let distinct_rate: f64 =  interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| (count as f64)/ (length as f64)).sum();
                    let prob_dr = interval_counts.iter().zip(interval_length.iter()).map(|(&count, &length)| continuous_poisson_probability(count as f32, length as f64, distinct_rate)).collect();
                    temp_prob = prob_dr;
                }
                _ => {
                    panic!("{:?} rate is not defined in the program", rate);
                }
            }
        }

        let tlen = t.len.get();
        let mut prob_vec = vec![0.0; tlen + 1];
        let mut bin_start = 0;
        for (i, bin_length) in start_elemtns.iter().enumerate() {
            let bin_end = bin_start + (*bin_length).floor() as usize;
        
            let start_index = bin_start;
            let end_index = if i != start_elemtns.len() {bin_end} else {tlen + 1};

            prob_vec[start_index as usize..end_index as usize]
                .iter_mut()
                .for_each(|v| *v = temp_prob[i as usize]);

            bin_start = bin_end;
        }
        let cdf: Vec<f64> = prob_vec.iter().scan(0.0, |acc, &prob| {
            *acc += prob;
            Some(*acc)
        }).collect();
        t.coverage_prob = cdf.iter().map(|&x| x as f32).collect();
    
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

