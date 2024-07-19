use crate::util::oarfish_types::TranscriptInfo;
use itertools::izip;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use statrs::distribution::{Normal, ContinuousCDF};

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

pub fn sigmoid_probability(
    interval_count: &[f32],
    interval_length: &[f32],
    distinct_rate: f64,
    tname: &String,
    tlen: f64,
    bin_len: usize,
) -> Vec<f64> {
    let interval_counts = interval_count;
    let interval_lengths = interval_length;
    let count_sum = interval_counts.iter().sum::<f32>();

    if count_sum == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    if distinct_rate == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    let probabilities: Vec<f64> = interval_counts
        .iter()
        .zip(interval_lengths.iter())
        .map(|(&count, &length)| {
            if count == 0.0 || length == 0.0 {
                0.0
            } else {
                if count <0.0 || length < 0.0 || distinct_rate <0.0 {
                eprintln!("count: {:?}, length: {:?}, rate:{:?}", count, length, distinct_rate);
                }
                (count as f64) / (length as f64 * distinct_rate)
            }
        })
        .collect();

    //###########################################################################################
    // obtain the maximum value of 709 on each bin count
    let avg_cov_prob: f64 = probabilities.iter().sum::<f64>() / (tlen / bin_len as f64);
    let bin_diff_avg: Vec<f64> = probabilities.iter().map(|&prob| avg_cov_prob - prob).collect();
    let sigmoid_prob: Vec<f64> = bin_diff_avg.iter().map(|&x| sigmoid(x)).collect();
    //###########################################################################################
    let sum = sigmoid_prob.iter().sum::<f64>();
    let normalized_prob: Vec<f64> = sigmoid_prob.iter().map(|&prob| {
        let normalized = prob / sum;
        if normalized.is_nan() {
            eprintln!("Warning: Division resulted in NaN. prob: {}, sum: {}", prob, sum);
            eprintln!("interval_counts = {:?}", interval_counts);
            eprintln!("Warning: result: {:?}", sigmoid_prob);
            panic!("prob_function, normalized_prob is not valid!");
        }
        normalized
    }).collect();

    normalized_prob
}




pub fn normal_probability(
    interval_count: &[f32],
    interval_length: &[f32],
    distinct_rate: f64,
    tname: &String,
) -> Vec<f64> {
    let interval_counts = interval_count;
    let interval_lengths = interval_length;
    let count_sum = interval_counts.iter().sum::<f32>();
    const ZERO_THRESH: f64 = 1e-20;

    if count_sum == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    if distinct_rate == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    let probabilities: Vec<f64> = interval_counts
        .iter()
        .zip(interval_lengths.iter())
        .map(|(&count, &length)| {
            if count == 0.0 || length == 0.0 {
                0.0
            } else {
                if count <0.0 || length < 0.0 || distinct_rate <0.0 {
                eprintln!("count: {:?}, length: {:?}, rate:{:?}", count, length, distinct_rate);
                }
                (count as f64) / (length as f64 * distinct_rate)
            }
        })
        .collect();

    // compute the quantities (in the numerator and denominator) that we will
    // use to compute the binomial probabilities.
    let sum_vec = count_sum as f64;

    let normal_prob : Vec<f64> = probabilities.iter().zip(interval_counts.iter()).map(|(&prob, &count)| {
        let mean = if prob > ZERO_THRESH { sum_vec * prob } else { sum_vec * ZERO_THRESH };
        if mean.is_nan() || mean.is_infinite() {
            eprintln!("mean is: {:?}", mean);
            eprintln!("prob and sum_vec and count is: {:?}\t {:?}\t {:?}", prob, sum_vec, count);
            panic!("Incorrect result. normal_probability function provides nan or infinite values for mean");
        }

        let std = if ((prob <= ZERO_THRESH) || ((1.0 - prob) <= ZERO_THRESH)) {f64::sqrt(sum_vec * ZERO_THRESH * (1.0 - ZERO_THRESH)) as f64} else { f64::sqrt(sum_vec * prob * (1.0 - prob)) as f64};
        if std.is_nan() || std.is_infinite() {
            eprintln!("std is: {:?}", std);
            eprintln!("prob and sum_vec and count is: {:?}\t {:?}\t {:?}", prob, sum_vec, count);
            panic!("Incorrect result. normal_probability function provides nan or infinite values for std");
        }

        let n_p = Normal::new(mean, std).unwrap();

        let result = (n_p.cdf(count as f64 + 0.5) - n_p.cdf(count as f64 - 0.5));

        eprintln!("prob: {}, count: {}, mean: {}, std: {}, result:{}", prob, count, mean, std, result);

        result
    }).collect();



    // Compute the sum
    let sum: f64 = normal_prob
        .iter().sum();

    // Normalize the probabilities
    let normalized_prob: Vec<f64> = normal_prob
        .iter()
        .map(|&normal_p| {
            let normalized = normal_p / sum;
            if normalized.is_nan() {
                eprintln!("Warning: Division resulted in NaN. normal_p: {}, sum: {}, normalized: {}, interval_counts: {:?}", normal_p, sum, normalized, interval_counts);
                panic!("prob_function, normalized_prob is not valid!");
            } else if normalized.is_infinite() {
                eprintln!("Warning: Division resulted in inf. normal_p: {}, sum: {}, normalized: {}", normal_p, sum, normalized);
                panic!("prob_function, normalized_prob is not valid!");
            }
            normalized
        })
        .collect();

    if normalized_prob.iter().any(|&prob| prob < 0.0 || prob > 1.0) {
        eprintln!("Warning: Normalized probability out of bounds. Normalized probabilities: {:?}", normalized_prob);
        panic!("prob_function, normalized_prob out of valid range!");
    }

    let normalized_sum = normalized_prob.iter().sum::<f64>();
    if(count_sum != 0.0 && normalized_sum == 0.0){
        panic!("warning in binomial_probability function: Numerical Instability. the count in each bin is non-zero but the resultant probability is zero.");
    }


    normalized_prob
}

pub fn binomial_probability(
    interval_count: &[f32],
    interval_length: &[f32],
    distinct_rate: f64,
    tname: &String,
) -> Vec<f64> {
    let interval_counts = interval_count;
    let interval_lengths = interval_length;
    let count_sum = interval_counts.iter().sum::<f32>();
    const ZERO_THRESH: f64 = 1e-20;

    if count_sum == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    if distinct_rate == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    let probabilities: Vec<f64> = interval_counts
        .iter()
        .zip(interval_lengths.iter())
        .map(|(&count, &length)| {
            if count == 0.0 || length == 0.0 {
                0.0
            } else {
                if count <0.0 || length < 0.0 || distinct_rate <0.0 {
                eprintln!("count: {:?}, length: {:?}, rate:{:?}", count, length, distinct_rate);
                }
                (count as f64) / (length as f64 * distinct_rate)
            }
        })
        .collect();

    // compute the quantities (in the numerator and denominator) that we will
    // use to compute the binomial probabilities.
    //###########################################################################################
    // obtain the maximum value of 709 on each bin count
    let max_count = interval_counts.iter().cloned().fold(f32::NAN, f32::max);
    //let complementary_count: Vec<f32> = interval_counts.iter().map(|&count_val| count_sum - count_val).collect();
    //let max_comp_count = complementary_count.iter().cloned().fold(f32::NAN, f32::max);
    //let max_val = max_count.max(max_comp_count as f32);
    let max_val = max_count;
    let interval_count_modified: Vec<f32> = interval_counts.iter().map(|&count_val| if count_val == max_val {709.0} else {count_val * 709.0 / max_val}).collect();
    let sum_vec = interval_count_modified.iter().sum::<f32>();
    //###########################################################################################
    //let sum_vec = count_sum;
    let log_numerator1: f64 = ln_gamma(sum_vec as f64 + 1.0);
    let log_denominator: Vec<f64> = interval_count_modified
        .iter()
        .map(|&count| ln_gamma(count as f64 + 1.0) + ln_gamma((sum_vec - count) as f64 + 1.0))
        .collect();

    let (log_numerator2, log_numerator3) : (Vec<f64>, Vec<f64>) = probabilities.iter().zip(interval_count_modified.iter()).map(|(&prob, &count)| {
        let num2 = if prob > ZERO_THRESH { prob.ln() * (count as f64) } else { ZERO_THRESH.ln() * (count as f64) };
        if num2.is_nan() || num2.is_infinite() {
            eprintln!("num2 is: {:?}", num2);
            eprintln!("prob and sum_vec and count is: {:?}\t {:?}\t {:?}", prob, sum_vec, count);
            panic!("Incorrect result. multinomial_probability function provides nan or infinite values for log_numerator3");
        }

        let num3 = if (1.0 - prob) > ZERO_THRESH {(1.0 - prob).ln() * (sum_vec - count) as f64} else { ZERO_THRESH.ln() * (sum_vec - count) as f64};
        if num3.is_nan() || num3.is_infinite() {
            eprintln!("num3 is: {:?}", num3);
            eprintln!("prob and sum_vec and count is: {:?}\t {:?}\t {:?}", prob, sum_vec, count);
            panic!("Incorrect result. multinomial_probability function provides nan or infinite values for log_numerator3");
        }

        (num2, num3)
    }).unzip();

    let result: Vec<f64> = izip!(log_denominator.clone(), log_numerator2.clone(), log_numerator3.clone()).map(
        |(denom, num2, num3)| {
            let res = (log_numerator1 - denom + num2 +num3).exp();
            if res.is_nan() || res.is_infinite(){ // || (res == 0.0 && *count != 0.0) {
                let len = probabilities.len();
                eprintln!("{log_numerator1}, {denom}, {num2}, {num3}, {res}, {sum_vec}, {len}");
                eprintln!("interval_counts = {:?}", interval_count_modified);
                eprintln!("probabilities = {:?}", probabilities);
                let t: Vec<f64> = interval_counts.iter().map(|b| sum_vec as f64 - *b as f64).collect();
                eprintln!("sum_vec - count = {:?}", t);
                panic!("Incorrect result. multinomial_probability function provides nan or infinite values for result");
            }
            res
        }).collect();


//==============================================================================================
    //// Normalize the probabilities by dividing each element by the sum
    ////let normalized_prob: Vec<f64> = result.iter().map(|&prob| prob / sum).collect();
    let sum = result.iter().sum::<f64>();
    let normalized_prob: Vec<f64> = result.iter().map(|&prob| {
        let normalized = prob / sum;
        if normalized.is_nan() {
            eprintln!("Warning: Division resulted in NaN. prob: {}, sum: {}", prob, sum);
            eprintln!("interval_counts = {:?}", interval_count_modified);
            eprintln!("Warning: result: {:?}", result);
            panic!("prob_function, normalized_prob is not valid!");
        }
        normalized
    }).collect();


    //new method of log values
    //let log_values: Vec<f64> = izip!(log_denominator, log_numerator2, log_numerator3)
    //.map(|(denom, num2, num3)| log_numerator1 - denom + num2 + num3)
    //.collect();
//
    //// Find the maximum log value to use for stability
    //let max_log_value = log_values
    //    .iter()
    //    .cloned()
    //    .fold(f64::NEG_INFINITY, f64::max);
//
    //// Compute the exponential sum in a numerically stable way
    //let exp_sum: f64 = log_values
    //    .iter()
    //    .map(|&log_val| (log_val - max_log_value).exp())
    //    .sum();
//
    //// Normalize the probabilities
    //let normalized_prob: Vec<f64> = log_values
    //    .iter()
    //    .map(|&log_val| {
    //        let normalized = (log_val - max_log_value).exp() / exp_sum;
    //        if normalized.is_nan() {
    //            eprintln!("Warning: Division resulted in NaN. log_val: {}, max_log_value: {}, exp_sum: {}", log_val, max_log_value, exp_sum);
    //            panic!("prob_function, normalized_prob is not valid!");
    //        } else if normalized.is_infinite() {
    //            eprintln!("Warning: Division resulted in inf. log_val: {}, max_log_value: {}, exp_sum: {}", log_val, max_log_value, exp_sum);
    //            panic!("prob_function, normalized_prob is not valid!");
    //        }
    //        normalized
    //    })
    //    .collect();
//
    //if normalized_prob.iter().any(|&prob| prob < 0.0 || prob > 1.0) {
    //    eprintln!("Warning: Normalized probability out of bounds. Normalized probabilities: {:?}", normalized_prob);
    //    panic!("prob_function, normalized_prob out of valid range!");
    //}
//
    //let normalized_sum = normalized_prob.iter().sum::<f64>();
    //if(count_sum != 0.0 && normalized_sum == 0.0){
    //    panic!("warning in binomial_probability function: Numerical Instability. the count in each bin is non-zero but the resultant probability is zero.");
    //}


    normalized_prob
}

pub fn binomial_continuous_prob(txps: &mut Vec<TranscriptInfo>, threads: usize, txps_name: &Vec<String>, bin_len: usize) {
    use tracing::info;
    use tracing::info_span;

    let _log_span = info_span!("binomial_continuous_prob").entered();
    info!("computing coverage probabilities");
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();

    txps.par_iter_mut().enumerate().for_each(|(i, t)| {
        let temp_prob: Vec<f64> = if t.coverage_bins.len() != 0 {
            /*
            let bin_counts: Vec<f32>;
            let bin_lengths: Vec<f32>;
            let _num_discarded_read_temp: usize;
            let _bin_coverage: Vec<f64>;
            (
                bin_counts,
                bin_lengths,
                _num_discarded_read_temp,
                _bin_coverage,
            ) = bin_transcript_normalize_counts(t, bins); //binning the transcript length and obtain the counts and length vectors
                                                          //==============================================================================================
            */
            //let constant = t.num_read / 10.0;
            //let constant = 1.0;
            let constant = t.num_read / 100.0;
            t.coverage_bins.iter_mut().for_each(|elem| *elem += constant);
            let (bin_counts, bin_lengths) = t.get_normalized_counts_and_lengths();

            let distinct_rate: f64 = bin_counts
                .iter()
                .zip(bin_lengths.iter())
                .map(|(&count, &length)| (count as f64) / (length as f64))
                .sum();
            let tname = &txps_name[i];
            let tlen = t.lenf;
            binomial_probability(&bin_counts, &bin_lengths, distinct_rate, tname)
            //sigmoid_probability(&bin_counts, &bin_lengths, distinct_rate, tname, tlen, bin_len)
            //normal_probability(&bin_counts, &bin_lengths, distinct_rate, tname)
        } else {
            std::unimplemented!("coverage model with 0 bins is not currently implemented");
        };

        t.coverage_prob = temp_prob;
    });
    info!("done");
}
