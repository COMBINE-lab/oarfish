use crate::util::oarfish_types::TranscriptInfo;
use itertools::izip;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use tracing::{error, warn};

pub fn binomial_probability(
    interval_count: &[f32],
    interval_length: &[f32],
    distinct_rate: f64,
) -> Vec<f64> {
    let interval_counts = interval_count;
    let interval_lengths = interval_length;
    let count_sum = interval_counts.iter().sum::<f32>();
    const ZERO_THRESH: f64 = 1e-20;
    // @zzare: Where does the magic number 709 come from?
    const MAX_SCALE_COUNT: f64 = 709_f64;

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
                if count < 0.0 || length < 0.0 || distinct_rate < 0.0 {
                    warn!(
                        "count: {:?}, length: {:?}, rate:{:?}",
                        count, length, distinct_rate
                    );
                }
                (count as f64) / (length as f64 * distinct_rate)
            }
        })
        .collect();

    // compute the quantities (in the numerator and denominator) that we will
    // use to compute the binomial probabilities.
    //###########################################################################################
    // obtain the maximum value over the bins.
    // NOTE: f32::max(f32::NAN, x) = x for all x != f32::NAN.
    let max_count = interval_counts.iter().cloned().fold(f32::NAN, f32::max);
    assert!(
        !max_count.is_nan(),
        "Max bin count of NAN was encountered. Please report this issue on GitHub"
    );

    //let complementary_count: Vec<f32> = interval_counts.iter().map(|&count_val| count_sum - count_val).collect();
    //let max_comp_count = complementary_count.iter().cloned().fold(f32::NAN, f32::max);
    //let max_val = max_count.max(max_comp_count as f32);

    let max_val = max_count;
    let interval_count_modified: Vec<f32> = interval_counts
        .iter()
        .map(|&count_val| {
            if count_val == max_val {
                MAX_SCALE_COUNT as f32
            } else {
                (((count_val as f64) * MAX_SCALE_COUNT) / (max_val as f64)) as f32
            }
        })
        .collect();
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
            error!("num2 is: {:?}", num2);
            error!("prob and sum_vec and count is: {:?}\t {:?}\t {:?}", prob, sum_vec, count);
            panic!("Incorrect result. multinomial_probability function provides nan or infinite values for log_numerator3");
        }

        let num3 = if (1.0 - prob) > ZERO_THRESH {(1.0 - prob).ln() * (sum_vec - count) as f64} else { ZERO_THRESH.ln() * (sum_vec - count) as f64};
        if num3.is_nan() || num3.is_infinite() {
            error!("num3 is: {:?}", num3);
            error!("prob and sum_vec and count is: {:?}\t {:?}\t {:?}", prob, sum_vec, count);
            panic!("Incorrect result. multinomial_probability function provides nan or infinite values for log_numerator3");
        }

        (num2, num3)
    }).unzip();

    let result: Vec<f64> = izip!(log_denominator.clone(), log_numerator2.clone(), log_numerator3.clone()).map(
        |(denom, num2, num3)| {
            let res = (log_numerator1 - denom + num2 +num3).exp();
            if res.is_nan() || res.is_infinite(){ // || (res == 0.0 && *count != 0.0) {
                let len = probabilities.len();
                error!("{log_numerator1}, {denom}, {num2}, {num3}, {res}, {sum_vec}, {len}");
                error!("interval_counts = {:?}", interval_count_modified);
                error!("probabilities = {:?}", probabilities);
                let t: Vec<f64> = interval_counts.iter().map(|b| sum_vec as f64 - *b as f64).collect();
                error!("sum_vec - count = {:?}", t);
                panic!("Incorrect result. multinomial_probability function provides nan or infinite values for result");
            }
            res
        }).collect();

    //==============================================================================================
    //// Normalize the probabilities by dividing each element by the sum
    ////let normalized_prob: Vec<f64> = result.iter().map(|&prob| prob / sum).collect();
    let sum = result.iter().sum::<f64>();
    let normalized_prob: Vec<f64> = result
        .iter()
        .map(|&prob| {
            let normalized = prob / sum;
            if normalized.is_nan() {
                error!(
                    "Warning: Division resulted in NaN. prob: {}, sum: {}",
                    prob, sum
                );
                error!("interval_counts = {:?}", interval_count_modified);
                error!("Warning: result: {:?}", result);
                panic!("prob_function, normalized_prob is not valid!");
            }
            normalized
        })
        .collect();

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

pub fn binomial_continuous_prob(txps: &mut [TranscriptInfo], bins: &u32, threads: usize) {
    use tracing::info;
    use tracing::info_span;

    let _log_span = info_span!("binomial_continuous_prob").entered();
    info!("computing coverage probabilities");

    let compute_txp_coverage_probs = |_i: usize, t: &mut TranscriptInfo| {
        let temp_prob: Vec<f64> = if *bins != 0 {
            let min_cov = t.total_weight / 100.;
            t.coverage_bins.iter_mut().for_each(|elem| *elem += min_cov);
            let (bin_counts, bin_lengths) = t.get_normalized_counts_and_lengths(false);

            let distinct_rate: f64 = bin_counts
                .iter()
                .zip(bin_lengths.iter())
                .map(|(&count, &length)| (count as f64) / (length as f64))
                .sum();
            binomial_probability(&bin_counts, &bin_lengths, distinct_rate)
        } else {
            std::unimplemented!("coverage model with 0 bins is not currently implemented");
        };
        t.coverage_prob = temp_prob;
    };

    // if we are requesting only a single thread, then don't bother with
    // the overhead of e.g. creating a thread pool and doing parallel
    // iteration, etc.
    if threads > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            txps.par_iter_mut()
                .enumerate()
                .for_each(|(_i, t)| compute_txp_coverage_probs(_i, t));
        });
    } else {
        txps.iter_mut()
            .enumerate()
            .for_each(|(_i, t)| compute_txp_coverage_probs(_i, t));
    }
    info!("done");
}
