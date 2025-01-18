use crate::util::oarfish_types::TranscriptInfo;
use itertools::izip;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use tracing::{error, warn};

/// implements a scaled (by `a`) logistic function
/// that is clamped (>= 1e-8, <= 0.99999)
fn logistic(x: f64, a: f64) -> f64 {
    let result = 1.0 / (1.0 + (-a * x).exp());
    result.clamp(1e-8, 0.99999)
}

pub fn logstic_function(interval_count: &[f32], interval_length: &[f32]) -> Vec<f64> {
    let interval_counts = interval_count;
    let interval_lengths = interval_length;
    let count_sum = interval_counts.iter().map(|&x| x as f64).sum::<f64>();

    if count_sum == 0.0 {
        return vec![0.0; interval_counts.len()];
    }

    // compute the expected value of the counts
    let expected_count = count_sum / interval_counts.len() as f64;

    let difference: Vec<f64> = interval_counts
        .iter()
        .map(|&count| (expected_count - (count as f64)) / expected_count as f64)
        .collect();

    let logistic_prob: Vec<f64> = difference
        .iter()
        .map(|&diff| logistic(diff, 4.0) as f64)
        .collect();

    //==============================================================================================
    //// Normalize the probabilities by dividing each element by the sum
    ////let normalized_prob: Vec<f64> = result.iter().map(|&prob| prob / sum).collect();
    /*
    let sum = logistic_prob.iter().sum::<f64>();
    let normalized_prob: Vec<f64> = logistic_prob
        .iter()
        .map(|&prob| {
            let normalized = prob / sum;
            if normalized.is_nan() {
                error!(
                    "Warning: Division resulted in NaN. prob: {}, sum: {}",
                    prob, sum
                );
                //error!("interval_counts = {:?}", interval_count_modified);
                //error!("Warning: result: {:?}", result);
                //panic!("prob_function, normalized_prob is not valid!");
            }
            normalized as f64
        })
        .collect();

    normalized_prob
    */
    logistic_prob
}

pub fn logistic_prob(txps: &mut [TranscriptInfo], bins: &u32, threads: usize) {
    use tracing::info;
    use tracing::info_span;

    let _log_span = info_span!("logistic_prob").entered();
    info!("computing coverage probabilities");

    let compute_txp_coverage_probs = |_i: usize, t: &mut TranscriptInfo| {
        let temp_prob: Vec<f64> = if *bins != 0 {
            let min_cov = t.total_weight / 100.;
            t.coverage_bins.iter_mut().for_each(|elem| *elem += min_cov);
            let (bin_counts, bin_lengths) = t.get_normalized_counts_and_lengths();

            logstic_function(&bin_counts, &bin_lengths)
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
