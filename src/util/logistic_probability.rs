use crate::util::oarfish_types::TranscriptInfo;
use rayon::prelude::*;
use tracing::{info, instrument, warn};

/// implements a scaled (by `a`) logistic function
/// that is clamped (>= 1e-8, <= 0.99999)
fn logistic(x: f64, a: f64) -> f64 {
    let result = 1.0 / (1.0 + (-a * x).exp());
    result.clamp(1e-8, 0.99999)
}

pub fn logstic_function(
    growth_rate: f64,
    interval_count: &[f32],
    interval_length: &[f32],
    verbose: bool,
) -> Vec<f64> {
    let interval_counts = interval_count;
    let _interval_lengths = interval_length;
    let count_sum = interval_counts.iter().map(|&x| x as f64).sum::<f64>();

    if count_sum <= 1e-8 {
        return vec![0.0; interval_counts.len()];
    }

    // compute the expected value of the counts; just the sum of
    // bin counts divided by the number of bins.
    let expected_count = count_sum / interval_counts.len() as f64;
    if verbose {
        eprintln!("expected count: {expected_count}");
    }
    let logistic_prob: Vec<f64> = interval_counts
        .iter()
        //.enumerate()
        .map(|&count| {
            /*
            let mut n = 0.0;
            let mut ts = 0.0;
            for ix in i.saturating_sub(2)..(i + 2).min(interval_counts.len()) {
                let d = 0.5 * ((i as isize) - (ix as isize)).abs() as f64;
                let w = (d * d).exp() / (2. * PI).sqrt();
                ts += w * interval_counts[ix] as f64;
                n += w;
            }
            let expected_count = ts / (n as f64);
            */
            let diff = (expected_count - (count as f64)) / expected_count;
            if verbose {
                eprintln!("diff : {diff}");
            }
            logistic(diff, growth_rate)
        })
        .collect();
    logistic_prob
}

#[instrument(skip(txps, header))]
pub fn logistic_prob(
    header: &noodles_sam::header::Header,
    txps: &mut [TranscriptInfo],
    growth_rate: f64,
    bin_width: &u32,
    threads: usize,
) {
    info!("computing coverage probabilities");
    let compute_txp_coverage_probs = |_i: usize, t: &mut TranscriptInfo| {
        let verbose = header.reference_sequences().get_index(_i).expect("foo").0
            == bstr::B(b"SIRV303")
            || header.reference_sequences().get_index(_i).expect("foo").0 == bstr::B(b"SIRV301");

        let temp_prob: Vec<f64> = if *bin_width != 0 {
            assert!(!t.coverage_bins.is_empty());
            let min_cov = t.total_weight / 100.;
            t.coverage_bins.iter_mut().for_each(|elem| *elem += min_cov);
            let (bin_counts, bin_lengths) = t.get_normalized_counts_and_lengths(verbose);
            logstic_function(growth_rate, &bin_counts, &bin_lengths, verbose)
        } else {
            std::unimplemented!("coverage model with 0 bin width is not currently implemented");
        };
        t.coverage_prob = temp_prob;
        /*
        if header.reference_sequences().get_index(_i).expect("foo").0 == bstr::B(b"SIRV303") {
            eprintln!("SIRV303: {:#?}", t.coverage_prob);
        } else if header.reference_sequences().get_index(_i).expect("foo").0 == bstr::B(b"SIRV301")
        {
            eprintln!("SIRV301: {:#?}", t.coverage_prob);
        }
        */
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
