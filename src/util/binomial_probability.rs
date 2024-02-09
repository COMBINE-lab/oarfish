use crate::util::count_function::bin_transcript_normalize_counts;
use crate::util::oarfish_types::TranscriptInfo;
use itertools::izip;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;

pub fn binomial_probability(
    interval_count: &[f32],
    interval_length: &[f32],
    distinct_rate: f64,
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
                //eprintln!("count: {:?}, length: {:?}, rate:{:?}", count, length, distinct_rate);
                (count as f64) / (length as f64 * distinct_rate)
            }
        })
        .collect();

    // compute the quantities (in the numerator and denominator) that we will
    // use to compute the binomial probabilities.
    let sum_vec = count_sum;
    let log_numerator1: f64 = ln_gamma(sum_vec as f64 + 1.0);
    let log_denominator: Vec<f64> = interval_counts
        .iter()
        .map(|&count| ln_gamma(count as f64 + 1.0) + ln_gamma((sum_vec - count) as f64 + 1.0))
        .collect();

    let (log_numerator2, log_numerator3) : (Vec<f64>, Vec<f64>) = probabilities.iter().zip(interval_counts.iter()).map(|(&prob, &count)| {
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

    let result: Vec<f64> = izip!(log_denominator, log_numerator2, log_numerator3).map(
        |(denom, num2, num3)| {
            let res = (log_numerator1 - denom + num2 +num3).exp();
            if res.is_nan() || res.is_infinite() {
                eprintln!("{log_numerator1}, {denom}, {num2}, {num3}");
                eprintln!("interval_counts = {:?}", &interval_counts);
                panic!("Incorrect result. multinomial_probability function provides nan or infinite values for result");
            }
            res
        }).collect();

    // Compute the sum of probabilities
    let sum: f64 = result.iter().sum();

    // Normalize the probabilities by dividing each element by the sum
    let normalized_prob: Vec<f64> = result.iter().map(|&prob| prob / sum).collect();

    normalized_prob
}

pub fn binomial_continuous_prob(txps: &mut Vec<TranscriptInfo>, bins: &u32, threads: usize) {
    use tracing::info;
    use tracing::info_span;

    let _log_span = info_span!("binomial_continuous_prob").entered();
    info!("computing coverage probabilities");

    rayon::ThreadPoolBuilder::new()
        .num_threads(1)//threads)
        .build()
        .unwrap();

    txps.par_iter_mut().enumerate().for_each(|(_i, t)| {
        let temp_prob: Vec<f64>;

        if *bins != 0 {
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
            let (bin_counts, bin_lengths) = t.get_normalized_counts_and_lengths();

            let distinct_rate: f64 = bin_counts
                .iter()
                .zip(bin_lengths.iter())
                .map(|(&count, &length)| (count as f64) / (length as f64))
                .sum();
            temp_prob = binomial_probability(&bin_counts, &bin_lengths, distinct_rate);
        } else {
            std::unimplemented!("coverage model with 0 bins is not currently implemented");
        }

        t.coverage_prob = temp_prob;
    });
}
