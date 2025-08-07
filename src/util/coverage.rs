use crate::util::logistic_probability::logistic_prob;
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::EMInfo;
use itertools::*;

use super::constants::EM_DENOM_THRESH;

pub fn update_coverage(
    emi: &mut EMInfo,
    counts: &[f64],
    header: &noodles_sam::header::Header,
    growth_rate: f64,
    bin_width: u32,
    threads: usize,
) -> anyhow::Result<()> {
    // set coverage based on current abundance estimates
    let mut txps = Vec::<usize>::new();
    let mut txp_probs = Vec::<f64>::new();
    for (alns, probs, _coverage_probs) in emi.eq_map.iter() {
        let mut denom = 0.0_f64;

        for (a, p) in izip!(alns, probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = 1.0;
            denom += counts[target_id] * prob * cov_prob;
        }

        txps.clear();
        txp_probs.clear();

        let mut denom2 = 0.0_f64;

        for (a, p) in izip!(alns, probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = 1.0;
            let nprob = ((counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
            txps.push(target_id);
            txp_probs.push(nprob);
            denom2 += nprob;
        }

        for p in txp_probs.iter_mut() {
            *p /= denom2;
        }

        for (a, p) in izip!(alns, txp_probs.iter()) {
            let tid = a.ref_id as usize;
            let prob = *p;
            if *p >= EM_DENOM_THRESH {
                emi.txp_info[tid].add_interval(a.start, a.end, prob);
            }
        }
    }

    //obtaining the Cumulative Distribution Function (CDF) for each transcript
    logistic_prob(header, emi.txp_info, growth_rate, &bin_width, threads);
    //Normalize the probabilities for the records of each read
    normalize_read_probs(emi.eq_map, emi.txp_info, &bin_width);
    Ok(())
}
