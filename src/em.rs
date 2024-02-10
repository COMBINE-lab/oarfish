use itertools::izip;
use tracing::info;
use num_format::{Locale, ToFormattedString};
use crate::util::oarfish_types::{
    EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};

/// Performs one iteration of the EM algorithm by looping over all
/// alignments and computing their estimated probability of being
/// the true alignment (using the abunance estimates from `prev_counts`).
/// Then, `curr_counts` is computed by summing over the expected assignment
/// likelihood for all reads mapping to each target.
#[inline]
fn m_step(
    eq_map: &InMemoryAlignmentStore,
    _tinfo: &mut [TranscriptInfo],
    model_coverage: bool,
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) {
    const DENOM_THRESH: f64 = 1e-30_f64;
    for (alns, probs, coverage_probs) in eq_map.iter() {
        let mut denom = 0.0_f64;
        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };

            denom += prev_count[target_id] * prob * cov_prob;
        }

        // If this read can be assigned
        if denom > DENOM_THRESH {
            // Loop over all possible assignment locations and proportionally
            // allocate the read according to our model and current parameter
            // estimates.
            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                let inc = (prev_count[target_id] * prob * cov_prob) / denom;
                curr_counts[target_id] += inc;
            }
        }
    }
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
pub fn em(em_info: &mut EMInfo) -> Vec<f64> {
    let eq_map = em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &mut Vec<TranscriptInfo> = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let convergence_thresh = em_info.convergence_thresh;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    // initialize the estimated counts for the EM procedure
    let mut prev_counts: Vec<f64>;

    if let Some(ref init_counts) = em_info.init_abundances {
        // initalize with the short-read quantification
        prev_counts = init_counts.clone();
    } else {
        // uniform, length normalized abundance
        let avg = total_weight / (tinfo.len() as f64);
        prev_counts = vec![avg; tinfo.len()];
    }

    let mut curr_counts = vec![0.0f64; tinfo.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;
    let mut _fl_prob = 0.5f64;

    // for up to the maximum number of iterations
    while niter < max_iter {
        // allocate the fragments and compute the new counts
        m_step(
            eq_map,
            tinfo,
            fops.model_coverage,
            &mut prev_counts,
            &mut curr_counts,
        );

        // compute the relative difference in the parameter estimates
        // between the current and previous rounds
        for i in 0..curr_counts.len() {
            if prev_counts[i] > 1e-8 {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = rel_diff.max(rd);
            }
        }

        // swap the current and previous abundances
        std::mem::swap(&mut prev_counts, &mut curr_counts);
        // clear out the new abundances
        curr_counts.fill(0.0_f64);

        // if the maximum relative difference is small enough
        // and we've done at least 10 rounds of the EM, then
        // exit (early stop).
        if (rel_diff < convergence_thresh) && (niter > 50) {
            break;
        }
        // increment the iteration and, if this iteration
        // is a multiple of 10, print out  the maximum relative
        // difference we observed.
        niter += 1;
        if niter % 10 == 0 {
            info!(
                "iteration {}; rel diff {}",
                niter.to_formatted_string(&Locale::en),
                rel_diff
            );
        }
        rel_diff = 0.0_f64;
    }

    // set very small abundances to 0
    prev_counts.iter_mut().for_each(|x| {
        if *x < 1e-8 {
            *x = 0.0
        }
    });
    // perform one more EM round, since we just zeroed out
    // very small abundances
    m_step(
        eq_map,
        tinfo,
        fops.model_coverage,
        &mut prev_counts,
        &mut curr_counts,
    );
    //  return the final estimated abundances
    curr_counts
}


