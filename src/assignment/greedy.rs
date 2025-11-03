use super::connected_components::TranscriptConnectedComponentLabeling;
use crate::{
    assignment,
    util::oarfish_types::{AlnInfo, EMInfo},
};
use itertools::*;

/// solve the assignment using the greedy heuristic. Assign each read to it's highest posterior
/// probability alignment, unless the corresponding transcript is full
pub fn solve<'a, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a, F: Fn() -> I>(
    em_info: &EMInfo,
    counts: &[f64],
    make_iter: F,
) -> Vec<u32> {
    let fops = &em_info.eq_map.filter_opts;
    let mut assignments = Vec::<u32>::new();
    let current_assignments = vec![0; em_info.txp_info.len()];
    let model_coverage = fops.model_coverage;

    for (alns, probs, coverage_probs) in em_info.eq_map.iter() {
        let mut max_idx = alns[0].ref_id as usize;
        let mut max_prob = 0f64;
        for (idx, (a, p, cp)) in izip!(alns, probs, coverage_probs).enumerate() {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;

            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            let cond_prob = counts[target_id] * prob * cov_prob;
            // if this is the best probability yet, and we can fit this read, then assign it.
            if cond_prob > max_prob
                && current_assignments[target_id] < (counts[target_id] as usize) + 1
            {
                max_prob = cond_prob;
                max_idx = idx;
            }
        }
        assignments.push(alns[max_idx].ref_id);
    }

    assignments
}
