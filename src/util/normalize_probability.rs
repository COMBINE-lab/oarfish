use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use rayon::prelude::*;
use tracing::{error, info, instrument};

/// The coverage probability of an alignment to `txp` spanning `[start, end)`: the coverage
/// probability of its evidence bins, averaged by how much of each it covers.
#[inline]
fn expected_cov_prob(txp: &TranscriptInfo, start: f64, end: f64) -> f64 {
    let tlen: f64 = txp.len.get() as f64;
    let coverage_probability: &[f64] = &txp.coverage_prob;
    let (mut total_weight, mut cov_prob) = (0.0f64, 0.0f64);
    crate::util::oarfish_types::for_each_evidence_bin(
        tlen,
        coverage_probability.len(),
        start,
        end,
        |b, frac| {
            total_weight += frac;
            cov_prob += frac * coverage_probability[b];
        },
    );
    if total_weight <= 0.0 {
        // A degenerate span: no evidence either way.
        total_weight = 1.0;
        cov_prob = 0.5;
    }

    if cov_prob.is_nan() || cov_prob.is_infinite() {
        error!("cov_prob: {:?}", cov_prob);
        error!("tlen: {:?}, bins: {}", tlen, coverage_probability.len());
        error!("start_aln: {}, end_aln: {}", start, end);
        panic!("Error: Invalid result. normalize_read_probs function.");
    }
    cov_prob / total_weight
}

/// Sets each alignment's coverage probability ([`expected_cov_prob`]), normalized over its
/// read's alignments. Reads are independent, and are processed in parallel.
#[instrument(skip(store, txp_info))]
pub fn normalize_read_probs(store: &mut InMemoryAlignmentStore, txp_info: &[TranscriptInfo]) {
    info!("normalizing read probabilities");
    let n = store.total_len();
    let mut normalized_coverage_prob: Vec<f64> = Vec::with_capacity(n);
    let out = normalized_coverage_prob.as_mut_ptr() as usize;
    let st: &InMemoryAlignmentStore = store;
    (0..st.len()).into_par_iter().for_each(|r| {
        let range = st.read_range(r);
        let (alns, _, _) = st.read(r);
        // SAFETY: slots `range` belong to read r alone, and the reads' ranges cover 0..n.
        let dst = unsafe {
            std::slice::from_raw_parts_mut((out as *mut f64).add(range.start), range.len())
        };
        let mut nprob_sum = 0.0f64;
        for (d, a) in dst.iter_mut().zip(alns) {
            *d = expected_cov_prob(&txp_info[a.ref_id as usize], a.start as f64, a.end as f64);
            nprob_sum += *d;
        }
        let sum: f64 = if nprob_sum > 0.0 { nprob_sum } else { 1.0 };
        dst.iter_mut().for_each(|d| *d /= sum);
    });
    // SAFETY: every slot was written above.
    unsafe { normalized_coverage_prob.set_len(n) };
    store.coverage_probabilities = normalized_coverage_prob;
    info!("done");
}
