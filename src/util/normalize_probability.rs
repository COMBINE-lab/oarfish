use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use rayon::prelude::*;
use tracing::{error, info, instrument};

#[inline]
fn check_and_divide(tlen: f64, nbins: usize, start: f64, end: f64, tw: f64, cp: f64) -> f64 {
    let (mut total_weight, mut cov_prob) = (tw, cp);
    if total_weight <= 0.0 {
        // A degenerate span: no evidence either way.
        total_weight = 1.0;
        cov_prob = 0.5;
    }
    if cov_prob.is_nan() || cov_prob.is_infinite() {
        error!("cov_prob: {:?}", cov_prob);
        error!("tlen: {:?}, bins: {}", tlen, nbins);
        error!("start_aln: {}, end_aln: {}", start, end);
        panic!("Error: Invalid result. normalize_read_probs function.");
    }
    cov_prob / total_weight
}

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
    check_and_divide(
        tlen,
        coverage_probability.len(),
        start,
        end,
        total_weight,
        cov_prob,
    )
}

/// The bins an alignment overlaps, split as `for_each_bin_overlap` computes them: the first
/// bin, a run of interior bins `first + 1 .. last - 1` that it covers fully (each has
/// overlap fraction exactly 1), and the last bin `last - 1`. `None` when the span covers
/// fewer than three bins, or when rounding leaves an interior bin not fully covered.
#[inline]
fn split_span(len: f64, nbins: usize, start: f64, end: f64) -> Option<(usize, usize, f64, f64)> {
    if nbins < 3 || len <= 0.0 {
        return None;
    }
    let w = len / nbins as f64;
    let (s, e) = (start.min(end).max(0.0), start.max(end).min(len));
    if e <= s {
        return None;
    }
    let first = ((s / w).floor() as usize).min(nbins - 1);
    let last = ((e / w).ceil() as usize).clamp(first + 1, nbins);
    if last - first < 3 {
        return None;
    }
    // Bin starts `b * w` and ends `min((b + 1) * w, len)` are monotone in `b`, so the
    // interior bins are all fully covered when the first and the last of them are. The negated
    // comparisons also send a NaN to the general path.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(s <= (first + 1) as f64 * w) || !(e >= ((last - 1) as f64 * w).min(len)) {
        return None;
    }
    // first bin
    let (bs, be) = (first as f64 * w, ((first + 1) as f64 * w).min(len));
    let ov = e.min(be) - s.max(bs);
    let f0 = if ov > 0.0 { ov / (be - bs) } else { 0.0 };
    // last bin: the one the evidence rule may leave out (interior bins are full)
    let b = last - 1;
    let (bs, be) = (b as f64 * w, ((b + 1) as f64 * w).min(len));
    let ov = e.min(be) - s.max(bs);
    const FULL: f64 = 1.0 - 1e-9;
    let fl = if ov > 0.0 { ov / (be - bs) } else { 0.0 };
    let fl = if fl >= FULL { fl } else { 0.0 };
    Some((first, last, f0, fl))
}

/// [`expected_cov_prob`], with the fully covered interior bins summed directly: their
/// overlap fraction is exactly 1, so this skips the per-bin overlap arithmetic and division
/// while adding the same terms in the same order (the result is bit-identical). Spans of
/// fewer than three bins, and spans where rounding leaves an interior bin short, take
/// [`expected_cov_prob`].
#[inline]
fn expected_cov_prob_fast(txp: &TranscriptInfo, start: f64, end: f64) -> f64 {
    let tlen: f64 = txp.len.get() as f64;
    let p: &[f64] = &txp.coverage_prob;
    match split_span(tlen, p.len(), start, end) {
        None => expected_cov_prob(txp, start, end),
        Some((first, last, f0, fl)) => {
            let mut tw = 0.0f64;
            let mut cp = 0.0f64;
            if f0 > 0.0 {
                tw += f0;
                cp += f0 * p[first];
            }
            for &x in &p[first + 1..last - 1] {
                tw += 1.0;
                cp += x;
            }
            if fl > 0.0 {
                tw += fl;
                cp += fl * p[last - 1];
            }
            check_and_divide(tlen, p.len(), start, end, tw, cp)
        }
    }
}

/// Sets each alignment's coverage probability ([`expected_cov_prob`], via
/// [`expected_cov_prob_fast`]), normalized over its
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
            let t = a.ref_id as usize;
            let (s, e) = (a.start as f64, a.end as f64);
            *d = expected_cov_prob_fast(&txp_info[t], s, e);
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
