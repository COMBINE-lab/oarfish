use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use tracing::{error, info, instrument};

#[instrument(skip(store, txp_info))]
pub fn normalize_read_probs(store: &mut InMemoryAlignmentStore, txp_info: &[TranscriptInfo]) {
    let mut normalize_probs_temp: Vec<f64> = vec![];
    let mut normalized_coverage_prob: Vec<f64> = vec![];

    info!("normalizing read probabilities");
    //iterate over all alignments in the bam file
    for (alns, _as_probs, _coverage_prob) in store.iter() {
        let mut nprob_sum = 0.0f64;
        //iterate over the alignments of a read
        for a in alns.iter() {
            let target_id: usize = a.ref_id as usize;
            let start_aln: f64 = a.start as f64;
            let end_aln: f64 = a.end as f64;
            let tlen: f64 = txp_info[target_id].len.get() as f64;
            let coverage_probability: &Vec<f64> = &txp_info[target_id].coverage_prob;
            // The read's weight on this transcript is the coverage probability of its evidence
            // bins, averaged by how much of each it covers.
            let (mut total_weight, mut cov_prob) = (0.0f64, 0.0f64);
            crate::util::oarfish_types::for_each_evidence_bin(
                tlen,
                coverage_probability.len(),
                start_aln,
                end_aln,
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
                error!("start_aln: {}, end_aln: {}", start_aln, end_aln);
                panic!("Error: Invalid result. normalize_read_probs function.");
            }
            let expected_cov_prob = cov_prob / (total_weight);
            nprob_sum += expected_cov_prob;
            normalize_probs_temp.push(expected_cov_prob);
        }
        let sum_normalize_probs_temp: f64 = if nprob_sum > 0.0 { nprob_sum } else { 1.0 };

        // extend the total normalized vector with the current normalized sums.
        normalized_coverage_prob.extend(
            normalize_probs_temp
                .iter()
                .map(|&prob| prob / sum_normalize_probs_temp),
        );

        normalize_probs_temp.clear();
    }
    store.coverage_probabilities = normalized_coverage_prob;
    info!("done");
}
