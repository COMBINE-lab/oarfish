use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use tracing::{error, info, instrument};

#[instrument(skip(store, txp_info))]
pub fn normalize_read_probs(
    store: &mut InMemoryAlignmentStore,
    txp_info: &[TranscriptInfo],
    bin_width: &u32,
) {
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
            let bin_length: f64 = *bin_width as f64; //txp_info[target_id].len.get() as f32;
            let start_bin: usize = (start_aln / bin_length) as usize;
            let end_bin: usize =
                ((end_aln / bin_length) as usize).min(coverage_probability.len() - 1);

            let bin_end = |i| -> f64 { ((bin_length * (i as f64)) + bin_length).min(tlen) };

            let bin_start = |i| -> f64 { (bin_length * (i as f64)).min(tlen) };

            let (total_weight, cov_prob): (f64, f64) = if start_bin == end_bin {
                let w = (end_aln - start_aln) / bin_length;
                (w, w * coverage_probability[start_bin])
            } else {
                (start_bin..end_bin).fold((0f64, 0f64), |acc, i| {
                    let w = if i == start_bin {
                        (bin_end(start_bin) - start_aln) / bin_length
                    } else if i == end_bin {
                        (end_aln - bin_start(end_bin)) / bin_length
                    } else {
                        1.0
                    };
                    (acc.0 + w, acc.1 + (w * coverage_probability[i]))
                })
            };

            if cov_prob.is_nan() || cov_prob.is_infinite() {
                error!("cov_prob: {:?}", cov_prob);
                error!("length: {:?}", ((end_aln - start_aln) / bin_length));
                error!("length2: {:?}", (tlen / bin_length));
                error!("start_bin: {}, end_bin: {}", start_bin, end_bin);
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
