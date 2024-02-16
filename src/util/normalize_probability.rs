use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use tracing::{info, instrument};

#[instrument(skip(store, txp_info))]
pub fn normalize_read_probs(
    store: &mut InMemoryAlignmentStore,
    txp_info: &[TranscriptInfo],
    num_bins: &u32,
) {
    let mut normalize_probs_temp: Vec<f64> = vec![];
    let mut normalized_coverage_prob: Vec<f64> = vec![];

    info!("normalizing read probabilities");
    //iterate over all alignments in the bam file
    for (alns, _as_probs, _coverage_prob) in store.iter() {
        //iterate over the alignments of a read
        for a in alns.iter() {
            let target_id: usize = a.ref_id as usize;
            let start_aln: f32 = a.start as f32;
            let end_aln: f32 = a.end as f32;
            let tlen: f32 = txp_info[target_id].len.get() as f32;
            let bin_length: f32 = tlen / *num_bins as f32;
            let start_bin: usize = (start_aln / bin_length) as usize;
            let end_bin: usize = (end_aln / bin_length) as usize;
            let coverage_probability: &Vec<f64> = &txp_info[target_id].coverage_prob;

            let cov_prob: f64 = (start_bin..end_bin)
                .map(|i| {
                    match (
                        i == start_bin,
                        i == end_bin,
                        end_bin == (*num_bins - 1) as usize,
                    ) {
                        (true, false, false) => {
                            let coverage_prob = coverage_probability[i]
                                / (((i + 1) as f32 * bin_length).floor()
                                    - (i as f32 * bin_length).floor())
                                    as f64;
                            ((start_aln as usize)..((i + 1) * bin_length.ceil() as usize))
                                .map(|_j| coverage_prob)
                                .sum()
                        }
                        (false, true, true) => {
                            let coverage_prob = coverage_probability[i]
                                / (tlen + 1.0 - (i as f32 * bin_length).ceil()) as f64;
                            ((i * bin_length.ceil() as usize)..=(end_aln as usize))
                                .map(|_j| coverage_prob)
                                .sum()
                        }
                        (false, true, false) => {
                            let coverage_prob = coverage_probability[i]
                                / (((i + 1) as f32 * bin_length).floor()
                                    - (i as f32 * bin_length).floor())
                                    as f64;
                            ((i * bin_length.ceil() as usize)..=(end_aln as usize))
                                .map(|_j| coverage_prob)
                                .sum()
                        }
                        _ => coverage_probability[i],
                    }
                })
                .sum();

            if cov_prob.is_nan() || cov_prob.is_infinite() {
                panic!("Error: Invalid result. normalize_read_probs function.");
            }
            normalize_probs_temp.push(cov_prob);
        }

        let nprob_sum = normalize_probs_temp.iter().sum::<f64>();
        let sum_normalize_probs_temp: f64 = if nprob_sum > 0.0 { nprob_sum } else { 1.0 };

        let normalized_prob_section: Vec<f64> = normalize_probs_temp
            .iter()
            .map(|&prob| prob / sum_normalize_probs_temp)
            .collect();

        normalized_coverage_prob.extend(normalized_prob_section);
        normalize_probs_temp.clear();
    }
    store.coverage_probabilities = normalized_coverage_prob;
    info!("done");
}
