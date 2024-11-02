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
    for (alns, _as_probs, _coverage_prob, _read_names) in store.iter() {
        //iterate over the alignments of a read
        for a in alns.iter() {
            let target_id: usize = a.ref_id as usize;
            let start_aln: f32 = a.start as f32;
            let end_aln: f32 = a.end as f32;
            let tlen: f32 = txp_info[target_id].len.get() as f32;
            let bin_length: f32 = *bin_width as f32; //txp_info[target_id].len.get() as f32;
            let start_bin: usize = (start_aln / bin_length) as usize;
            let end_bin: usize = (end_aln / bin_length) as usize;
            let coverage_probability: &Vec<f64> = &txp_info[target_id].coverage_prob;
            let num_bins: u32 = (tlen as f64 / bin_length as f64).ceil() as u32;

            let cov_prob: f64 = (start_bin..end_bin)
                .map(|i| {
                    let bin_floor_diff =
                        ((i + 1) as f32 * bin_length).floor() - (i as f32 * bin_length).floor();
                    let ceil_diff = tlen + 1.0 - (i as f32 * bin_length).ceil();

                    if bin_floor_diff <= 0.0
                        || ceil_diff <= 0.0
                        || !(coverage_probability[i]
                            .partial_cmp(&0.0)
                            .unwrap_or(std::cmp::Ordering::Less)
                            .is_ge())
                        || coverage_probability[i].is_nan()
                        || coverage_probability[i].is_infinite()
                    {
                        error!("bin_floor_diff: {:?}", bin_floor_diff);
                        error!("ceil_diff: {:?}", ceil_diff);
                        error!("coverage_probability[i]: {:?}", coverage_probability[i]);
                        error!("coverage_probability: {:?}", coverage_probability);
                        error!("tlen: {:?}", tlen);
                        error!("bin_length: {:?}", bin_length);
                        error!("num_bins: {:?}", num_bins);
                        panic!("first panic");
                    }

                    match (
                        i == start_bin,
                        i == end_bin,
                        end_bin == (num_bins - 1) as usize,
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
                let expected_cov_prob =
                    cov_prob / ((end_aln - start_aln) as f64 / bin_length as f64);
                let final_cov_prob = expected_cov_prob * (tlen as f64 / bin_length as f64);
                error!("cov_prob: {:?}", cov_prob);
                error!(
                    "length: {:?}",
                    ((end_aln - start_aln) as f64 / bin_length as f64)
                );
                error!("expected_cov_prob: {:?}", expected_cov_prob);
                error!("length2: {:?}", (tlen as f64 / bin_length as f64));
                error!("final_cov_prob: {:?}", final_cov_prob);
                error!("start_bin: {}, end_bin: {}", start_bin, end_bin);
                error!("start_aln: {}, end_aln: {}", start_aln, end_aln);
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
