use crate::util::oarfish_types::{TranscriptInfo, InMemoryAlignmentStore};
use noodles_sam as sam;
use noodles_bam as bam;

pub fn normalize_read_probs(store: &mut InMemoryAlignmentStore, txp_info: &Vec<TranscriptInfo>) {

    let mut normalize_probs_temp: Vec<f64> = vec![];
    let mut normalized_coverage_prob: Vec<f64> = vec![];

    //iterate over all alignments in the bam file
    for (alns, _as_probs, _coverage_prob) in store.iter() {
        //iterate over the alignments of a read
        for a in alns.iter() {
            let target_id = a.ref_id as usize;
            let start_aln = a.start as usize;
            let end_aln = a.end as usize;
            let cov_prob = (txp_info[target_id].coverage_prob[end_aln] - txp_info[target_id].coverage_prob[start_aln]) as f64;

            if cov_prob.is_nan() || cov_prob.is_infinite() {

                panic!("Error: Invalid result. normalize_read_probs function.");
            }
            normalize_probs_temp.push(cov_prob);
        }
        let sum_normalize_probs_temp: f64 = if normalize_probs_temp.iter().sum::<f64>() > 0.0 {normalize_probs_temp.iter().sum()} else {1.0};
        let normalized_prob_section: Vec<f64> = normalize_probs_temp.iter().map(|&prob| prob/sum_normalize_probs_temp).collect();

        normalized_coverage_prob.extend(normalized_prob_section);
        normalize_probs_temp.clear();
    }
    store.coverage_probabilities = normalized_coverage_prob;
}