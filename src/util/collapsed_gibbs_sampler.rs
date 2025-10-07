//use rand::prelude::*;
use rand::rng;
use rayon::prelude::*;
use rand::distr::Distribution;    
use rand_distr::weighted::WeightedIndex;
use itertools::izip;
use crate::util::oarfish_types::EMInfo;
use crate::TranscriptInfo;


//write the categorial function to choose the transcript index based on the probability
fn categorical_selection(weights: &[f64]) -> usize {
    let dist = WeightedIndex::new(weights).expect("invalid weights (must be non-negative and not all zero)");

    let mut rng = rng();

    let index = dist.sample(&mut rng);

    index
}

//normalize the probabilities
fn normalize_probability(prob: &[f64]) -> Vec<f64> {
    let nprob_sum: f64 = prob.iter().copied().sum();

    let sum_normalize_probs = if nprob_sum > 0.0 { nprob_sum } else { 1.0 };

    let normalized_prob: Vec<f64> = prob.iter()
        .map(|&x| x / sum_normalize_probs)
        .collect();

    normalized_prob
}


pub fn collapsed_sequential(
    em_info: &EMInfo
) -> (Vec<usize>, Vec<Vec<usize>>) {

    let eq_map = &em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = &em_info.txp_info;
    //let max_iter = em_info.max_iter;
    //let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    //constant for alpha value in the equation
    let alpha_dir = 1.0;
    let maximum_gibbs_iteration = 50;

    let mut read_candidate: Vec<Vec<usize>> = Vec::new();

    let mut counts = vec![0usize; tinfo.len()];

    let mut assignments: Vec<usize> = Vec::new();

    for (alns, probs, coverage_probs) in eq_map.iter(){
        let result: Vec<f64> = if !fops.model_coverage {
            probs.iter().map(|&p| p as f64).collect()
        } else {
            probs.iter().zip(coverage_probs.iter())
                 .map(|(&p, &c)| (p as f64) * c)
                 .collect()
        };
        let probs_vec = normalize_probability(&result);
        let idx = categorical_selection(&probs_vec);
        let tid = alns[idx].ref_id as usize;
        assignments.push(tid);
        counts[tid] += 1;
    }

    read_candidate.push(assignments.clone());

    for _iteration in 1..maximum_gibbs_iteration{
        
        eprintln!("the gibbs itration: {_iteration}");
        let mut initial_counts = counts.clone();
        let mut initial_assignments = assignments.clone();

        for (i, (alns, probs, coverage_probs)) in eq_map.iter().enumerate(){
            // Leave-one-out: remove current assignment for this read
            let cur_tid = initial_assignments[i];
            initial_counts[cur_tid] = initial_counts[cur_tid]
                .checked_sub(1)
                .expect("counts underflow");

            // Collapsed weights for this read: w_m ∝ P(r_n|m) * (alpha_dir + C_m^{(-n)})
            let mut weights = Vec::with_capacity(alns.len());
            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
        
                let cm = initial_counts[target_id] as f64; // C_m^{(-n)}
                weights.push(prob * cov_prob * (alpha_dir + cm));
            }

            // Normalize then draw a categorical sample
            let probs = normalize_probability(&weights);
            let idx = categorical_selection(&probs);
            let new_tid = alns[idx].ref_id as usize;

            // Update assignment and counts immediately (exact Gibbs)
            initial_assignments[i] = new_tid;
            initial_counts[new_tid] += 1;
        }

        read_candidate.push(initial_assignments);

    }

    //eprintln!("read_candidate0: {:?}", read_candidate[0]);
    //eprintln!("read_candidate49: {:?}", read_candidate[49]);

    //hard assignment such that choosing the transcript with highest frequesncy for each read
    let m = tinfo.len();
    let n = read_candidate[0].len();
    let hard_assignment = (0..n)
        .into_par_iter()
        .map(|read_idx| {
            let mut counts_assign = vec![0usize; m];
            for it in read_candidate.iter() {
                counts_assign[it[read_idx]] += 1;
            }
            // Argmax
            let (mut best_tid, mut best_cnt) = (0, 0);
            for (tid, &c) in counts_assign.iter().enumerate() {
                if c > best_cnt {
                    best_tid = tid;
                    best_cnt = c;
                }
            }
            best_tid
        }).collect();

    (hard_assignment, read_candidate)

}
