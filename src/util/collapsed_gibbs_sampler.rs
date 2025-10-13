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


pub fn collapsed_sequential_prior_prob(
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





pub fn collapsed_sequential_posterior_prob(
    em_info: &EMInfo,
    em_counts: &[f64],
) -> (Vec<usize>, Vec<Vec<usize>>) {

    eprintln!("running the posterior probability for initial counts");

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

    let mut txp_probs = Vec::<f64>::new();

    for (alns, probs, coverage_probs) in eq_map.iter() {

        let mut denom = 0.0_f64;

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
            denom += em_counts[target_id] * prob * cov_prob;
        }

        txp_probs.clear();

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
            let nprob = ((em_counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
            txp_probs.push(nprob);
        }

        let probs_vec = normalize_probability(&txp_probs);
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






//pub fn max_likelihood_posterior_prob(
//    em_info: &EMInfo,
//    em_counts: &[f64],
//) -> (Vec<usize>, Vec<Vec<usize>>) {
//
//    eprintln!("running the posterior probability for initial counts");
//
//    let eq_map = &em_info.eq_map;
//    let fops = &eq_map.filter_opts;
//    let tinfo: &[TranscriptInfo] = &em_info.txp_info;
//    //let max_iter = em_info.max_iter;
//    //let total_weight: f64 = eq_map.num_aligned_reads() as f64;
//
//    //constant for alpha value in the equation
//    let num_iteration = 50;
//
//    let mut read_candidate: Vec<Vec<usize>> = Vec::new();
//
//    let mut counts = vec![0usize; tinfo.len()];
//
//    let mut prob_idx: Vec<usize> = Vec::new();
//
//    let mut assignments: Vec<usize> = Vec::new();
//
//    let mut txp_probs = Vec::<f64>::new();
//
//    let mut prev_likelihood = 0.0_f64;
//    let mut curr_likelihood = 0.0_f64;
//
//    //obtain the initial counts and assignments
//    for (alns, probs, coverage_probs) in eq_map.iter() {
//
//        let mut denom = 0.0_f64;
//
//        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
//            let target_id = a.ref_id as usize;
//            let prob = *p as f64;
//            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
//            denom += em_counts[target_id] * prob * cov_prob;
//        }
//
//        txp_probs.clear();
//
//        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
//            let target_id = a.ref_id as usize;
//            let prob = *p as f64;
//            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
//            let nprob = ((em_counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
//            txp_probs.push(nprob);
//        }
//
//        let probs_vec = normalize_probability(&txp_probs);
//        let idx = categorical_selection(&probs_vec);
//        prob_idx.push(idx);
//        let tid = alns[idx].ref_id as usize;
//        assignments.push(tid);
//        counts[tid] += 1;
//    }
//
//    read_candidate.push(assignments.clone());
//
//    let mut initial_counts = counts.clone();
//    let mut initial_assignments = assignments.clone();
//    let mut initial_prob_idx = prob_idx.clone();
//
//    for _iteration in 1..num_iteration{
//        
//        eprintln!("number of itration: {_iteration}");
//
//        //obtain the assignmetn that provides the maximum likelihood for each read
//        for read_num in 0..len(assignments){
//            //obtain the likelihood for a different alignment of a read
//            let prob_idx_initial = initial_prob_idx[read_num];
//            let tid_initial = assignments[read_num];
//            let final_tid = tid_initial;
//
//            
//
//            for align_num in 0..len(eq_map.alignments[read_num]){
//
//                let new_tid = eq_map.alignments[read_num][align_num].ref_id as usize;
//                updated_counts[tid_initial] = initial_counts[tid_initial].checked_sub(1).expect("counts underflow");
//                updated_counts[new_tid] = updated_counts[new_tid].saturating_add(1);
//
//                for (i, (alns, probs, coverage_probs)) in eq_map.iter().enumerate(){
//                    let index = if read_num == i {
//                        align_num
//                    } else {
//                        initial_prob_idx[i]
//                    };
//                    let prob = *probs[index] as f64;
//                    let cov_prob = if fops.model_coverage { *coverage_probs[index] } else { 1.0 };
//                    curr_likelihood += np.log(prob) + np.log(cov_prob) + np.log(updated_counts[new_tid]);
//                }
//            
//                if curr_likelihood > prev_likelihood{
//                    prev_likelihood = curr_likelihood;
//                    final_tid = new_tid;
//                    final_prob_idx = align_num;
//                }
//            }
//
//            // Update assignment and counts
//            initial_assignments[read_num] = final_tid;
//            initial_counts[final_tid] += 1;
//            initial_prob_idx[final_tid] = final_prob_idx;
//
//        }
//
//        read_candidate.push(initial_assignments);
//
//    }
//
//    //hard assignment is the result of the last iteration of the above loop
//    let hard_assignment = initial_assignments;
//
//    (hard_assignment, read_candidate)
//
//}







pub fn optimized_max_likelihood_posterior_prob(
    em_info: &EMInfo,
    em_counts: &[f64],
) -> (Vec<usize>, Vec<Vec<usize>>) {

    eprintln!("optimized maximum likelihood model");
    eprintln!("running the posterior probability for initial counts");

    let eq_map = &em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = &em_info.txp_info;
    //let max_iter = em_info.max_iter;
    //let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    let mut read_candidate: Vec<Vec<usize>> = Vec::new();

    let mut counts = vec![0usize; tinfo.len()];

    let mut prob_idx: Vec<usize> = Vec::new();

    let mut assignments: Vec<usize> = Vec::new();

    let mut txp_probs = Vec::<f64>::new();

    //let mut current_likelihood = 0.0_f64;

    //let total_n = eq_map.alignments.len() as f64;


    //obtain the initial counts and assignments
    for (alns, probs, coverage_probs) in eq_map.iter() {

        let mut denom = 0.0_f64;

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
            denom += em_counts[target_id] * prob * cov_prob;
        }

        txp_probs.clear();

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
            let nprob = ((em_counts[target_id] * prob * cov_prob) / denom).clamp(0.0, 1.0);
            txp_probs.push(nprob);

        }

        let probs_vec = normalize_probability(&txp_probs);
        //??????????????
        //choose the initial assinments based on the categorical prob or choosing the one with highest prob?
        let idx = categorical_selection(&probs_vec);
        prob_idx.push(idx);
        let tid = alns[idx].ref_id as usize;
        assignments.push(tid);
        counts[tid] += 1;

        //compute the likelihood for the initial assignment
        //let prob = probs[idx] as f64;
        //let cov_prob = if fops.model_coverage { coverage_probs[idx] } else { 1.0 };
        //current_likelihood += prob.ln() + cov_prob.ln() + (counts[tid] as f64 * (counts[tid] as f64 / total_n).ln());
    }

    read_candidate.push(assignments.clone());

    let mut initial_counts = counts.clone();
    let mut initial_assignments = assignments.clone();
    let mut initial_prob_idx = prob_idx.clone();

    let num_iteration = 50;
    let epsilon = 1e-8_f64;

    for _iteration in 1..num_iteration{
        
        eprintln!("number of itration: {_iteration}");
        let mut sweep_delta = 0.0_f64;

        //obtain the assignmets that provides the maximum likelihood for each read
        for (read_num, (alns, probs, coverage_probs)) in eq_map.iter().enumerate() {

            //obtain the likelihood for a different alignment of a read
            let prob_idx_initial = initial_prob_idx[read_num];
            let tid_initial = assignments[read_num];
            let mut best_delta = 0.0;
            let mut best_tid = tid_initial;
            let mut best_prob_idx = prob_idx_initial;
            

            for (align_num, (a, _p, _cp)) in izip!(alns, probs, coverage_probs).enumerate(){

                if align_num == prob_idx_initial{
                    continue;
                }

                let new_tid = a.ref_id as usize;
                
                let cu = initial_counts[tid_initial] as f64;
                let cv = initial_counts[new_tid] as f64;

                let delta = ((cv + 1.0) * (cv + 1.0).ln()) +
                            ((cu - 1.0) * (cu - 1.0).ln()) -
                            (cv * (cv).ln()) -
                            (cu * (cu).ln());

            
                if delta > best_delta{
                    best_delta = delta;
                    best_tid = new_tid;
                    best_prob_idx = align_num;
                }
            }

            // Update assignment and counts
            if best_prob_idx != prob_idx_initial && best_delta > 0.0 {
                initial_assignments[read_num] = best_tid;
                initial_counts[best_tid] += 1;
                initial_counts[tid_initial] -= 1;
                initial_prob_idx[read_num] = best_prob_idx;
                //current_likelihood += best_delta;
                sweep_delta += best_delta;
            }

        }

        read_candidate.push(initial_assignments.clone());
        if sweep_delta <= epsilon {
            break;
        }

    }

    //hard assignment is the result of the last iteration of the above loop
    let hard_assignment = initial_assignments;

    (hard_assignment, read_candidate)

}

