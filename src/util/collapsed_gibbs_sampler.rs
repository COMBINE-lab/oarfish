//use rand::prelude::*;
use rand::rng;
use rand::rngs::StdRng;
use rand::{SeedableRng};
use rand::seq::SliceRandom;
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
        //let idx = categorical_selection(&probs_vec);
        let idx = probs_vec
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(i, _)| i)
                    .unwrap();
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


fn ln_safe(x: f64) -> f64 {
    const EPS: f64 = 1e-300; // small, to avoid -inf
    (x.max(EPS)).ln()
}



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
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    let mut read_candidate: Vec<Vec<usize>> = Vec::new();

    let mut counts = vec![0usize; tinfo.len()];

    let mut prob_idx: Vec<usize> = Vec::new();

    let mut assignments: Vec<usize> = Vec::new();

    let mut txp_probs = Vec::<f64>::new();

    let alpha_dir = 1.0;

    let mut current_likelihood = 0.0_f64;

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
        //let idx = categorical_selection(&probs_vec);
        let idx = probs_vec
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(i, _)| i)
                    .unwrap();

        prob_idx.push(idx);
        let tid = alns[idx].ref_id as usize;
        assignments.push(tid);
        counts[tid] += 1;

        //compute the likelihood for the initial assignment
        let prob = probs[idx] as f64;
        let cov_prob = if fops.model_coverage { coverage_probs[idx] } else { 1.0 };
        current_likelihood += prob.ln() + cov_prob.ln(); //+ (counts[tid] as f64 * (counts[tid] as f64 / total_weight).ln());
    }

    let n_reads = assignments.len() as f64; // <- true N
    for &c in &counts {
        if c > 0 {
            let cf = c as f64;
            current_likelihood += cf * (cf / n_reads).ln();
        }
    }



    read_candidate.push(assignments.clone());

    let mut initial_counts = counts.clone();
    let mut initial_assignments = assignments.clone();
    let mut initial_prob_idx = prob_idx.clone();

    let num_iteration = 1000;
    let epsilon = 1e-8_f64;
    for iteration in 1..num_iteration{
        
        let mut sweep_delta = 0.0_f64;
        let mut moved = 0usize;

        //obtain the assignmets that provides the maximum likelihood for each read
        for (read_num, (alns, probs, coverage_probs)) in eq_map.iter().enumerate() {

            //obtain the likelihood for a different alignment of a read
            let prob_idx_initial = initial_prob_idx[read_num];
            let tid_initial = initial_assignments[read_num];
            let mut best_delta = 0.0;
            let mut best_tid = tid_initial;
            let mut best_prob_idx = prob_idx_initial;
            

            for (align_num, (a, p, cp)) in izip!(alns, probs, coverage_probs).enumerate(){

                if align_num == prob_idx_initial{
                    continue;
                }

                let new_tid = a.ref_id as usize;
                
                let cu = initial_counts[tid_initial] as f64;
                let pu = probs[prob_idx_initial] as f64;
                let cpu = if fops.model_coverage { coverage_probs[prob_idx_initial] } else { 1.0 };

                let cv = initial_counts[new_tid] as f64;
                let pv = *p as f64;
                let cpv = if fops.model_coverage { *cp } else { 1.0 };

                let delta = ((cv + 1.0) * (cv + 1.0).ln()) +
                            ((cu - 1.0) * (cu - 1.0).ln()) -
                            (cv * cv.ln()) -
                            (cu * cu.ln()) + 
                            (pv.ln() + cpv.ln()) - 
                            (pu.ln() + cpu.ln());

                //let delta = ((cv + 1.0) * ln_safe(cv + 1.0)) +
                //            ((cu - 1.0) * ln_safe(cu - 1.0)) -
                //            (cv * ln_safe(cv)) -
                //            (cu * ln_safe(cu)) + 
                //            (ln_safe(pv) + ln_safe(cpv)) - 
                //            (ln_safe(pu) + ln_safe(cpu));

                //let delta = ((alpha_dir + cv).ln()) -
                //            ((alpha_dir + cu - 1.0).ln()) + 
                //            ((pv).ln() + (cpv).ln()) - 
                //            ((pu).ln() + (cpu).ln());

            
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
                current_likelihood += best_delta;
                sweep_delta += best_delta;
                moved += 1;
            }

        }

        eprintln!("number of itration: {iteration}             sweep_delta: {sweep_delta},          likelihood: {current_likelihood},           #moved: {moved}");

        read_candidate.push(initial_assignments.clone());
        if sweep_delta <= epsilon {
            break;
        }

    }

    //hard assignment is the result of the last iteration of the above loop
    let hard_assignment = initial_assignments;

    //(hard_assignment, read_candidate)

    //let m = tinfo.len();
    //let n = read_candidate[0].len();
    //let hard_assignment = (0..n)
    //    .into_par_iter()
    //    .map(|read_idx| {
    //        let mut counts_assign = vec![0usize; m];
    //        for it in read_candidate.iter() {
    //            counts_assign[it[read_idx]] += 1;
    //        }
    //        // Argmax
    //        let (mut best_tid, mut best_cnt) = (0, 0);
    //        for (tid, &c) in counts_assign.iter().enumerate() {
    //            if c > best_cnt {
    //                best_tid = tid;
    //                best_cnt = c;
    //            }
    //        }
    //        best_tid
    //    }).collect();
//
    (hard_assignment, read_candidate)

}




//###########################################################################################################
//compute the results for different orders of the reads and choose the one with the highest likelihood values

fn read_orders_indices(num_reads: usize, num_orders: usize, seed: u64) -> Vec<Vec<usize>> {

    if num_orders == 0 {
        return Vec::new();
    }

    let base_idx: Vec<usize> = (0..num_reads).collect();

    // Build entries 1..num_orders in parallel
    let rest: Vec<Vec<usize>> = (1..num_orders)
        .into_par_iter()
        .map(|i| {
            let mut idx = base_idx.clone();
            let mut rng = StdRng::seed_from_u64(
                seed ^ (i as u64).wrapping_mul(0x9E3779B97F4A7C15)
            );
            idx.shuffle(&mut rng);
            idx
        })
        .collect();

    // first is the original order, the rest are shuffled
    let mut out = Vec::with_capacity(num_orders);
    out.push(base_idx.clone());
    out.extend(rest);

    assert!(
        out.len() == num_orders,
        "expected {} orders, got {}",
        num_orders,
        out.len()
    );
    out
}

fn compute_delta_function(
    read_orders: &Vec<usize>, 
    num_iteration: usize, 
    em_info: &EMInfo, 
    initial_counts: &Vec<usize>, 
    initial_prob_idx: &Vec<usize>, 
    initial_assignments: &Vec<usize>, 
    initial_likelihood: &f64,
    epsilon: f64) 
    -> (Vec<usize>, f64) {

    let mut counts: Vec<usize> = initial_counts.to_vec();
    let mut prob_idx: Vec<usize> = initial_prob_idx.to_vec();
    let mut assignments: Vec<usize> = initial_assignments.to_vec();
    let mut current_likelihood: f64 = *initial_likelihood;

    let per_read: Vec<_> = em_info.eq_map.iter().collect();

    for iteration in 0..num_iteration{
        
        let mut sweep_delta = 0.0_f64;
        let mut moved = 0_usize;

        //obtain the assignmets that provides the maximum likelihood for each read
        for &read_idx in read_orders.iter(){

            //obtain the likelihood for a different alignment of a read
            let prob_idx_initial = prob_idx[read_idx];
            let tid_initial = assignments[read_idx];

            let mut best_delta = 0.0;
            let mut index_delta: Vec<(usize, f64)> = vec![];
            let mut best_tid = tid_initial;
            let mut best_prob_idx = prob_idx_initial;

            let (alns, probs, covs) = per_read[read_idx];

            for (align_num, (a, p, cp)) in izip!(alns, probs, covs).enumerate(){
                if align_num == prob_idx_initial{
                    index_delta.push((align_num, 0.0));
                    continue;
                }
                let new_tid = a.ref_id as usize;
                let cu = counts[tid_initial] as f64;
                let pu = probs[prob_idx_initial] as f64;
                let cpu = if em_info.eq_map.filter_opts.model_coverage { covs[prob_idx_initial] } else { 1.0 };
                let cv = counts[new_tid] as f64;
                let pv = *p as f64;
                let cpv = if em_info.eq_map.filter_opts.model_coverage { *cp } else { 1.0 };

                let delta = ((cv + 1.0) * (cv + 1.0).ln()) +
                            ((cu - 1.0) * (cu - 1.0).ln()) -
                            (cv * cv.ln()) -
                            (cu * cu.ln()) + 
                            (pv.ln() + cpv.ln()) - 
                            (pu.ln() + cpu.ln());
                
                if !delta.is_finite() {
                    panic!("delta became non-finite at iteration");
                }
                index_delta.push((align_num, delta));
                //if delta > best_delta{
                //    best_delta = delta;
                //    best_tid = new_tid;
                //    best_prob_idx = align_num;
                //}
            }
            //check if it has positive or zero values in the delta vector
            let has_pos  = index_delta.iter().any(|&(_,d)| d > 0.0);
            let has_zero = index_delta.iter().any(|&(_,d)| d == 0.0);
            //use categorical function between those with poisitive move
            if !has_pos && !has_zero {
                continue;
            }
            
            let delta_vec: Vec<f64> = if has_pos {
                index_delta.iter().map(|&(_, d)| if d > 0.0 { d } else { 0.0 }).collect()
            } else {
                let zeros = index_delta.iter().filter(|&&(_, d)| d == 0.0).count() as f64;
                index_delta.iter().map(|&(_, d)| if d == 0.0 { 1.0 / zeros } else { 0.0 }).collect()
            };

            let probs = normalize_probability(&delta_vec);           // &[f64] -> Vec<f64>
            let idx  = categorical_selection(&probs); 

            
            let new_tid = alns[index_delta[idx].0].ref_id as usize;
            best_delta = index_delta[idx].1;
            best_tid = new_tid;
            best_prob_idx = index_delta[idx].0;
            

            // Update assignment and counts
            if best_prob_idx != prob_idx_initial && best_delta > 0.0 {
                assignments[read_idx] = best_tid;
                counts[best_tid] += 1;
                counts[tid_initial] -= 1;
                prob_idx[read_idx] = best_prob_idx;
                current_likelihood += best_delta;
                sweep_delta += best_delta;
                moved += 1;
            }

        }

        //read_candidate.push(initial_assignments.clone());
        if sweep_delta <= epsilon {
            break;
        }

    }

    (assignments, current_likelihood)
}




pub fn greedy_likelihood_assignment(
    em_info: &EMInfo,
    em_counts: &[f64],
) -> Vec<usize> {

    eprintln!("greedy likelihood model");
    eprintln!("running the posterior probability for initial counts");

    let eq_map = &em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = &em_info.txp_info;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    let mut counts = vec![0usize; tinfo.len()];
    let mut prob_idx: Vec<usize> = Vec::new();
    let mut assignments: Vec<usize> = Vec::new();
    let mut txp_probs = Vec::<f64>::new();
    let mut current_likelihood = 0.0_f64;

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
        let idx = probs_vec
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(i, _)| i)
                    .unwrap();

        prob_idx.push(idx);
        let tid = alns[idx].ref_id as usize;
        assignments.push(tid);
        counts[tid] += 1;

        //compute the likelihood for the initial assignment
        let prob = probs[idx] as f64;
        let cov_prob = if fops.model_coverage { coverage_probs[idx] } else { 1.0 };
        current_likelihood += prob.ln() + cov_prob.ln(); //+ (counts[tid] as f64 * (counts[tid] as f64 / total_weight).ln());
    }

    let n_reads = assignments.len() as f64; // <- true N
    for &c in &counts {
        if c > 0 {
            let cf = c as f64;
            current_likelihood += cf * (cf / n_reads).ln();
        }
    }


    //compute the different orders of read indices
    let num_reads = n_reads as usize;
    let num_orders = 50;
    let seed = 49;

    let read_orders: Vec<Vec<usize>> = read_orders_indices(num_reads, num_orders, seed);

    //compute the assignments and likelihood for each indices of the reads
    let num_iteration: usize = 50;
    let epsilon = 1e-8_f64;

    let results: Vec<(Vec<usize>, f64)> = read_orders
        .par_iter()                            
        .map(|indices_order| {
            compute_delta_function(
                indices_order,
                num_iteration,
                &em_info,
                &counts,
                &prob_idx,
                &assignments,
                &current_likelihood,
                epsilon,
            )
        })
        .collect();

    let (read_candidate, likelihoods): (Vec<Vec<usize>>, Vec<f64>) = results.into_iter().unzip();
    
    //update the index that have the highest likelihood value
    //let idx = likelihoods
    //                .iter()
    //                .enumerate()
    //                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
    //                .map(|(i, _)| i)
    //                .unwrap();

    //hard assignment is the result of the last iteration of the above loop
    //let hard_assignment = read_candidate[idx].clone();

    //hard assignment is obtained based on the maximum occurences in results for all read indices order
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

    hard_assignment

}