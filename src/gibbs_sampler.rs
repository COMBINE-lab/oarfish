use core::f64;

use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use crate::util::probs::LogSpace;
use itertools::*;
use rand::Rng;
use rand::distr::{Distribution, Uniform};
use rand_distr::weighted::WeightedIndex;
use rv::prelude::*;

pub fn update_chain<'a, DFn, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a>(
    eq_map_iter: I,
    tinfo: &[TranscriptInfo],
    model_coverage: bool,
    density_fn: DFn,
    chain_state: &mut ChainState,
) where
    DFn: Fn(usize, usize) -> f64,
{
    let a_dir = 1f64;
    let a_act = 2f64;
    let _b_act = 1f64;

    let mut rng = rand::rng();
    let mut txp_ids = Vec::new();
    let mut txp_probs = Vec::new();

    let tot_reads = (chain_state.assigned_ids.len() - 1) as f64;
    let num_txps = chain_state.counts.len() as f64;

    for (read_id, (alns, probs, coverage_probs)) in eq_map_iter.enumerate() {
        txp_ids.clear();
        txp_probs.clear();

        let assigned_tid = chain_state.assigned_ids[read_id];
        // remove this read from this transcript
        chain_state.counts[assigned_tid as usize] -= 1;
        let mut denom = 0.0_f64;

        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let txp_len = tinfo[target_id].lenf as usize;
            let aln_len = a.alignment_span() as usize;

            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            let dens_prob = density_fn(txp_len, aln_len);

            let sprob = (prob * cov_prob * dens_prob)
                * (a_act + tot_reads)
                * ((a_dir + chain_state.counts[a.ref_id as usize] as f64)
                    / (num_txps * a_dir + tot_reads));
            txp_ids.push(a.ref_id);
            txp_probs.push(sprob);
            denom += sprob;
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            let cat = WeightedIndex::new(&txp_probs).unwrap();
            let s = cat.sample(&mut rng);
            let new_assigned_tid = txp_ids[s];
            // set the assignment for this read
            chain_state.assigned_ids[read_id] = new_assigned_tid;
            // update the count accordingly
            chain_state.counts[new_assigned_tid as usize] += 1;
        } else {
            chain_state.counts[assigned_tid as usize] += 1;
        }
    }
}

fn run_sampler<'a, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a, F: Fn() -> I>(
    em_info: &'a EMInfo,
    ml_abundances: &[f64], // abundances obtained from the EM algorithm
    make_iter: F,
) {
    let mut cs = init_counts_and_assignments(em_info, ml_abundances, make_iter);
    let make_iter = || em_info.eq_map.iter();
    let tinfo = &em_info.txp_info;
    let fops = &em_info.eq_map.filter_opts;
    let model_coverage = fops.model_coverage;
    let density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.,
        }
    };

    let mut best_log_assignments: Vec<u32>;
    let mut best_log_likelihood = crate::util::probs::MIN_LOG_P.get_linear();
    for _ in 0..100 {
        update_chain(make_iter(), tinfo, model_coverage, density_fn, &mut cs);
        let ll = cs.log_likelihood();
        if ll > best_log_likelihood {
            best_log_likelihood = ll;
            best_log_assignments = cs.assigned_ids.clone();
        }
    }
}

struct ChainState {
    pub counts: Vec<u64>, // the current (integer) assignment of counts to each transcript
    pub assigned_ids: Vec<u32>, // the current assignment of each read to a transcript of origin
}

impl ChainState {
    pub fn log_likelihood(&self) -> f64 {
        let mut ll = crate::util::probs::MIN_LOG_P;
        let tot: u64 = self.counts.iter().sum();
        let norm = LogSpace::new_from_linear(1f64 / (tot as f64));
        for c in &self.counts {
            let log_c = LogSpace::new_from_linear(*c as f64);
            // get the actual f64 out so we can apply the identity that
            // log_b(x^p) = p * log_b(x)
            let ll_contrib = log_c.get_ln() * (log_c * norm).get_ln();
            ll += LogSpace::new_from_ln(ll_contrib);
        }

        ll.get_ln()
    }
}

fn init_counts_and_assignments<
    'a,
    I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a,
    F: Fn() -> I,
>(
    em_info: &'a EMInfo,
    ml_abundances: &[f64], // abundances obtained from the EM algorithm
    make_iter: F,
) -> ChainState {
    let mut counts = vec![0; em_info.num_txps()];
    let mut assigned_ids: Vec<u32> = Vec::new();
    let tinfo = &em_info.txp_info;
    let fops = &em_info.eq_map.filter_opts;
    let model_coverage = fops.model_coverage;
    let density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.,
        }
    };

    let mut rng = rand::rng();
    let mut txp_ids = Vec::new();
    let mut txp_probs = Vec::new();
    for (alns, probs, coverage_probs) in make_iter() {
        txp_ids.clear();
        txp_probs.clear();
        let mut denom = 0.0_f64;
        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let txp_len = tinfo[target_id].lenf as usize;
            let aln_len = a.alignment_span() as usize;

            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            let dens_prob = density_fn(txp_len, aln_len);

            let sprob = ml_abundances[target_id] * prob * cov_prob * dens_prob;
            txp_ids.push(a.ref_id);
            txp_probs.push(sprob);
            denom += sprob;
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            let cat = WeightedIndex::new(&txp_probs).unwrap();
            let s = cat.sample(&mut rng);
            assigned_ids.push(txp_ids[s]);
        }
    }

    ChainState {
        counts,
        assigned_ids,
    }
}
