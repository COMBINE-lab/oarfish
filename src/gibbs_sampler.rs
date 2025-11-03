use core::f64;

use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use crate::util::probs::LogSpace;
use itertools::*;
use rand::Rng;
use rand::distr::Distribution;
use rand_distr::weighted::WeightedIndex;
use std::default::Default;
use tracing::{info, warn};

/// The parameters to be used for the Gibbs sampler. The
/// default values are taken from the BitSeq[1] paper, though
/// currently the "noise" transcript estimation is not performed.
/// reference:
/// 1) Glaus, P., Honkela, A., & Rattray, M. (2012).
///    Identifying differentially expressed transcripts from RNA-seq data with biological variation.
///    Bioinformatics, 28(13), 1721-1728.
pub(crate) struct SamplerParams {
    pub a_dir: f64,
    pub a_act: f64,
    pub b_act: f64,
    pub niter: u64,
}

impl Default for SamplerParams {
    fn default() -> Self {
        Self {
            a_dir: 1f64,
            a_act: 2f64,
            b_act: 2f64,
            niter: 100u64,
        }
    }
}

impl SamplerParams {
    pub fn new_bitseq_with_iter(niter: u64) -> Self {
        Self {
            niter,
            ..Default::default()
        }
    }
}

/// make a pass over all of the alignments and, given the current state of the chain, sample
/// a new assignment for each read based on the complete conditional of the rest of the
/// assignments.
///
/// This implements the collapsed Gibbs sampler from BitSeq (though currently without a noise
/// transcript term).
///
/// The return value is the log-likelihood of the newly-sampled chain.
fn update_chain<'a, DFn, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a>(
    sampler_params: &SamplerParams,
    eq_map_iter: I,
    tinfo: &[TranscriptInfo],
    model_coverage: bool,
    density_fn: DFn,
    chain_state: &mut ChainState,
) where
    DFn: Fn(usize, usize) -> f64,
{
    let a_dir = sampler_params.a_dir;
    let a_act = sampler_params.a_act;
    let _b_act = sampler_params.b_act;

    let mut rng = rand::rng();
    let mut txp_ids = Vec::new();
    let mut txp_probs = Vec::new();
    let mut txp_log_cond_probs = Vec::new();

    // -1 because we will use this to normalize each step in the chain
    // where we have removed 1 read
    let tot_reads = (chain_state.assigned_ids.len() - 1) as f64;
    // total number of transcripts in the transcriptome (not just expressed)
    let num_txps = chain_state.counts.len() as f64;

    let mut aln_idx = 0usize;
    for (read_id, (alns, probs, coverage_probs)) in eq_map_iter.enumerate() {
        // if this is a uniquely-aligned read, don't bother with
        // all of this
        if alns.len() == 1 {
            chain_state.inc_sample_counts(aln_idx);
            aln_idx += 1;
            continue;
        }
        // clear out our temporary storage vectors
        txp_ids.clear();
        txp_probs.clear();
        txp_log_cond_probs.clear();

        // the ID assigned to this transcript before the reassignment
        let assigned_tid: usize = chain_state.assigned_ids[read_id] as usize;
        // remove this read from this transcript's counts
        chain_state.counts[assigned_tid] -= 1;
        // the normalizing constant
        let mut denom = 0.0_f64;
        let mut idx = 0;
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
            let cond_prob = prob * cov_prob * dens_prob;
            let sprob = cond_prob
                * (a_act + tot_reads)
                * ((a_dir + chain_state.counts[a.ref_id as usize] as f64)
                    / (num_txps * a_dir + tot_reads));

            txp_log_cond_probs.push(LogSpace::new_from_linear(cond_prob));
            txp_ids.push(a.ref_id);
            txp_probs.push(sprob);
            denom += sprob;
            idx += 1;
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            let cat = WeightedIndex::new(&txp_probs).unwrap();
            //let s = max_sprob_idx;
            let s = cat.sample(&mut rng);
            let new_assigned_tid = txp_ids[s] as usize;
            chain_state.inc_sample_counts(aln_idx + s);
            // set the assignment for this read
            chain_state.assigned_ids[read_id] = new_assigned_tid as u32;
            // update the count accordingly
            chain_state.counts[new_assigned_tid] += 1;
            // update the conditional probability
            chain_state.log_cond_probs[read_id] = txp_log_cond_probs[s];
        } else {
            // the read doesn't move
            chain_state.counts[assigned_tid] += 1;
            // conditional probability remains the same
        }
        aln_idx += idx;
    }
}

fn greedy_pass<'a, DFn, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a>(
    sampler_params: &SamplerParams,
    eq_map_iter: I,
    tinfo: &[TranscriptInfo],
    model_coverage: bool,
    density_fn: DFn,
    chain_state: &mut ChainState,
) -> f64
where
    DFn: Fn(usize, usize) -> f64,
{
    let a_dir = sampler_params.a_dir;
    let a_act = sampler_params.a_act;
    let _b_act = sampler_params.b_act;

    let mut rng = rand::rng();

    let tot_reads = (chain_state.assigned_ids.len()) as f64;
    let log_tot_reads = tot_reads.ln();
    // total number of transcripts in the transcriptome (not just expressed)
    let num_txps = chain_state.counts.len() as f64;

    let mut total_delta = 0.;
    let mut aln_idx = 0usize;
    for (read_id, (alns, probs, coverage_probs)) in eq_map_iter.enumerate() {
        // if this is a uniquely-aligned read, don't bother with
        // all of this
        if alns.len() == 1 {
            //assignments[read_id] = chain_state.assigned_ids[read_id];
            continue;
        }

        // the ID assigned to this transcript before the reassignment
        let assigned_tid: usize = chain_state.assigned_ids[read_id] as usize;
        let assigned_log_cond_prob = chain_state.log_cond_probs[read_id];

        let mut best_delta = 0.;
        let mut best_log_cond_prob = assigned_log_cond_prob;
        let mut best_tid = assigned_tid;

        let mut denom = 0.0_f64;
        let mut idx = 0;
        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;

            // don't need to consider this case (not moving the read)
            if target_id == assigned_tid {
                continue;
            }

            let txp_len = tinfo[target_id].lenf as usize;
            let aln_len = a.alignment_span() as usize;

            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            let dens_prob = density_fn(txp_len, aln_len);
            let cond_prob = prob * cov_prob * dens_prob;

            let log_cond_prob = LogSpace::new_from_linear(cond_prob);

            let c_u = chain_state.counts[assigned_tid] as f64;
            let c_v = chain_state.counts[target_id] as f64;
            let log_c_u = c_u.ln();
            let log_c_v = c_v.ln();
            let log_theta_u = log_c_u - log_tot_reads;
            let log_theta_v = log_c_v - log_tot_reads;
            let log_theta_u_prime = (c_u - 1.).ln() - log_tot_reads;
            let log_theta_v_prime = (c_v + 1.).ln() - log_tot_reads;
            let delta = (-assigned_log_cond_prob.get_ln() - c_u * log_theta_u - c_v * log_theta_v)
                + (log_cond_prob.get_ln()
                    + (c_u - 1.) * log_theta_u_prime
                    + (c_v + 1.) * log_theta_v_prime);
            if delta > best_delta {
                best_delta = delta;
                best_log_cond_prob = log_cond_prob;
                best_tid = target_id;
            }
        }

        if best_delta > 0. {
            println!(
                "moving {} from txp {} to {}; delta = {}",
                read_id, assigned_tid, best_tid, best_delta
            );
            chain_state.counts[best_tid] += 1;
            chain_state.counts[assigned_tid] -= 1;
            chain_state.log_cond_probs[read_id] = best_log_cond_prob;
            chain_state.assigned_ids[read_id] = best_tid as u32;
        }
        total_delta += best_delta;
    }

    total_delta
}

/// run the sampler and return the highest log-likelihood assignments
pub fn run_collapsed_sampler<
    'a,
    I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a,
    F: Fn() -> I,
>(
    sampler_params: SamplerParams,
    em_info: &'a EMInfo,
    ml_abundances: &[f64], // abundances obtained from the EM algorithm
    make_iter: F,
) -> Vec<u32> {
    // use the EM abundance estimates to initialize discrete assignments
    let with_frequencies = true;
    let mut cs = init_counts_and_assignments(em_info, ml_abundances, make_iter, with_frequencies);
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

    // we'll keep track of the best log-likelihood we've seen so far as
    // well as the assignments that gave rise to it.
    let mut best_log_assignments: Vec<u32> = Vec::new();
    let mut best_log_likelihood = f64::NEG_INFINITY;
    let bar = indicatif::ProgressBar::new(sampler_params.niter);
    bar.set_style(
        indicatif::ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:20.green/blue} {msg} {human_pos:>12}",
        )
        .unwrap()
        .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
    );
    bar.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

    let mut best_iter = 0;
    let mut prev_ll = cs.log_likelihood();

    let mut best_ll = prev_ll;

    let mut cs_backup = cs.clone();
    let mut best_assign = cs.assigned_ids.clone();
    let mut resample_prob = 0.5_f64;

    for iter in 0..sampler_params.niter {
        bar.inc(1);
        update_chain(
            &sampler_params,
            make_iter(),
            tinfo,
            model_coverage,
            density_fn,
            &mut cs,
        );

        /*
        let total_delta = greedy_pass(
            &sampler_params,
            make_iter(),
            tinfo,
            model_coverage,
            density_fn,
            &mut cs,
        );
        */

        /*
        if iter > 1 && iter % 50 == 0 {
            let current_ll = cs.log_likelihood();
            if best_ll < current_ll {
                best_ll = current_ll;
                cs_backup = cs.clone();
                best_assign = cs.assigned_ids.clone();
                resample_prob /= 2.;

                resample_prob = resample_prob.clamp(0.01, 0.99);
                update_chain(
                    &sampler_params,
                    make_iter(),
                    tinfo,
                    model_coverage,
                    density_fn,
                    &mut cs,
                    resample_prob,
                );
            } else {
                cs = cs_backup.clone();
                resample_prob *= 2.;
            }
        }
        */

        let ll = cs.log_likelihood(); //_under_abundances(ml_abundances);
        //
        if iter > 1 && (ll < prev_ll) {
            warn!("prev_ll = {prev_ll} > ll = {ll}!!");
        }
        bar.set_message(format!("log_likelihood {} (at iteration {})", ll, iter));

        prev_ll = ll;

        /*
        let ll = cs.log_likelihood(); //_under_abundances(ml_abundances);
        if ll > best_log_likelihood && iter > 10 {
            best_log_likelihood = ll;
            best_log_assignments = cs.assigned_ids.clone();
            best_iter = iter;
            bar.set_message(format!(
                "best log-likelihood so far {} (at iteration {})",
                best_log_likelihood, best_iter
            ));
        }
        */
    }
    bar.finish();

    best_log_assignments = cs
        .extract_most_frequent_assignments(make_iter())
        .expect("should be a valid chain state");

    best_log_assignments
}

#[derive(Clone)]
struct ChainState {
    pub counts: Vec<u64>, // the current (integer) assignment of counts to each transcript
    pub assigned_ids: Vec<u32>, // the current assignment of each read to a transcript of origin
    pub log_cond_probs: Vec<LogSpace>, // conditional assignment probabilites for each currently assigned read (in log space)
    pub tot_counts: u64,
    pub frequencies: Option<Vec<u16>>,
}

impl ChainState {
    #[allow(dead_code)]
    fn log_likelihood_under_abundances(&self, abundances: &[f64]) -> f64 {
        let mut ll = LogSpace::new_from_linear(1f64);
        let norm = 1f64 / (self.tot_counts as f64);
        for (t_assign, c_prob) in self.assigned_ids.iter().zip(self.log_cond_probs.iter()) {
            ll *= LogSpace::new_from_linear((abundances[*t_assign as usize] + 1f64) * norm)
                * (*c_prob);
        }
        ll.get_ln()
    }

    fn log_likelihood(&self) -> f64 {
        let mut ll = LogSpace::new_from_linear(1f64);
        let norm = 1f64 / (self.tot_counts as f64);
        for (t_assign, c_prob) in self.assigned_ids.iter().zip(self.log_cond_probs.iter()) {
            ll *=
                LogSpace::new_from_linear((self.counts[*t_assign as usize] as f64) * norm + 1e-15)
                    * (*c_prob);
        }
        ll.get_ln()
    }

    fn inc_sample_counts(&mut self, aln_idx: usize) {
        if let Some(freqs) = self.frequencies.as_mut() {
            freqs[aln_idx] += 1;
        }
    }

    fn extract_most_frequent_assignments<
        'a,
        I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a,
    >(
        &mut self,
        eq_map_iter: I,
    ) -> anyhow::Result<Vec<u32>> {
        if let Some(sample_counts) = self.frequencies.as_mut() {
            let num_reads = self.tot_counts as usize;
            let mut most_frequent_assignments = vec![0u32; num_reads];
            let mut aln_idx = 0usize;
            for (read_id, (alns, probs, coverage_probs)) in eq_map_iter.enumerate() {
                // if this is a uniquely-aligned read, don't bother with
                // all of this
                if alns.len() == 1 {
                    most_frequent_assignments[read_id] = alns[0].ref_id;
                    aln_idx += 1;
                    continue;
                }

                let mut max_freq = 0_u16;
                let mut max_freq_tid = 0_u32;
                let mut idx = 0;
                for (a, _p, _cp) in izip!(alns, probs, coverage_probs) {
                    let tid_ct = sample_counts[aln_idx + idx];
                    if tid_ct > max_freq {
                        max_freq = tid_ct;
                        max_freq_tid = a.ref_id;
                    }
                    idx += 1;
                }
                most_frequent_assignments[read_id] = max_freq_tid;
                aln_idx += idx;
            }
            Ok(most_frequent_assignments)
        } else {
            anyhow::bail!(
                "Cannot extract most freuent assignment from a ChainState that doesn't track assignment frequencies"
            );
        }
    }
}

/// Obtain an initial set of counts and assignments from the maximum likelihood abundances
/// and use these to create an initial `ChainState` from which we can run the Gibbs sampler.
fn init_counts_and_assignments<
    'a,
    I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a,
    F: Fn() -> I,
>(
    em_info: &'a EMInfo,
    ml_abundances: &[f64], // abundances obtained from the EM algorithm
    make_iter: F,
    with_frequencies: bool,
) -> ChainState {
    let mut counts = vec![0; em_info.num_txps()];
    let mut assigned_ids: Vec<u32> = Vec::new();
    let mut log_cond_probs: Vec<LogSpace> = Vec::new();
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
    let mut txp_log_cond_probs = Vec::new();

    let mut tot_counts = 0;
    let mut num_aln = 0;

    for (alns, probs, coverage_probs) in make_iter() {
        txp_ids.clear();
        txp_probs.clear();
        txp_log_cond_probs.clear();

        let mut denom = 0.0_f64;
        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            num_aln += 1;
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let txp_len = tinfo[target_id].lenf as usize;
            let aln_len = a.alignment_span() as usize;

            let prob = *p as f64;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            let dens_prob = density_fn(txp_len, aln_len);

            let cond_prob = prob * cov_prob * dens_prob;
            let sprob = cond_prob * ml_abundances[target_id];
            txp_ids.push(a.ref_id);
            txp_probs.push(sprob);
            txp_log_cond_probs.push(cond_prob);
            denom += sprob;
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            let cat = WeightedIndex::new(&txp_probs).unwrap();
            let s = cat.sample(&mut rng);
            assigned_ids.push(txp_ids[s]);
            counts[txp_ids[s] as usize] += 1;
            log_cond_probs.push(LogSpace::new_from_linear(txp_log_cond_probs[s]));
            tot_counts += 1;
        }
    }

    let frequencies = if with_frequencies {
        Some(vec![0; num_aln])
    } else {
        None
    };

    ChainState {
        counts,
        assigned_ids,
        log_cond_probs,
        tot_counts,
        frequencies,
    }
}
