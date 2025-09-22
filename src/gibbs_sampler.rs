use core::f64;

use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use crate::util::probs::LogSpace;
use itertools::*;
use rand::distr::Distribution;
use rand_distr::weighted::WeightedIndex;
use std::default::Default;
use tracing::info;

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
    fn new_bitseq_with_iter(niter: u64) -> Self {
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
    let mut txp_cond_probs = Vec::new();

    // -1 because we will use this to normalize each step in the chain
    // where we have removed 1 read
    let tot_reads = (chain_state.assigned_ids.len() - 1) as f64;
    // total number of transcripts in the transcriptome (not just expressed)
    let num_txps = chain_state.counts.len() as f64;

    for (read_id, (alns, probs, coverage_probs)) in eq_map_iter.enumerate() {
        // if this is a uniquely-aligned read, don't bother with
        // all of this
        if alns.len() == 1 {
            continue;
        }
        // clear out our temporary storage vectors
        txp_ids.clear();
        txp_probs.clear();
        txp_cond_probs.clear();

        // the ID assigned to this transcript before the reassignment
        let assigned_tid: usize = chain_state.assigned_ids[read_id] as usize;
        // remove this read from this transcript's counts
        chain_state.counts[assigned_tid] -= 1;
        // the normalizing constant
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
            let cond_prob = prob * cov_prob * dens_prob;
            let sprob = cond_prob
                * (a_act + tot_reads)
                * ((a_dir + chain_state.counts[a.ref_id as usize] as f64)
                    / (num_txps * a_dir + tot_reads));

            txp_cond_probs.push(LogSpace::new_from_linear(cond_prob));
            txp_ids.push(a.ref_id);
            txp_probs.push(sprob);
            denom += sprob;
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            let cat = WeightedIndex::new(&txp_probs).unwrap();
            let s = cat.sample(&mut rng);
            let new_assigned_tid = txp_ids[s] as usize;
            // set the assignment for this read
            chain_state.assigned_ids[read_id] = new_assigned_tid as u32;
            // update the count accordingly
            chain_state.counts[new_assigned_tid] += 1;
            // update the conditional probability
            chain_state.cond_probs[read_id] = txp_cond_probs[s];
        } else {
            // the read doesn't move
            chain_state.counts[assigned_tid] += 1;
            // conditional probability remains the same
        }
    }
}

/// run the sampler and return the highest log-likelihood assignments
pub fn run_sampler<
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

    for _ in 0..sampler_params.niter {
        //info!("posterior Gibbs sampler iteration {i}");
        bar.inc(1);
        update_chain(
            &sampler_params,
            make_iter(),
            tinfo,
            model_coverage,
            density_fn,
            &mut cs,
        );
        let ll = cs.log_likelihood();
        if ll > best_log_likelihood {
            best_log_likelihood = ll;
            best_log_assignments = cs.assigned_ids.clone();
            //info!("new best log-likelihood = {}", best_log_likelihood);
            bar.set_message(format!(
                "best log-likelihood so far {}",
                best_log_likelihood
            ));
        }
    }
    bar.finish();

    best_log_assignments
}

struct ChainState {
    pub counts: Vec<u64>, // the current (integer) assignment of counts to each transcript
    pub assigned_ids: Vec<u32>, // the current assignment of each read to a transcript of origin
    pub cond_probs: Vec<LogSpace>, // conditional assignment probabilites for each currently assigned read (in log space)
    pub tot_counts: u64,
}

impl ChainState {
    fn log_likelihood(&self) -> f64 {
        let mut ll = LogSpace::new_from_linear(1f64);
        let norm = 1f64 / (self.tot_counts as f64);
        for (t_assign, c_prob) in self.assigned_ids.iter().zip(self.cond_probs.iter()) {
            ll *= LogSpace::new_from_linear(self.counts[*t_assign as usize] as f64 * norm)
                * (*c_prob);
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
    let mut cond_probs: Vec<LogSpace> = Vec::new();
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
    let mut txp_cond_probs = Vec::new();

    let mut tot_counts = 0;

    for (alns, probs, coverage_probs) in make_iter() {
        txp_ids.clear();
        txp_probs.clear();
        txp_cond_probs.clear();

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

            let cond_prob = prob * cov_prob * dens_prob;
            let sprob = ml_abundances[target_id];
            txp_ids.push(a.ref_id);
            txp_probs.push(sprob);
            txp_cond_probs.push(cond_prob);
            denom += sprob;
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            let cat = WeightedIndex::new(&txp_probs).unwrap();
            let s = cat.sample(&mut rng);
            assigned_ids.push(txp_ids[s]);
            counts[txp_ids[s] as usize] += 1;
            cond_probs.push(LogSpace::new_from_linear(txp_cond_probs[s]));
            tot_counts += 1;
        }
    }

    ChainState {
        counts,
        assigned_ids,
        cond_probs,
        tot_counts,
    }
}
