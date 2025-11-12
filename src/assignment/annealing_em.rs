use std::sync::atomic::Ordering;

use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use atomic_float::AtomicF64;
use itertools::izip;
use num_format::{Locale, ToFormattedString};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use tracing::{info, span, trace};

type EqIterateT<'a> = (&'a [AlnInfo], &'a [f32], &'a [f64]);

#[allow(dead_code)]
/// Configuration for temperature schedule
#[derive(Debug, Clone)]
pub struct AnnealingConfig {
    pub n_iterations: usize,
    pub schedule_type: ScheduleType,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum ScheduleType {
    /// Two-phase: constant T=1.0, then exponential decay
    TwoPhase {
        warm_up_frac: f64,
        decay_rate: f64, // Multiplier for exponential (default: 4.0)
    },
    /// Exponential decay throughout
    Exponential { t_start: f64, t_end: f64 },
    /// Linear decay
    Linear { t_start: f64, t_end: f64 },
}

#[allow(dead_code)]
impl AnnealingConfig {
    pub fn two_phase(n_iterations: usize, warm_up_frac: f64) -> Self {
        Self {
            n_iterations,
            schedule_type: ScheduleType::TwoPhase {
                warm_up_frac,
                decay_rate: 4.0,
            },
        }
    }

    pub fn generate_schedule(&self) -> Vec<f64> {
        match &self.schedule_type {
            ScheduleType::TwoPhase {
                warm_up_frac,
                decay_rate,
            } => {
                let warm_up_iters = (self.n_iterations as f64 * warm_up_frac).round() as usize;
                let hardening_iters = self.n_iterations - warm_up_iters;

                let mut schedule = Vec::with_capacity(self.n_iterations);
                schedule.extend(std::iter::repeat_n(1.0, warm_up_iters));

                if hardening_iters > 0 {
                    for i in 0..hardening_iters {
                        let t = i as f64 / hardening_iters as f64;
                        schedule.push((-decay_rate * t).exp());
                    }
                }

                schedule
            }
            ScheduleType::Exponential { t_start, t_end } => {
                temperature_schedule_exponential(self.n_iterations, *t_start, *t_end)
            }
            ScheduleType::Linear { t_start, t_end } => {
                temperature_schedule_linear(self.n_iterations, *t_start, *t_end)
            }
        }
    }
}

/// Generates a two-phase temperature schedule for deterministic annealing.
///
/// # Arguments
/// * `n_iterations` - Total number of EM iterations
/// * `warm_up_frac` - Fraction of iterations to spend in standard EM (phase 1)
///
/// # Returns
/// A vector of temperature values, one per iteration
///
/// # Example
/// ```
/// let schedule = temperature_schedule_two_phase(100, 0.7);
/// // First 70 iterations: T = 1.0 (standard EM)
/// // Last 30 iterations: T decays exponentially toward 0
/// ```
#[allow(dead_code)]
pub fn temperature_schedule_two_phase(n_iterations: usize, warm_up_frac: f64) -> Vec<f64> {
    assert!(
        (0.0..=1.0).contains(&warm_up_frac),
        "warm_up_frac must be between 0 and 1"
    );
    assert!(n_iterations > 0, "n_iterations must be positive");

    let warm_up_iters = (n_iterations as f64 * warm_up_frac).round() as usize;
    let hardening_iters = n_iterations - warm_up_iters;

    let mut schedule = Vec::with_capacity(n_iterations);

    // Phase 1: Standard EM with T = 1.0
    schedule.extend(std::iter::repeat_n(1.0, warm_up_iters));

    // Phase 2: Exponential decay
    // T = exp(-4 * t) where t goes from 0 to 1
    if hardening_iters > 0 {
        for i in 0..hardening_iters {
            let t = i as f64 / hardening_iters as f64;
            let temperature = (-4.0 * t).exp();
            schedule.push(temperature);
        }
    }

    schedule
}

/// Alternative: Exponential decay schedule throughout all iterations
pub fn temperature_schedule_exponential(n_iterations: usize, t_start: f64, t_end: f64) -> Vec<f64> {
    assert!(n_iterations > 0, "n_iterations must be positive");
    assert!(
        t_start > 0.0 && t_end > 0.0,
        "temperatures must be positive"
    );

    let ratio = t_end / t_start;
    (0..n_iterations)
        .map(|i| {
            let progress = i as f64 / (n_iterations - 1) as f64;
            t_start * ratio.powf(progress)
        })
        .collect()
}

/// Linear decay schedule (simple but less commonly used)
pub fn temperature_schedule_linear(n_iterations: usize, t_start: f64, t_end: f64) -> Vec<f64> {
    assert!(n_iterations > 0, "n_iterations must be positive");

    (0..n_iterations)
        .map(|i| {
            let progress = i as f64 / (n_iterations - 1) as f64;
            t_start + (t_end - t_start) * progress
        })
        .collect()
}

/// Performs one iteration of the EM algorithm by looping over all
/// alignments and computing their estimated probability of being
/// the true alignment (using the abunance estimates from `prev_counts`).
/// Then, `curr_counts` is computed by summing over the expected assignment
/// likelihood for all reads mapping to each target.
#[inline]
fn m_step_par<DFn>(
    eq_iterates: &[EqIterateT],
    tinfo: &[TranscriptInfo],
    model_coverage: bool,
    density_fn: DFn,
    temp: f64,
    prev_count: &mut [AtomicF64],
    curr_counts: &mut [AtomicF64],
) where
    DFn: Fn(usize, usize) -> f64 + Sync,
{
    let inv_temp = 1.0 / temp;
    // for (alns, probs, coverage_probs) in eq_map.iter() {
    eq_iterates.par_iter().for_each_with(
        &curr_counts,
        |curr_counts, (alns, probs, coverage_probs)| {
            let mut denom = 0.0_f64;
            for (a, p, cp) in izip!(*alns, *probs, *coverage_probs) {
                // Compute the probability of assignment of the
                // current read based on this alignment and the
                // target's estimated abundance.
                let target_id = a.ref_id as usize;

                let txp_len = tinfo[target_id].lenf as usize;
                let aln_len = a.alignment_span() as usize;

                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                let dens_prob = density_fn(txp_len, aln_len);

                denom += prev_count[target_id].load(Ordering::Relaxed)
                    * (prob * cov_prob * dens_prob).powf(inv_temp);
            }

            // If this read can be assigned
            if denom > constants::EM_DENOM_THRESH {
                // Loop over all possible assignment locations and proportionally
                // allocate the read according to our model and current parameter
                // estimates.
                for (a, p, cp) in izip!(*alns, *probs, *coverage_probs) {
                    let target_id = a.ref_id as usize;

                    let txp_len = tinfo[target_id].lenf as usize;
                    let aln_len = a.alignment_span() as usize;

                    let prob = *p as f64;
                    let cov_prob = if model_coverage { *cp } else { 1.0 };
                    let dens_prob = density_fn(txp_len, aln_len);
                    let inc = (prev_count[target_id].load(Ordering::Relaxed)
                        * (prob * cov_prob * dens_prob).powf(inv_temp))
                        / denom;
                    //curr_counts[target_id] += inc;
                    curr_counts[target_id].fetch_add(inc, Ordering::AcqRel);
                }
            }
        },
    );
}

/// Performs one iteration of the EM algorithm by looping over all
/// alignments and computing their estimated probability of being
/// the true alignment (using the abunance estimates from `prev_counts`).
/// Then, `curr_counts` is computed by summing over the expected assignment
/// likelihood for all reads mapping to each target.
#[inline]
fn m_step<'a, DFn, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])>>(
    eq_map_iter: I,
    tinfo: &[TranscriptInfo],
    model_coverage: bool,
    density_fn: DFn,
    temp: f64,
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) where
    DFn: Fn(usize, usize) -> f64,
{
    // don't let temp become too small for numerical stability
    let temp = temp.max(0.001);
    let inv_temp = 1.0 / temp;
    for (alns, probs, coverage_probs) in eq_map_iter {
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

            denom += prev_count[target_id] * (prob * cov_prob * dens_prob).powf(inv_temp);
        }

        // If this read can be assigned
        if denom > constants::EM_DENOM_THRESH {
            // Loop over all possible assignment locations and proportionally
            // allocate the read according to our model and current parameter
            // estimates.
            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;
                let txp_len = tinfo[target_id].lenf as usize;
                let aln_len = a.alignment_span() as usize;

                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                let dens_prob = density_fn(txp_len, aln_len);

                let inc =
                    (prev_count[target_id] * (prob * cov_prob * dens_prob).powf(inv_temp)) / denom;
                curr_counts[target_id] += inc;
            }
        }
    }
}

#[allow(dead_code)]
/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
pub fn solve_par(
    em_info: &EMInfo,
    schedule: AnnealingConfig,
    counts: &[f64],
    nthreads: usize,
) -> Vec<u32> {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    let eq_map = &em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let convergence_thresh = em_info.convergence_thresh;
    let eq_iterates: Vec<EqIterateT> = eq_map.iter().collect();

    // initialize the estimated counts for the EM procedure
    let mut prev_counts = Vec::with_capacity(counts.len());
    prev_counts.extend_from_slice(counts);

    let mut curr_counts: Vec<AtomicF64> = vec![0.0f64; tinfo.len()]
        .iter()
        .map(|x| AtomicF64::new(*x))
        .collect();

    let mut prev_counts: Vec<AtomicF64> = prev_counts.iter().map(|x| AtomicF64::new(*x)).collect();

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;
    let mut _fl_prob = 0.5f64;

    let density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.,
        }
    };

    let schedule_iterates = schedule.generate_schedule();

    pool.install(|| {
        // for up to the maximum number of iterations
        while niter < max_iter {
            // allocate the fragments and compute the new counts
            m_step_par(
                &eq_iterates,
                tinfo,
                fops.model_coverage,
                density_fn,
                schedule_iterates[niter as usize],
                &mut prev_counts,
                &mut curr_counts,
            );

            // compute the relative difference in the parameter estimates
            // between the current and previous rounds
            for i in 0..curr_counts.len() {
                if prev_counts[i].load(Ordering::Relaxed) > constants::MIN_READ_THRESH {
                    let cc = curr_counts[i].load(Ordering::Relaxed);
                    let pc = prev_counts[i].load(Ordering::Relaxed);
                    let rd = (cc - pc) / pc;
                    rel_diff = rel_diff.max(rd);
                }
            }

            // swap the current and previous abundances
            std::mem::swap(&mut prev_counts, &mut curr_counts);

            // clear out the new abundances
            curr_counts
                .par_iter()
                .for_each(|x| x.store(0.0f64, Ordering::Relaxed)); //fill(0.0_f64);

            // if the maximum relative difference is small enough
            // and we've done at least 10 rounds of the EM, then
            // exit (early stop).
            if (rel_diff < convergence_thresh) && (niter > 1) {
                break;
            }
            // increment the iteration and, if this iteration
            // is a multiple of 10, print out  the maximum relative
            // difference we observed.
            niter += 1;
            if niter % 10 == 0 {
                if niter % 100 == 0 {
                    info!(
                        "iteration {}; rel diff {}",
                        niter.to_formatted_string(&Locale::en),
                        rel_diff
                    );
                } else {
                    trace!(
                        "iteration {}; rel diff {}",
                        niter.to_formatted_string(&Locale::en),
                        rel_diff
                    );
                }
            }
            rel_diff = 0.0_f64;
        }
    });

    //  return the final estimated abundances
    let prev_counts = prev_counts
        .iter()
        .map(|x| x.load(Ordering::Relaxed))
        .collect::<Vec<f64>>();

    let nreads = em_info.eq_map.len();
    let mut assignments = Vec::<u32>::with_capacity(nreads);
    let inv_temp = 1.0 / (0.001f64.max(schedule_iterates[niter as usize]));
    // go through one more time and do the hard assignments
    for (alns, probs, coverage_probs) in eq_map.iter() {
        let mut highest_prob = 0.0_f64;
        let mut best_txp = 0;
        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let txp_len = tinfo[target_id].lenf as usize;
            let aln_len = a.alignment_span() as usize;

            let prob = *p as f64;
            let cov_prob = if fops.model_coverage { *cp } else { 1.0 };
            let dens_prob = density_fn(txp_len, aln_len);
            let post_prob = prev_counts[target_id] * (prob * cov_prob * dens_prob).powf(inv_temp);
            if post_prob > highest_prob {
                highest_prob = post_prob;
                best_txp = target_id;
            }
        }
        assignments.push(best_txp as u32);
    }

    //  return the final estimated abundances
    assignments
}

/// The code that actually performs the EM loop in the single-threaded context.
/// The parameters are
/// `em_info` : an [EMInfo] struct that contains the relevant parameters and data
/// `make_iter`: a function that returns an iterator over the alignments and conditional
/// probabilites
/// `do_log`: `true` if logging information is written about this EM run and false otherwise
///
/// returns:
/// A [Vec<f64>] represented the expected read assignments to each transcript.
pub fn do_em<'a, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a, F: Fn() -> I>(
    em_info: &'a EMInfo,
    make_iter: F,
    schedule: AnnealingConfig,
    do_log: bool,
) -> Vec<u32> {
    let eq_map = &em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let convergence_thresh = em_info.convergence_thresh;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;
    let model_coverage = &eq_map.filter_opts.model_coverage;

    // initialize the estimated counts for the EM procedure
    let mut prev_counts: Vec<f64>;
    let mut curr_counts: Vec<f64> = vec![0.0f64; tinfo.len()];

    if let Some(ref init_counts) = em_info.init_abundances {
        // initalize with the short-read quantification
        prev_counts = init_counts.clone();
    } else {
        // uniform, length normalized abundance
        let avg = total_weight / (tinfo.len() as f64);
        prev_counts = vec![avg; tinfo.len()];
    }

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;
    let mut _fl_prob = 0.5f64;

    let density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.,
        }
    };

    let schedule_iterates = schedule.generate_schedule();

    // for up to the maximum number of iterations
    while niter < max_iter {
        // allocate the fragments and compute the new counts
        m_step(
            make_iter(),
            tinfo,
            fops.model_coverage,
            density_fn,
            schedule_iterates[niter as usize],
            &mut prev_counts,
            &mut curr_counts,
        );

        // compute the relative difference in the parameter estimates
        // between the current and previous rounds
        for i in 0..curr_counts.len() {
            if prev_counts[i] > constants::MIN_READ_THRESH {
                let cc = curr_counts[i];
                let pc = prev_counts[i];
                let rd = (cc - pc) / pc;
                rel_diff = rel_diff.max(rd);
            }
        }

        // swap the current and previous abundances
        std::mem::swap(&mut prev_counts, &mut curr_counts);

        // clear out the new abundances
        curr_counts.fill(0.0_f64);

        // if the maximum relative difference is small enough
        // and we've done at least 10 rounds of the EM, then
        // exit (early stop).
        if (rel_diff < convergence_thresh) && (niter > 50) {
            break;
        }
        // increment the iteration and, if this iteration
        // is a multiple of 10, print out  the maximum relative
        // difference we observed.
        niter += 1;
        if do_log && (niter % 10 == 0) {
            if niter % 100 == 0 {
                info!(
                    "iteration {}; rel diff {}",
                    niter.to_formatted_string(&Locale::en),
                    rel_diff
                );
            } else {
                trace!(
                    "iteration {}; rel diff {}",
                    niter.to_formatted_string(&Locale::en),
                    rel_diff
                );
            }
        }
        rel_diff = 0.0_f64;
    }

    let nreads = em_info.eq_map.len();
    let mut assignments = Vec::<u32>::with_capacity(nreads);
    let inv_temp = 1.0 / (0.001f64.max(schedule_iterates[niter as usize]));
    // go through one more time and do the hard assignments
    for (alns, probs, coverage_probs) in make_iter() {
        let mut highest_prob = 0.0_f64;
        let mut best_txp = 0;
        for (a, p, cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let txp_len = tinfo[target_id].lenf as usize;
            let aln_len = a.alignment_span() as usize;

            let prob = *p as f64;
            let cov_prob = if *model_coverage { *cp } else { 1.0 };
            let dens_prob = density_fn(txp_len, aln_len);
            let post_prob = prev_counts[target_id] * (prob * cov_prob * dens_prob).powf(inv_temp);
            if post_prob > highest_prob {
                highest_prob = post_prob;
                best_txp = target_id;
            }
        }
        assignments.push(best_txp as u32);
    }

    //  return the final estimated abundances
    assignments
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
#[allow(dead_code)]
pub fn solve(em_info: &EMInfo, schedule: AnnealingConfig) -> Vec<u32> {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    // function that produces an iterator over the
    // scored alignments on demand
    let make_iter = || em_info.eq_map.iter();

    do_em(em_info, make_iter, schedule, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_phase_schedule() {
        let schedule = temperature_schedule_two_phase(100, 0.7);

        assert_eq!(schedule.len(), 100);

        // First 70 iterations should be 1.0
        for v in schedule.iter().take(70) {
            assert_eq!(*v, 1.0);
        }

        // Last 30 should decay
        assert!(schedule[70] < 1.0);
        assert!(schedule[99] < schedule[70]);
        assert!(schedule[99] > 0.0);
    }

    #[test]
    fn test_exponential_schedule() {
        let schedule = temperature_schedule_exponential(10, 1.0, 0.01);

        assert_eq!(schedule.len(), 10);
        assert!((schedule[0] - 1.0).abs() < 1e-10);
        assert!((schedule[9] - 0.01).abs() < 1e-10);

        // Should be monotonically decreasing
        for i in 1..schedule.len() {
            assert!(schedule[i] < schedule[i - 1]);
        }
    }

    #[test]
    fn test_linear_schedule() {
        let schedule = temperature_schedule_linear(11, 1.0, 0.0);

        assert_eq!(schedule.len(), 11);
        assert_eq!(schedule[0], 1.0);
        assert_eq!(schedule[10], 0.0);
        assert!((schedule[5] - 0.5).abs() < 1e-10);
    }
}
