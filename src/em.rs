use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use itertools::izip;
use num_format::{Locale, ToFormattedString};
use rand::rng as trng;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use rayon::prelude::IntoParallelRefMutIterator;
use tracing::{info, span, trace};

use crate::bootstrap;

struct PreparedEq<'a> {
    alignments: &'a [AlnInfo],
    weight_start: usize,
}

fn prepare_equivalence_classes<'a>(
    em_info: &'a EMInfo<'_, '_, '_>,
) -> (Vec<PreparedEq<'a>>, Vec<f64>) {
    let model_coverage = em_info.eq_map.filter_opts.model_coverage;
    let density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.0,
        }
    };
    let mut weights = Vec::with_capacity(em_info.eq_map.alignments.len());
    let mut classes = Vec::with_capacity(em_info.eq_map.len());
    for (alns, probs, coverage_probs) in em_info.eq_map.iter() {
        let weight_start = weights.len();
        weights.extend(izip!(alns, probs, coverage_probs).map(|(a, p, cp)| {
            let target_id = a.ref_id as usize;
            let txp_len = em_info.txp_info[target_id].lenf as usize;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            *p as f64 * cov_prob * density_fn(txp_len, a.alignment_span() as usize)
        }));
        classes.push(PreparedEq {
            alignments: alns,
            weight_start,
        });
    }
    (classes, weights)
}

#[inline]
fn m_step_prepared(
    eq_iterates: &[PreparedEq],
    weights: &[f64],
    prev_count: &[f64],
    curr_counts: &mut [f64],
) {
    for eq in eq_iterates {
        let alns = eq.alignments;
        if let [alignment] = alns {
            curr_counts[alignment.ref_id as usize] += 1.0;
            continue;
        }
        let eq_weights = &weights[eq.weight_start..eq.weight_start + alns.len()];
        let denom: f64 = alns
            .iter()
            .zip(eq_weights)
            .map(|(a, weight)| prev_count[a.ref_id as usize] * weight)
            .sum();
        if denom > constants::EM_DENOM_THRESH {
            for (a, weight) in alns.iter().zip(eq_weights) {
                let target_id = a.ref_id as usize;
                curr_counts[target_id] += (prev_count[target_id] * weight) / denom;
            }
        }
    }
}

/// Performs one iteration of the EM algorithm by looping over all
/// alignments and computing their estimated probability of being
/// the true alignment (using the abunance estimates from `prev_counts`).
/// Then, `curr_counts` is computed by summing over the expected assignment
/// likelihood for all reads mapping to each target.
#[inline]
fn m_step_par(
    eq_iterates: &[PreparedEq],
    weights: &[f64],
    prev_count: &[f64],
    curr_counts: &mut [f64],
    shards: &mut [Vec<f64>],
) {
    let chunk_size = eq_iterates.len().div_ceil(shards.len().max(1));
    shards
        .par_iter_mut()
        .enumerate()
        .for_each(|(shard, local)| {
            local.fill(0.0);
            let start = shard * chunk_size;
            let end = ((shard + 1) * chunk_size).min(eq_iterates.len());
            for eq in &eq_iterates[start..end] {
                let alns = eq.alignments;
                // A read with one candidate contributes exactly one count: its
                // alignment/coverage weight cancels between numerator and
                // denominator. Avoid evaluating that weight twice.
                if let [alignment] = alns {
                    local[alignment.ref_id as usize] += 1.0;
                    continue;
                }
                let eq_weights = &weights[eq.weight_start..eq.weight_start + alns.len()];
                let mut denom = 0.0_f64;
                for (a, weight) in alns.iter().zip(eq_weights) {
                    let target_id = a.ref_id as usize;
                    denom += prev_count[target_id] * weight;
                }
                if denom > constants::EM_DENOM_THRESH {
                    for (a, weight) in alns.iter().zip(eq_weights) {
                        let target_id = a.ref_id as usize;
                        local[target_id] += (prev_count[target_id] * weight) / denom;
                    }
                }
            }
        });
    curr_counts
        .par_iter_mut()
        .enumerate()
        .for_each(|(target_id, count)| {
            *count = shards.iter().map(|local| local[target_id]).sum();
        });
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
    prev_count: &[f64],
    curr_counts: &mut [f64],
) where
    DFn: Fn(usize, usize) -> f64,
{
    for (alns, probs, coverage_probs) in eq_map_iter {
        if let [alignment] = alns {
            curr_counts[alignment.ref_id as usize] += 1.0;
            continue;
        }
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

            denom += prev_count[target_id] * prob * cov_prob * dens_prob;
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

                let inc = (prev_count[target_id] * prob * cov_prob * dens_prob) / denom;
                curr_counts[target_id] += inc;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct EMResult {
    pub counts: Vec<f64>,
    pub evaluations: u32,
    pub converged: bool,
}

fn convergence_distance(previous: &[f64], current: &[f64]) -> f64 {
    previous
        .iter()
        .zip(current)
        .filter(|(p, _)| **p > constants::MIN_READ_THRESH)
        .map(|(p, c)| (c - p).abs() / c.abs().max(constants::EM_DENOM_THRESH))
        .fold(0.0, f64::max)
}

fn run_driver(
    em_info: &EMInfo,
    mut fixed_point: impl FnMut(&[f64], &mut [f64]),
    do_log: bool,
) -> EMResult {
    const MIN_EVAL: u32 = 50;
    let n = em_info.txp_info.len();
    let avg = if n == 0 {
        0.0
    } else {
        em_info.eq_map.num_aligned_reads() as f64 / n as f64
    };
    let has_initial_abundances = em_info
        .init_abundances
        .as_ref()
        .is_some_and(|v| v.len() == n);
    let mut counts = em_info
        .init_abundances
        .as_ref()
        .filter(|v| v.len() == n)
        .cloned()
        .unwrap_or_else(|| vec![avg; n]);
    // Exact zeros are absorbing states for multiplicative EM updates.  A
    // warm-start estimate may contain zeros even though the new likelihood
    // supports those transcripts, so retain a tiny positive route back into
    // the active set.  This also makes accelerator choice independent of the
    // warm-start's zero pattern.
    if has_initial_abundances {
        counts
            .iter_mut()
            .for_each(|x| *x = x.max(10.0 * constants::MIN_READ_THRESH));
    }
    let mut next = vec![0.0; n];
    let mut evaluations = 0;
    // Reserve one evaluation for the historical post-threshold redistribution
    // step so --max-em-iter is a strict total M-step budget.
    let iteration_budget = em_info.max_iter.saturating_sub(1);

    let converged = match em_info.accel {
        crate::prog_opts::EmAccel::None => {
            let mut done = false;
            while evaluations < iteration_budget {
                fixed_point(&counts, &mut next);
                evaluations += 1;
                let d = convergence_distance(&counts, &next);
                std::mem::swap(&mut counts, &mut next);
                if do_log && evaluations.is_multiple_of(10) {
                    if evaluations.is_multiple_of(100) {
                        info!(
                            "iteration {}; rel diff {}",
                            evaluations.to_formatted_string(&Locale::en),
                            d
                        );
                    } else {
                        trace!(
                            "iteration {}; rel diff {}",
                            evaluations.to_formatted_string(&Locale::en),
                            d
                        );
                    }
                }
                if evaluations >= MIN_EVAL.min(iteration_budget)
                    && d.is_finite()
                    && d < em_info.convergence_thresh
                {
                    done = true;
                    break;
                }
            }
            done
        }
        crate::prog_opts::EmAccel::Squarem => crate::em_accel::squarem(
            &mut fixed_point,
            &mut counts,
            &mut next,
            iteration_budget,
            MIN_EVAL.min(iteration_budget),
            &mut evaluations,
            em_info.convergence_thresh,
            convergence_distance,
        ),
        crate::prog_opts::EmAccel::Daarem => crate::em_accel::daarem(
            &mut fixed_point,
            &mut counts,
            &mut next,
            iteration_budget,
            MIN_EVAL.min(iteration_budget),
            &mut evaluations,
            em_info.convergence_thresh,
            convergence_distance,
        ),
    };

    counts.iter_mut().for_each(|x| {
        if *x < constants::MIN_READ_THRESH {
            *x = 0.0;
        }
    });
    if evaluations < em_info.max_iter {
        fixed_point(&counts, &mut next);
        evaluations += 1;
    } else {
        next.copy_from_slice(&counts);
    }
    if do_log {
        info!(evaluations, converged, accelerator = ?em_info.accel, "EM completed");
    }
    EMResult {
        counts: next,
        evaluations,
        converged,
    }
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
#[allow(dead_code)]
pub fn em(em_info: &EMInfo, _nthreads: usize) -> EMResult {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    let (eq_iterates, weights) = prepare_equivalence_classes(em_info);
    let mut fixed_point = |src: &[f64], dst: &mut [f64]| {
        dst.fill(0.0);
        m_step_prepared(&eq_iterates, &weights, src, dst);
    };
    run_driver(em_info, &mut fixed_point, true)
}

pub fn do_bootstrap(em_info: &EMInfo) -> Vec<f64> {
    let mut rng = trng();
    let n = em_info.eq_map.len();
    let inds = bootstrap::get_sample_inds(n, &mut rng);

    // to not sample the indices but instead just
    // run with all reads sampled once
    /*
    let mut inds = Vec::with_capacity(n);
    for i in 0..n { inds.push(i); }
    */

    // function that produces an iterator over the
    // scored alignments on demand
    let make_iter = || em_info.eq_map.random_sampling_iter(&inds);

    let density_fn = |x, y| match em_info.kde_model {
        Some(ref kde) => kde[(x, y)],
        None => 1.0,
    };
    let mut fixed_point = |src: &[f64], dst: &mut [f64]| {
        dst.fill(0.0);
        m_step(
            make_iter(),
            em_info.txp_info,
            em_info.eq_map.filter_opts.model_coverage,
            density_fn,
            src,
            dst,
        );
    };
    run_driver(em_info, &mut fixed_point, false).counts
}

pub fn bootstrap(em_info: &EMInfo, num_boot: u32, nthreads: usize) -> Vec<Vec<f64>> {
    let span = span!(tracing::Level::INFO, "bootstrap");
    let _guard = span.enter();

    info!("will collection {num_boot} bootstraps");

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    pool.install(|| {
        (0..num_boot)
            .into_par_iter()
            .map(|i| {
                let span = span!(tracing::Level::INFO, "bootstrap");
                let _guard = span.enter();
                info!("evaluating bootstrap replicate {}", i);
                do_bootstrap(em_info)
            })
            .collect()
    })
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
pub fn em_par(em_info: &EMInfo, nthreads: usize) -> EMResult {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    // Pack invariant alignment, coverage and density terms once. They are
    // otherwise recomputed twice per alignment on every EM iteration.
    let (eq_iterates, weights) = prepare_equivalence_classes(em_info);
    let mut curr_counts = vec![0.0; em_info.txp_info.len()];
    // Reuse private dense accumulators across all fixed-point evaluations.
    // Capping the shard count avoids excessive memory use on high-core hosts.
    let shard_count = nthreads.clamp(1, 64);
    let mut shards = vec![vec![0.0; em_info.txp_info.len()]; shard_count];
    let mut fixed_point = |src: &[f64], dst: &mut [f64]| {
        m_step_par(&eq_iterates, &weights, src, &mut curr_counts, &mut shards);
        dst.copy_from_slice(&curr_counts);
    };
    pool.install(|| run_driver(em_info, &mut fixed_point, true))
}
