use std::sync::atomic::Ordering;

use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use atomic_float::AtomicF64;
use itertools::izip;
use num_format::{Locale, ToFormattedString};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use tracing::{info, span, trace};

use crate::bootstrap;

type EqIterateT<'a> = (&'a [AlnInfo], &'a [f32], &'a [f64]);

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
    prev_count: &mut [AtomicF64],
    curr_counts: &mut [AtomicF64],
) where
    DFn: Fn(usize, usize) -> f64 + Sync,
{
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

                denom +=
                    prev_count[target_id].load(Ordering::Relaxed) * prob * cov_prob * dens_prob;
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
                        * prob
                        * cov_prob
                        * dens_prob)
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
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) where
    DFn: Fn(usize, usize) -> f64,
{
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
    do_log: bool,
) -> Vec<f64> {
    let eq_map = em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let convergence_thresh = em_info.convergence_thresh;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

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

    // for up to the maximum number of iterations
    while niter < max_iter {
        // allocate the fragments and compute the new counts
        m_step(
            make_iter(),
            tinfo,
            fops.model_coverage,
            density_fn,
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
        if do_log && niter.is_multiple_of(10) {
            if niter.is_multiple_of(100) {
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

    // set very small abundances to 0
    for x in &mut prev_counts {
        if *x < constants::MIN_READ_THRESH {
            *x = 0.0;
        }
    }
    // perform one more EM round, since we just zeroed out
    // very small abundances
    m_step(
        make_iter(),
        tinfo,
        fops.model_coverage,
        density_fn,
        &mut prev_counts,
        &mut curr_counts,
    );
    //  return the final estimated abundances
    curr_counts
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
#[allow(dead_code)]
pub fn em(em_info: &EMInfo, _nthreads: usize) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    // function that produces an iterator over the
    // scored alignments on demand
    let make_iter = || em_info.eq_map.iter();

    do_em(em_info, make_iter, true)
}

pub fn do_bootstrap(em_info: &EMInfo, seed: Option<u64>, rep_idx: u64) -> Vec<f64> {
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    let mut rng: StdRng = match seed {
        Some(s) => {
            // Derive a deterministic per-replicate seed so that parallel execution
            // yields the same results regardless of scheduling.
            let derived = s.wrapping_add(0x9E37_79B9_7F4A_7C15u64.wrapping_mul(rep_idx + 1));
            StdRng::seed_from_u64(derived)
        }
        None => StdRng::from_os_rng(),
    };

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

    do_em(em_info, make_iter, false)
}

pub fn bootstrap(
    em_info: &EMInfo,
    num_boot: u32,
    nthreads: usize,
    seed: Option<u64>,
) -> Vec<Vec<f64>> {
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
                do_bootstrap(em_info, seed, i as u64)
            })
            .collect()
    })
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
pub fn em_par(em_info: &EMInfo, nthreads: usize) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    let eq_map = em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let convergence_thresh = em_info.convergence_thresh;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;
    let eq_iterates: Vec<EqIterateT> = eq_map.iter().collect();
    // initialize the estimated counts for the EM procedure
    let prev_counts: Vec<f64>;
    let mut curr_counts: Vec<AtomicF64> = vec![0.0f64; tinfo.len()]
        .iter()
        .map(|x| AtomicF64::new(*x))
        .collect();

    if let Some(ref init_counts) = em_info.init_abundances {
        // initalize with the short-read quantification
        prev_counts = init_counts.clone();
    } else {
        // uniform, length normalized abundance
        let avg = total_weight / (tinfo.len() as f64);
        prev_counts = vec![avg; tinfo.len()];
    }
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

    pool.install(|| {
        // for up to the maximum number of iterations
        while niter < max_iter {
            // allocate the fragments and compute the new counts
            m_step_par(
                &eq_iterates,
                tinfo,
                fops.model_coverage,
                density_fn,
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
            if niter.is_multiple_of(10) {
                if niter.is_multiple_of(100) {
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

        // set very small abundances to 0
        prev_counts.iter_mut().for_each(|x| {
            if x.load(Ordering::Relaxed) < constants::MIN_READ_THRESH {
                x.store(0.0, Ordering::Relaxed);
            }
        });

        // perform one more EM round, since we just zeroed out
        // very small abundances
        m_step_par(
            &eq_iterates,
            tinfo,
            fops.model_coverage,
            density_fn,
            &mut prev_counts,
            &mut curr_counts,
        );
    });
    //  return the final estimated abundances
    curr_counts
        .iter()
        .map(|x| x.load(Ordering::Relaxed))
        .collect::<Vec<f64>>()
}
