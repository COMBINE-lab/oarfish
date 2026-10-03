//! Equivalence-class EM for when no per-read model term is active.
//!
//! Without the coverage model and the KDE length model, a read enters the EM only through the
//! multiset of its (target, alignment-score probability) pairs. Reads with the same multiset
//! contribute identical terms, so they are collapsed into one class weighted by its read count.
//! Probabilities are keyed by their exact bits, so the class EM computes the same estimates as
//! the per-read EM (up to floating-point summation order) while touching each distinct class once
//! per iteration.

use std::sync::atomic::Ordering;

use atomic_float::AtomicF64;
use num_format::{Locale, ToFormattedString};
use rand::rng as trng;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use rustc_hash::FxHashMap;
use tracing::{info, span, trace};

use crate::bootstrap;
use crate::util::constants;
use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};

/// Reads collapsed by their (target, probability) multisets.
pub struct EqClasses {
    /// Targets of class `c`: `targets[offsets[c]..offsets[c + 1]]`, sorted.
    targets: Vec<u32>,
    /// Alignment-score probabilities, parallel to `targets`.
    probs: Vec<f32>,
    offsets: Vec<usize>,
    /// Reads per class.
    counts: Vec<f64>,
    /// Class of each read in the store (for bootstrap resampling).
    read_class: Vec<u32>,
}

impl EqClasses {
    pub fn from_store(store: &InMemoryAlignmentStore) -> Self {
        let mut ids: FxHashMap<Vec<(u32, u32)>, u32> = FxHashMap::default();
        let mut eqc = EqClasses {
            targets: Vec::new(),
            probs: Vec::new(),
            offsets: vec![0],
            counts: Vec::new(),
            read_class: Vec::with_capacity(store.len()),
        };
        let mut key: Vec<(u32, u32)> = Vec::new();
        for (alns, probs, _) in store.iter() {
            key.clear();
            key.extend(alns.iter().zip(probs).map(|(a, p)| (a.ref_id, p.to_bits())));
            key.sort_unstable();
            let next = eqc.counts.len() as u32;
            let c = *ids.entry(key.clone()).or_insert_with(|| {
                eqc.targets.extend(key.iter().map(|&(t, _)| t));
                eqc.probs
                    .extend(key.iter().map(|&(_, p)| f32::from_bits(p)));
                eqc.offsets.push(eqc.targets.len());
                eqc.counts.push(0.0);
                next
            });
            eqc.counts[c as usize] += 1.0;
            eqc.read_class.push(c);
        }
        eqc
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    fn class(&self, c: usize) -> (&[u32], &[f32]) {
        let r = self.offsets[c]..self.offsets[c + 1];
        (&self.targets[r.clone()], &self.probs[r])
    }
}

/// One EM round: distribute each class's reads over its targets in proportion to
/// `prev[t] * prob`, accumulating into `curr`.
fn m_step(eqc: &EqClasses, counts: &[f64], prev: &[AtomicF64], curr: &[AtomicF64]) {
    (0..eqc.len()).into_par_iter().for_each(|c| {
        let n = counts[c];
        if n == 0.0 {
            return;
        }
        let (ts, ps) = eqc.class(c);
        let denom: f64 = ts
            .iter()
            .zip(ps)
            .map(|(&t, &p)| prev[t as usize].load(Ordering::Relaxed) * p as f64)
            .sum();
        if denom > constants::EM_DENOM_THRESH {
            let scale = n / denom;
            for (&t, &p) in ts.iter().zip(ps) {
                let inc = prev[t as usize].load(Ordering::Relaxed) * p as f64 * scale;
                curr[t as usize].fetch_add(inc, Ordering::AcqRel);
            }
        }
    });
}

/// The EM of [`crate::em::em_par`] over classes with read `counts` (the class counts, or a
/// bootstrap resample of them): same initialization, convergence rule and final pruning round.
fn run(
    eqc: &EqClasses,
    counts: &[f64],
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    do_log: bool,
) -> Vec<f64> {
    let total_weight: f64 = counts.iter().sum();
    let init = match init_abundances {
        Some(init) => init.clone(),
        None => vec![total_weight / (tinfo.len() as f64); tinfo.len()],
    };
    let mut prev: Vec<AtomicF64> = init.into_iter().map(AtomicF64::new).collect();
    let mut curr: Vec<AtomicF64> = (0..tinfo.len()).map(|_| AtomicF64::new(0.0)).collect();

    let mut niter = 0_u32;
    while niter < max_iter {
        m_step(eqc, counts, &prev, &curr);
        let mut rel_diff = 0.0_f64;
        for (p, c) in prev.iter().zip(&curr) {
            let pc = p.load(Ordering::Relaxed);
            if pc > constants::MIN_READ_THRESH {
                rel_diff = rel_diff.max((c.load(Ordering::Relaxed) - pc) / pc);
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        curr.par_iter()
            .for_each(|x| x.store(0.0, Ordering::Relaxed));
        if rel_diff < convergence_thresh && niter > 1 {
            break;
        }
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
    }

    // Zero very small abundances, then one more round, as the per-read EM does.
    for x in &prev {
        if x.load(Ordering::Relaxed) < constants::MIN_READ_THRESH {
            x.store(0.0, Ordering::Relaxed);
        }
    }
    m_step(eqc, counts, &prev, &curr);
    curr.iter().map(|x| x.load(Ordering::Relaxed)).collect()
}

/// Abundance estimates (expected reads per target).
pub fn em(
    eqc: &EqClasses,
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    nthreads: usize,
) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        run(
            eqc,
            &eqc.counts,
            tinfo,
            max_iter,
            convergence_thresh,
            init_abundances,
            true,
        )
    })
}

/// Bootstrap replicates: each resamples the reads with replacement (as the per-read bootstrap
/// does) and runs the EM on the resulting class counts.
pub fn bootstrap(
    eqc: &EqClasses,
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    num_boot: u32,
    nthreads: usize,
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
                info!("evaluating bootstrap replicate {}", i);
                let mut rng = trng();
                let mut counts = vec![0.0_f64; eqc.len()];
                for r in bootstrap::get_sample_inds(eqc.read_class.len(), &mut rng) {
                    counts[eqc.read_class[r] as usize] += 1.0;
                }
                run(
                    eqc,
                    &counts,
                    tinfo,
                    max_iter,
                    convergence_thresh,
                    init_abundances,
                    false,
                )
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The class EM equals a plain per-read EM on the same reads.
    #[test]
    fn class_em_matches_per_read_em() {
        // Reads as (target, prob) lists; several repeat, in different orders.
        let reads: Vec<Vec<(u32, f32)>> = vec![
            vec![(0, 1.0), (1, 0.5)],
            vec![(1, 0.5), (0, 1.0)],
            vec![(1, 1.0), (2, 1.0)],
            vec![(0, 1.0)],
            vec![(2, 0.25), (0, 1.0), (1, 1.0)],
            vec![(1, 1.0), (2, 1.0)],
            vec![(2, 1.0)],
        ];
        let mut eqc = EqClasses {
            targets: vec![],
            probs: vec![],
            offsets: vec![0],
            counts: vec![],
            read_class: vec![],
        };
        let mut ids: FxHashMap<Vec<(u32, u32)>, u32> = FxHashMap::default();
        for r in &reads {
            let mut key: Vec<(u32, u32)> = r.iter().map(|&(t, p)| (t, p.to_bits())).collect();
            key.sort_unstable();
            let next = eqc.counts.len() as u32;
            let c = *ids.entry(key.clone()).or_insert_with(|| {
                eqc.targets.extend(key.iter().map(|&(t, _)| t));
                eqc.probs
                    .extend(key.iter().map(|&(_, p)| f32::from_bits(p)));
                eqc.offsets.push(eqc.targets.len());
                eqc.counts.push(0.0);
                next
            });
            eqc.counts[c as usize] += 1.0;
            eqc.read_class.push(c);
        }
        assert_eq!(eqc.len(), 5);
        let tinfo: Vec<TranscriptInfo> = (0..3)
            .map(|_| TranscriptInfo::with_len(std::num::NonZeroUsize::new(1000).unwrap()))
            .collect();
        let got = run(&eqc, &eqc.counts, &tinfo, 200, 1e-12, None, false);

        // Per-read EM with the same initialization, rounds and final pruning round.
        let mut prev = vec![reads.len() as f64 / 3.0; 3];
        let step = |prev: &[f64]| {
            let mut curr = vec![0.0; 3];
            for r in &reads {
                let d: f64 = r.iter().map(|&(t, p)| prev[t as usize] * p as f64).sum();
                for &(t, p) in r {
                    curr[t as usize] += prev[t as usize] * p as f64 / d;
                }
            }
            curr
        };
        for niter in 0..200u32 {
            let curr = step(&prev);
            let rd = prev
                .iter()
                .zip(&curr)
                .filter(|(p, _)| **p > constants::MIN_READ_THRESH)
                .map(|(p, c)| (c - p) / p)
                .fold(0.0_f64, f64::max);
            prev = curr;
            if rd < 1e-12 && niter > 1 {
                break;
            }
        }
        for x in &mut prev {
            if *x < constants::MIN_READ_THRESH {
                *x = 0.0;
            }
        }
        let want = step(&prev);
        for (g, w) in got.iter().zip(&want) {
            assert!((g - w).abs() < 1e-9, "{got:?} vs {want:?}");
        }
        assert!((got.iter().sum::<f64>() - reads.len() as f64).abs() < 1e-9);
    }
}
