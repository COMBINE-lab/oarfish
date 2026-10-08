//! Equivalence-class EM for when no per-read model term is active.
//!
//! Without the coverage model and the KDE length model, a read enters the EM only through the
//! multiset of its (target, alignment-score probability) pairs. Reads with the same multiset
//! contribute identical terms, so they are collapsed into one class weighted by its read count.
//! Probabilities are keyed by their exact bits, so the class EM computes the same estimates as
//! the per-read EM (up to floating-point summation order) while touching each distinct class once
//! per iteration. Each iteration uses the CSR layout of `crate::em_read`: one pass over classes
//! computes each class's scale, and one pass over transcripts sums the scales through a
//! transpose, so no update is shared between threads.

use num_format::{Locale, ToFormattedString};
use rand::rng as trng;
use rayon::prelude::*;
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
    /// Transpose: transcript `t`'s (class, probability) entries are
    /// `t_class[t_off[t]..t_off[t + 1]]` and the same range of `t_prob`.
    t_off: Vec<usize>,
    t_class: Vec<u32>,
    t_prob: Vec<f32>,
}

impl EqClasses {
    pub fn from_store(store: &InMemoryAlignmentStore, ntx: usize) -> Self {
        Self::from_reads(
            store.iter().map(|(alns, probs, _)| {
                alns.iter()
                    .zip(probs)
                    .map(|(a, &p)| (a.ref_id, p))
            }),
            store.len(),
            ntx,
        )
    }

    /// Classes from each read's (target, probability) pairs, over `ntx` transcripts.
    fn from_reads<R, I>(reads: R, nreads: usize, ntx: usize) -> Self
    where
        R: Iterator<Item = I>,
        I: Iterator<Item = (u32, f32)>,
    {
        let mut ids: FxHashMap<Vec<(u32, u32)>, u32> = FxHashMap::default();
        let mut eqc = EqClasses {
            targets: Vec::new(),
            probs: Vec::new(),
            offsets: vec![0],
            counts: Vec::new(),
            read_class: Vec::with_capacity(nreads),
            t_off: Vec::new(),
            t_class: Vec::new(),
            t_prob: Vec::new(),
        };
        let mut key: Vec<(u32, u32)> = Vec::new();
        for read in reads {
            key.clear();
            key.extend(read.map(|(t, p)| (t, p.to_bits())));
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
        eqc.transpose(ntx);
        eqc
    }

    /// Fill `t_off`, `t_class` and `t_prob` by a counting sort of the entries on target.
    fn transpose(&mut self, ntx: usize) {
        let mut off = vec![0_usize; ntx + 1];
        for &t in &self.targets {
            off[t as usize + 1] += 1;
        }
        for t in 0..ntx {
            off[t + 1] += off[t];
        }
        let mut next = off[..ntx].to_vec();
        let mut class = vec![0_u32; self.targets.len()];
        let mut prob = vec![0_f32; self.targets.len()];
        for c in 0..self.len() {
            let (ts, ps) = self.class(c);
            for (&t, &p) in ts.iter().zip(ps) {
                let k = &mut next[t as usize];
                class[*k] = c as u32;
                prob[*k] = p;
                *k += 1;
            }
        }
        (self.t_off, self.t_class, self.t_prob) = (off, class, prob);
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
/// `prev[t] * prob`, writing `curr`. `scale[c]` receives class `c`'s reads over its normalizer.
fn m_step(eqc: &EqClasses, counts: &[f64], prev: &[f64], scale: &mut [f64], curr: &mut [f64]) {
    scale.par_iter_mut().enumerate().for_each(|(c, s)| {
        let n = counts[c];
        *s = 0.0;
        if n == 0.0 {
            return;
        }
        let (ts, ps) = eqc.class(c);
        let denom: f64 = ts
            .iter()
            .zip(ps)
            .map(|(&t, &p)| prev[t as usize] * p as f64)
            .sum();
        if denom > constants::EM_DENOM_THRESH {
            *s = n / denom;
        }
    });
    curr.par_iter_mut()
        .with_min_len(1024)
        .enumerate()
        .for_each(|(t, x)| {
            let r = eqc.t_off[t]..eqc.t_off[t + 1];
            let sum: f64 = eqc.t_class[r.clone()]
                .iter()
                .zip(&eqc.t_prob[r])
                .map(|(&c, &p)| p as f64 * scale[c as usize])
                .sum();
            *x = prev[t] * sum;
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
    let mut prev = init;
    let mut curr = vec![0.0_f64; tinfo.len()];
    let mut scale = vec![0.0_f64; eqc.len()];

    let mut niter = 0_u32;
    while niter < max_iter {
        m_step(eqc, counts, &prev, &mut scale, &mut curr);
        let mut rel_diff = 0.0_f64;
        for (&pc, &c) in prev.iter().zip(&curr) {
            if pc > constants::MIN_READ_THRESH {
                rel_diff = rel_diff.max((c - pc) / pc);
            }
        }
        std::mem::swap(&mut prev, &mut curr);
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
    for x in &mut prev {
        if *x < constants::MIN_READ_THRESH {
            *x = 0.0;
        }
    }
    m_step(eqc, counts, &prev, &mut scale, &mut curr);
    curr
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
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), 3);
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
