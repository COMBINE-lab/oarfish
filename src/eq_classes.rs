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
//!
//! Without the presence model (which couples all transcripts through its base rate), the same
//! rounds run on the component layout of [`crate::em_components`] instead.

use num_format::{Locale, ToFormattedString};
use rand::rng as trng;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use tracing::{info, span, trace};

use crate::bootstrap;
use crate::em_components::{self, Counts, Plan};
use crate::presence::{Presence, PresenceOpts};
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
            store
                .iter()
                .map(|(alns, probs, _)| alns.iter().zip(probs).map(|(a, &p)| (a.ref_id, p))),
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

    /// The component layout of the classes (see [`crate::em_components`]).
    fn plan(&self, ntx: usize) -> Plan<f32> {
        Plan::new(
            |k| self.targets[k],
            |k| self.probs[k],
            &self.offsets,
            ntx,
            false,
        )
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

/// The initial abundances of [`run`]: `init_abundances`, or the reads spread evenly.
fn initial(counts: &[f64], ntx: usize, init_abundances: Option<&Vec<f64>>) -> Vec<f64> {
    match init_abundances {
        Some(init) => init.clone(),
        None => vec![counts.iter().sum::<f64>() / ntx as f64; ntx],
    }
}

/// The EM of [`crate::em::em_par`] over classes with read `counts` (the class counts, or a
/// bootstrap resample of them): same initialization, convergence rule and final pruning round.
#[allow(clippy::too_many_arguments)]
fn run(
    eqc: &EqClasses,
    counts: &[f64],
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    presence: Option<PresenceOpts>,
    do_log: bool,
) -> Vec<f64> {
    let mut prev = initial(counts, tinfo.len(), init_abundances);
    let mut curr = vec![0.0_f64; tinfo.len()];
    let mut scale = vec![0.0_f64; eqc.len()];
    // With the presence model the EM sees `q * prev` (see `crate::presence`).
    let mut pres = presence.map(|o| Presence::new(tinfo.len(), o));
    let mut eff = if pres.is_some() {
        vec![0.0_f64; tinfo.len()]
    } else {
        Vec::new()
    };

    let mut niter = 0_u32;
    while niter < max_iter {
        match &mut pres {
            Some(p) => {
                if p.due(niter) {
                    let reads = crate::presence::ReadAlns {
                        off: &eqc.offsets,
                        target: &|k| eqc.targets[k],
                    };
                    p.refresh(&prev, reads, |k| eqc.probs[k] as f64, |c| counts[c]);
                }
                p.effective(&prev, &mut eff);
                m_step(eqc, counts, &eff, &mut scale, &mut curr);
            }
            None => m_step(eqc, counts, &prev, &mut scale, &mut curr),
        }
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
    match &pres {
        Some(p) => {
            if do_log {
                info!(
                    "presence model: {} transcripts called absent (q < 0.5){}",
                    p.n_absent().to_formatted_string(&Locale::en),
                    p.base_rate()
                        .map_or(String::new(), |pi| format!("; base rate {pi:.3}"))
                );
            }
            p.effective(&prev, &mut eff);
            m_step(eqc, counts, &eff, &mut scale, &mut curr);
        }
        None => m_step(eqc, counts, &prev, &mut scale, &mut curr),
    }
    curr
}

/// Abundance estimates (expected reads per target).
pub fn em(
    eqc: &EqClasses,
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    presence: Option<PresenceOpts>,
    nthreads: usize,
) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        if presence.is_none() && em_components::enabled() {
            let t0 = std::time::Instant::now();
            let plan = eqc.plan(tinfo.len());
            info!(
                "{} connected components, laid out in {:.2}s",
                plan.n_components().to_formatted_string(&Locale::en),
                t0.elapsed().as_secs_f64()
            );
            let init = initial(&eqc.counts, tinfo.len(), init_abundances);
            return em_components::run(
                &plan,
                Counts::Of(&eqc.counts),
                &init,
                max_iter,
                convergence_thresh,
                1,
            );
        }
        run(
            eqc,
            &eqc.counts,
            tinfo,
            max_iter,
            convergence_thresh,
            init_abundances,
            presence,
            true,
        )
    })
}

/// Bootstrap replicates: each resamples the reads with replacement (as the per-read bootstrap
/// does) and runs the EM on the resulting class counts.
#[allow(clippy::too_many_arguments)]
pub fn bootstrap(
    eqc: &EqClasses,
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    presence: Option<PresenceOpts>,
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
        // One component layout serves every replicate; only the class counts differ.
        let plan = (presence.is_none() && em_components::enabled()).then(|| eqc.plan(tinfo.len()));
        (0..num_boot)
            .into_par_iter()
            .map(|i| {
                info!("evaluating bootstrap replicate {}", i);
                let mut rng = trng();
                let mut counts = vec![0.0_f64; eqc.len()];
                for r in bootstrap::get_sample_inds(eqc.read_class.len(), &mut rng) {
                    counts[eqc.read_class[r] as usize] += 1.0;
                }
                if let Some(plan) = &plan {
                    let init = initial(&counts, tinfo.len(), init_abundances);
                    return em_components::run(
                        plan,
                        Counts::Of(&counts),
                        &init,
                        max_iter,
                        convergence_thresh,
                        1,
                    );
                }
                run(
                    eqc,
                    &counts,
                    tinfo,
                    max_iter,
                    convergence_thresh,
                    init_abundances,
                    presence,
                    false,
                )
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reads over `ntx` transcripts in groups of `group` (each read's targets drawn from one
    /// group, so the groups are separate components), with some single-target reads, from a
    /// fixed pseudo-random sequence. The last transcript is never a target.
    fn grouped_reads(nreads: usize, ntx: usize, group: usize) -> Vec<Vec<(u32, f32)>> {
        let mut x: u64 = 0x9e3779b97f4a7c15;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let ngroups = (ntx - 1) / group;
        (0..nreads)
            .map(|_| {
                let g = (next() as usize) % ngroups;
                let k = 1 + (next() as usize) % 3;
                let mut ts: Vec<u32> = (0..k)
                    .map(|_| (g * group + (next() as usize) % group) as u32)
                    .collect();
                ts.sort_unstable();
                ts.dedup();
                ts.into_iter()
                    .map(|t| (t, [1.0_f32, 0.5, 0.25][(next() as usize) % 3]))
                    .collect()
            })
            .collect()
    }

    fn tinfo(ntx: usize) -> Vec<TranscriptInfo> {
        (0..ntx)
            .map(|_| TranscriptInfo::with_len(std::num::NonZeroUsize::new(1000).unwrap()))
            .collect()
    }

    /// The component layout runs the class loop's rounds: the same estimates, up to summation
    /// order, whether the global criterion stops it inside a tile or the iteration cap does.
    #[test]
    fn component_em_matches_class_em() {
        let (ntx, reads) = (301, grouped_reads(5000, 301, 6));
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), ntx);
        let ti = tinfo(ntx);
        let plan = eqc.plan(ntx);
        assert!(plan.n_components() > 10);
        for (max_iter, thresh) in [(1000, 1e-3), (1000, 1e-9), (150, 1e-12)] {
            let want = run(&eqc, &eqc.counts, &ti, max_iter, thresh, None, None, false);
            let init = initial(&eqc.counts, ntx, None);
            let got =
                em_components::run(&plan, Counts::Of(&eqc.counts), &init, max_iter, thresh, 1);
            for (t, (g, w)) in got.iter().zip(&want).enumerate() {
                assert!(
                    (g - w).abs() <= 1e-9 * w.max(1.0),
                    "{max_iter} {thresh} t{t}: {g} vs {w}"
                );
            }
            assert_eq!(got[ntx - 1], 0.0);
        }
    }

    /// The component EM, a large component's parallel rounds included, gives the same bits on
    /// any number of threads.
    #[test]
    fn component_em_is_independent_of_threads() {
        // Two groups of 150 transcripts: components large enough to step in parallel.
        let (ntx, reads) = (301, grouped_reads(200_000, 301, 150));
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), ntx);
        let plan = eqc.plan(ntx);
        assert!(plan.has_large());
        let init = initial(&eqc.counts, ntx, None);
        let on = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| em_components::run(&plan, Counts::Of(&eqc.counts), &init, 300, 1e-6, 1))
        };
        let one = on(1);
        for threads in [3, 8] {
            assert!(
                one.iter()
                    .zip(on(threads))
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
        }
    }

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
        let got = run(&eqc, &eqc.counts, &tinfo, 200, 1e-12, None, None, false);

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
