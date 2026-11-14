//! Per-read EM without atomics, used in place of `em::em_par` (`Csr` by default).
//!
//! In the per-read EM (with or without the coverage model) each alignment's weight, its score
//! probability times its coverage and length-density terms, is fixed across iterations, so an
//! iteration only redistributes each read over its alignments in proportion to
//! `prev[t] * w`. `em::em_par` adds every share into the shared abundance vector with an
//! atomic compare-and-swap. Two layouts avoid that:
//!
//! - `Local`: reads are split into chunks; each chunk accumulates into its own abundance
//!   vector, and the vectors are summed in parallel over blocks of transcripts.
//! - `Csr`: pass 1, over reads, computes each read's normalizer; pass 2, over transcripts,
//!   sums `w / normalizer` over the transcript's alignments through a precomputed transpose.
//!   Each pass writes only its own slots, and the result does not depend on scheduling.
//!
//! `OARFISH_READ_EM=local` or `OARFISH_READ_EM=atomic` selects the other layouts. With `Csr`,
//! the per-read EM runs the same rounds on the component layout of [`crate::em_components`]
//! (`OARFISH_EM_COMPONENTS=0` keeps the CSR layout); the truncation EM keeps the CSR layout.

use std::time::Instant;

use itertools::izip;
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use tracing::{info, span, trace};

use crate::em_components;
use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use crate::util::run_warnings::EmConvergence;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Impl {
    Local,
    Csr,
}

impl Impl {
    /// The per-read EM layout: `csr` by default; `OARFISH_READ_EM=local` selects the
    /// thread-local layout and `OARFISH_READ_EM=atomic` the atomic EM (None).
    pub fn from_env() -> Option<Impl> {
        match std::env::var("OARFISH_READ_EM").ok().as_deref() {
            Some("atomic") => None,
            Some("local") => Some(Impl::Local),
            _ => Some(Impl::Csr),
        }
    }
}

/// An alignment's fixed weight in the per-read EM: its score probability times, when those
/// models are on, its coverage probability and its length density.
#[inline]
pub(crate) fn aln_weight(
    a: &AlnInfo,
    p: f32,
    cp: f64,
    cov: bool,
    kde: Option<&kders::kde::KDEModel>,
    tinfo: &[TranscriptInfo],
) -> f64 {
    p as f64
        * if cov { cp } else { 1.0 }
        * match kde {
            Some(k) => {
                k[(
                    tinfo[a.ref_id as usize].lenf as usize,
                    a.alignment_span() as usize,
                )]
            }
            None => 1.0,
        }
}

/// The alignments as the EM reads them: straight from the store, in its order (which `bulk`
/// arranges by each read's heaviest target under [`aln_weight`], see [`read_order`]), with
/// each weight computed where it is used. Only the length-density weights, a table lookup,
/// are tabulated.
pub(crate) struct Flat<'a> {
    /// Read `r` owns alignments `off[r]..off[r + 1]`.
    pub off: &'a [usize],
    alns: &'a [AlnInfo],
    p: &'a [f32],
    /// Coverage probabilities, when the coverage model is on.
    cp: Option<&'a [f64]>,
    /// Every weight, when the length density is on.
    w: Option<Vec<f64>>,
}

impl<'a> Flat<'a> {
    fn new(em_info: &'a EMInfo) -> Self {
        let store = &em_info.eq_map;
        let cov = store.filter_opts.model_coverage;
        let w = em_info.kde_model.as_ref().map(|k| {
            (0..store.len())
                .flat_map(|r| {
                    let (alns, probs, cps) = store.read(r);
                    izip!(alns, probs, cps)
                        .map(|(a, &p, &cp)| aln_weight(a, p, cp, cov, Some(k), em_info.txp_info))
                        .collect::<Vec<_>>()
                })
                .collect()
        });
        Flat {
            off: store.boundaries(),
            alns: &store.alignments,
            p: &store.as_probabilities,
            cp: cov.then_some(&store.coverage_probabilities[..]),
            w,
        }
    }

    fn nreads(&self) -> usize {
        self.off.len() - 1
    }

    fn n(&self) -> usize {
        self.alns.len()
    }

    #[inline(always)]
    fn target(&self, k: usize) -> u32 {
        self.alns[k].ref_id
    }

    /// Alignment `k`'s weight; equal, bit for bit, to [`aln_weight`] (`x * 1.0 == x`).
    #[inline(always)]
    fn w(&self, k: usize) -> f64 {
        match &self.w {
            Some(w) => w[k],
            None => self.p[k] as f64 * self.cp.map_or(1.0, |cp| cp[k]),
        }
    }
}

/// Read ids `0..nreads` sorted by `key(read)`, the read's heaviest target: the reads of one
/// transcript get nearby ids, so the per-transcript pass reads their normalizers from nearby
/// memory. The EM does not depend on read order.
pub(crate) fn read_order(nreads: usize, key: impl Fn(usize) -> u32 + Sync) -> Vec<u32> {
    let mut order: Vec<(u32, u32)> = (0..nreads)
        .into_par_iter()
        .map(|r| (key(r), r as u32))
        .collect();
    order.par_sort_unstable();
    order.into_iter().map(|(_, r)| r).collect()
}

/// Transcript-major transpose of `Flat`: transcript `t` owns slots `off[t]..off[t + 1]`, each
/// a (read, weight) pair.
struct Transpose {
    off: Vec<usize>,
    read: Vec<u32>,
    w: Vec<f64>,
}

impl Transpose {
    fn new(f: &Flat, ntx: usize) -> Self {
        let mut off = vec![0usize; ntx + 1];
        for k in 0..f.n() {
            off[f.target(k) as usize + 1] += 1;
        }
        for i in 0..ntx {
            off[i + 1] += off[i];
        }
        let mut fill = off.clone();
        let n = f.n();
        let mut read = vec![0u32; n];
        let mut w = vec![0f64; n];
        for r in 0..f.nreads() {
            for k in f.off[r]..f.off[r + 1] {
                let t = f.target(k) as usize;
                read[fill[t]] = r as u32;
                w[fill[t]] = f.w(k);
                fill[t] += 1;
            }
        }
        Transpose { off, read, w }
    }
}

/// Read chunks of roughly equal alignment counts, `per_thread` per worker.
fn chunks(f: &Flat, nthreads: usize, per_thread: usize) -> Vec<(usize, usize)> {
    let n = (nthreads * per_thread).max(1);
    let target = f.n().div_ceil(n).max(1);
    let mut out = Vec::with_capacity(n);
    let mut start = 0;
    for r in 0..f.nreads() {
        if f.off[r + 1] - f.off[start] >= target {
            out.push((start, r + 1));
            start = r + 1;
        }
    }
    if start < f.nreads() {
        out.push((start, f.nreads()));
    }
    out
}

/// Abundance estimates (expected reads per target), and how the EM loop stopped.
pub fn em(em_info: &EMInfo, nthreads: usize, how: Impl) -> (Vec<f64>, EmConvergence) {
    let span = span!(tracing::Level::INFO, "em_read");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| run(em_info, nthreads, how))
}

fn run(em_info: &EMInfo, nthreads: usize, how: Impl) -> (Vec<f64>, EmConvergence) {
    let t0 = Instant::now();
    let tinfo = em_info.txp_info;
    let ntx = tinfo.len();
    let f = Flat::new(em_info);
    let total_weight = em_info.eq_map.num_aligned_reads() as f64;
    let mut prev = match em_info.init_abundances {
        Some(ref v) => v.clone(),
        None => vec![total_weight / ntx as f64; ntx],
    };
    if how == Impl::Csr && em_components::enabled() {
        let plan = em_components::Plan::new_by_size(|k| f.target(k), |k| f.w(k), f.off, ntx, true);
        info!(
            "per-read EM (components): {} reads, {} alignments, {} connected components, setup {:.2}s",
            f.nreads().to_formatted_string(&Locale::en),
            f.n().to_formatted_string(&Locale::en),
            plan.n_components().to_formatted_string(&Locale::en),
            t0.elapsed().as_secs_f64()
        );
        let t1 = Instant::now();
        let out = em_components::run_counted(
            &plan,
            em_components::Counts::Unit,
            &prev,
            em_info.max_iter,
            em_info.convergence_thresh,
            1,
        );
        info!(
            "per-read EM (components): {:.2}s",
            t1.elapsed().as_secs_f64()
        );
        return out;
    }
    let mut curr = vec![0.0; ntx];
    // Layout-specific state, built once.
    let tr = (how == Impl::Csr).then(|| Transpose::new(&f, ntx));
    let ch = chunks(&f, nthreads, 2);
    let mut bufs: Vec<Vec<f64>> = if how == Impl::Local {
        (0..ch.len()).map(|_| vec![0.0; ntx]).collect()
    } else {
        Vec::new()
    };
    let mut inv = vec![0.0f64; f.nreads()];
    info!(
        "per-read EM ({how:?}): {} reads, {} alignments, setup {:.2}s",
        f.nreads().to_formatted_string(&Locale::en),
        f.n().to_formatted_string(&Locale::en),
        t0.elapsed().as_secs_f64()
    );
    let t1 = Instant::now();
    let mut step = |prev: &[f64], curr: &mut [f64]| match how {
        Impl::Local => step_local(&f, &ch, &mut bufs, prev, curr),
        Impl::Csr => step_csr(&f, tr.as_ref().unwrap(), &mut inv, prev, curr),
    };
    let mut niter = 0u32;
    let (mut converged, mut last_rel_diff) = (false, f64::INFINITY);
    while niter < em_info.max_iter {
        step(&prev, &mut curr);
        let rel_diff = prev
            .par_iter()
            .zip(curr.par_iter())
            .filter(|(p, _)| **p > constants::MIN_READ_THRESH)
            .map(|(p, c)| (c - p) / p)
            .reduce(|| 0.0, f64::max);
        std::mem::swap(&mut prev, &mut curr);
        last_rel_diff = rel_diff;
        if rel_diff < em_info.convergence_thresh && niter > 1 {
            converged = true;
            break;
        }
        niter += 1;
        if niter.is_multiple_of(100) {
            info!(
                "iteration {}; rel diff {}",
                niter.to_formatted_string(&Locale::en),
                rel_diff
            );
        } else if niter.is_multiple_of(10) {
            trace!(
                "iteration {}; rel diff {}",
                niter.to_formatted_string(&Locale::en),
                rel_diff
            );
        }
    }
    prev.iter_mut().for_each(|x| {
        if *x < constants::MIN_READ_THRESH {
            *x = 0.0;
        }
    });
    step(&prev, &mut curr);
    info!(
        "per-read EM ({how:?}): {niter} iterations in {:.2}s",
        t1.elapsed().as_secs_f64()
    );
    let conv = EmConvergence {
        niter,
        max_iter: em_info.max_iter,
        converged,
        rel_diff: last_rel_diff,
        thresh: em_info.convergence_thresh,
    };
    (curr, conv)
}

/// Bootstrap replicates of the per-read EM on the component layout: each resamples the reads
/// with replacement, as [`crate::em::bootstrap`] does, and runs that EM's rounds (at least 50
/// before it may stop) with each read counted as often as it was drawn. One layout serves every
/// replicate.
pub fn bootstrap(em_info: &EMInfo, num_boot: u32, nthreads: usize) -> Vec<Vec<f64>> {
    let span = span!(tracing::Level::INFO, "bootstrap");
    let _guard = span.enter();
    info!("will collect {num_boot} bootstraps");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        let ntx = em_info.txp_info.len();
        let f = Flat::new(em_info);
        let plan = em_components::Plan::new_by_size(|k| f.target(k), |k| f.w(k), f.off, ntx, true);
        let total_weight = em_info.eq_map.num_aligned_reads() as f64;
        let init = match em_info.init_abundances {
            Some(ref v) => v.clone(),
            None => vec![total_weight / ntx as f64; ntx],
        };
        (0..num_boot)
            .into_par_iter()
            .map(|i| {
                info!("evaluating bootstrap replicate {}", i);
                let mut mult = vec![0.0_f64; f.nreads()];
                for r in crate::bootstrap::get_sample_inds(f.nreads(), &mut rand::rng()) {
                    mult[r] += 1.0;
                }
                em_components::run(
                    &plan,
                    em_components::Counts::Of(&mult),
                    &init,
                    em_info.max_iter,
                    em_info.convergence_thresh,
                    50,
                )
            })
            .collect()
    })
}

/// One iteration: each chunk accumulates into its own buffer; buffers are summed per block of
/// transcripts.
fn step_local(
    f: &Flat,
    ch: &[(usize, usize)],
    bufs: &mut [Vec<f64>],
    prev: &[f64],
    curr: &mut [f64],
) {
    ch.par_iter()
        .zip(bufs.par_iter_mut())
        .for_each(|(&(r0, r1), acc)| {
            acc.fill(0.0);
            for r in r0..r1 {
                let (lo, hi) = (f.off[r], f.off[r + 1]);
                let mut denom = 0.0;
                for k in lo..hi {
                    denom += prev[f.target(k) as usize] * f.w(k);
                }
                if denom > constants::EM_DENOM_THRESH {
                    let inv = 1.0 / denom;
                    for k in lo..hi {
                        let t = f.target(k) as usize;
                        acc[t] += prev[t] * f.w(k) * inv;
                    }
                }
            }
        });
    const BLOCK: usize = 4096;
    let bufs = &*bufs;
    curr.par_chunks_mut(BLOCK).enumerate().for_each(|(b, out)| {
        let (base, len) = (b * BLOCK, out.len());
        out.fill(0.0);
        for acc in bufs {
            for (o, a) in out.iter_mut().zip(&acc[base..base + len]) {
                *o += a;
            }
        }
    });
}

/// One iteration: per-read normalizers, then per-transcript sums over the transpose.
fn step_csr(f: &Flat, tr: &Transpose, inv: &mut [f64], prev: &[f64], curr: &mut [f64]) {
    inv.par_iter_mut().enumerate().for_each(|(r, x)| {
        let mut denom = 0.0;
        for k in f.off[r]..f.off[r + 1] {
            denom += prev[f.target(k) as usize] * f.w(k);
        }
        *x = if denom > constants::EM_DENOM_THRESH {
            1.0 / denom
        } else {
            0.0
        };
    });
    let inv = &*inv;
    curr.par_iter_mut().enumerate().for_each(|(t, c)| {
        let mut s = 0.0;
        for i in tr.off[t]..tr.off[t + 1] {
            s += tr.w[i] * inv[tr.read[i] as usize];
        }
        *c = prev[t] * s;
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::em_tests::{
        Classes, MIXED_NTX, Read, RefEm, Rng, assert_close, check_invariants, em_info, mixed_reads,
        on, ref_em, store_of, tinfo,
    };

    /// A store of the mixed fixture with the coverage model on: every alignment also carries
    /// a coverage probability, so its EM weight is score probability times coverage
    /// probability. Returns the reads with those weights folded in, for the reference.
    fn coverage_store<'h>(
        h: &'h noodles_sam::Header,
        txps: &mut [TranscriptInfo],
        seed: u64,
    ) -> (
        crate::util::oarfish_types::InMemoryAlignmentStore<'h>,
        Vec<Read>,
    ) {
        let reads = mixed_reads(seed);
        let mut store = store_of(h, &reads, txps);
        store.filter_opts.model_coverage = true;
        let mut rng = Rng::new(seed + 1);
        store.coverage_probabilities = (0..store.total_len())
            .map(|_| 0.2 + rng.below(800) as f64 / 1000.0)
            .collect();
        let mut k = 0;
        let weighted: Vec<Read> = reads
            .iter()
            .map(|r| {
                r.iter()
                    .map(|&(t, p)| {
                        let w = p as f64 * store.coverage_probabilities[k];
                        k += 1;
                        (t, w as f32)
                    })
                    .collect()
            })
            .collect();
        (store, weighted)
    }

    /// The reference classes of a store's reads under their exact f64 weights.
    fn ref_classes(f: &Flat, ntx: usize) -> Classes {
        Classes {
            ntx,
            cls: (0..f.nreads())
                .map(|r| {
                    let e = (f.off[r]..f.off[r + 1])
                        .map(|k| (f.target(k), f.w(k)))
                        .collect();
                    (e, 1.0)
                })
                .collect(),
        }
    }

    /// `read_order` is a permutation of the reads sorted by key, ties by read id.
    #[test]
    fn read_order_sorts_by_key_then_read() {
        let mut rng = Rng::new(40);
        let keys: Vec<u32> = (0..5000).map(|_| rng.below(37) as u32).collect();
        let order = read_order(keys.len(), |r| keys[r]);
        let mut seen = vec![false; keys.len()];
        for &r in &order {
            assert!(!std::mem::replace(&mut seen[r as usize], true));
        }
        assert!(seen.iter().all(|&s| s));
        assert!(
            order
                .windows(2)
                .all(|w| (keys[w[0] as usize], w[0]) < (keys[w[1] as usize], w[1]))
        );
        assert!(read_order(0, |_| 0).is_empty());
    }

    /// Chunks cover the reads in order, each but the last holding at least its share of the
    /// alignments; the transpose lists each transcript's alignments in read order.
    #[test]
    fn chunks_and_transpose_cover_the_alignments() {
        let h = noodles_sam::Header::default();
        let ntx = MIXED_NTX;
        let mut txps = tinfo(&vec![1000; ntx]);
        let (store, _) = coverage_store(&h, &mut txps, 41);
        let empty = store_of(&h, &[], &mut txps);
        for st in [&store, &empty] {
            let emi = em_info(st, &txps, 10, 1e-6);
            let f = Flat::new(&emi);
            for (threads, per) in [(1, 1), (1, 2), (3, 2), (8, 2), (64, 4)] {
                let ch = chunks(&f, threads, per);
                let share = f.n().div_ceil(threads * per).max(1);
                let mut next = 0;
                for (i, &(a, b)) in ch.iter().enumerate() {
                    assert_eq!(a, next);
                    assert!(b > a);
                    if i + 1 < ch.len() {
                        assert!(f.off[b] - f.off[a] >= share);
                    }
                    next = b;
                }
                assert_eq!(next, f.nreads());
            }
            let tr = Transpose::new(&f, ntx);
            for t in 0..ntx {
                let want: Vec<(u32, u64)> = (0..f.nreads())
                    .flat_map(|r| {
                        (f.off[r]..f.off[r + 1])
                            .filter(|&k| f.target(k) == t as u32)
                            .map(move |k| (r as u32, k))
                    })
                    .map(|(r, k)| (r, f.w(k).to_bits()))
                    .collect();
                let got: Vec<(u32, u64)> = (tr.off[t]..tr.off[t + 1])
                    .map(|i| (tr.read[i], tr.w[i].to_bits()))
                    .collect();
                assert_eq!(got, want, "t{t}");
            }
        }
    }

    /// An alignment's weight is its score probability, times its coverage probability when the
    /// coverage model is on; `Flat::w` gives the same bits.
    #[test]
    fn alignment_weights() {
        let h = noodles_sam::Header::default();
        let mut txps = tinfo(&vec![1000; MIXED_NTX]);
        let (mut store, _) = coverage_store(&h, &mut txps, 42);
        for cov in [true, false] {
            store.filter_opts.model_coverage = cov;
            let emi = em_info(&store, &txps, 10, 1e-6);
            let f = Flat::new(&emi);
            for r in 0..store.len() {
                let (alns, probs, cps) = store.read(r);
                for (j, (a, (&p, &cp))) in alns.iter().zip(probs.iter().zip(cps)).enumerate() {
                    let w = aln_weight(a, p, cp, cov, None, &txps);
                    assert_eq!(w, if cov { p as f64 * cp } else { p as f64 });
                    assert_eq!(w.to_bits(), f.w(store.boundaries()[r] + j).to_bits());
                }
            }
        }
    }

    /// The thread-local and transposed rounds are the reference round, on any number of
    /// threads; the transposed layout's loop (which `OARFISH_EM_COMPONENTS=0` selects) driven
    /// for a fixed number of rounds is the reference's.
    #[test]
    fn local_and_csr_rounds_are_the_reference_round() {
        let h = noodles_sam::Header::default();
        let ntx = MIXED_NTX;
        let mut txps = tinfo(&vec![1000; ntx]);
        let (store, _) = coverage_store(&h, &mut txps, 43);
        let emi = em_info(&store, &txps, 10, 1e-6);
        let f = Flat::new(&emi);
        let cl = ref_classes(&f, ntx);
        let mut rng = Rng::new(44);
        let prev: Vec<f64> = (0..ntx).map(|_| rng.below(30) as f64 * 0.5).collect();
        let want = cl.step(&prev);
        let tr = Transpose::new(&f, ntx);
        for threads in [1, 4, 16] {
            on(threads, || {
                let ch = chunks(&f, threads, 2);
                let mut bufs: Vec<Vec<f64>> = (0..ch.len()).map(|_| vec![0.0; ntx]).collect();
                let mut curr = vec![f64::NAN; ntx];
                step_local(&f, &ch, &mut bufs, &prev, &mut curr);
                assert_close(&curr, &want, 1e-12, &format!("local {threads}"));
                let mut inv = vec![0.0; f.nreads()];
                let mut curr = vec![f64::NAN; ntx];
                step_csr(&f, &tr, &mut inv, &prev, &mut curr);
                assert_close(&curr, &want, 1e-12, &format!("csr {threads}"));
            });
        }
        // 80 rounds, the pruning and the last round, as `run` runs them with thresh 0.
        let init = cl.uniform();
        let (want, _) = ref_em(&cl, &init, RefEm::new(80, 0.0, 1));
        let mut inv = vec![0.0; f.nreads()];
        let (mut prev, mut curr) = (init.clone(), vec![0.0; ntx]);
        for _ in 0..80 {
            step_csr(&f, &tr, &mut inv, &prev, &mut curr);
            std::mem::swap(&mut prev, &mut curr);
        }
        prev.iter_mut()
            .filter(|x| **x < constants::MIN_READ_THRESH)
            .for_each(|x| *x = 0.0);
        step_csr(&f, &tr, &mut inv, &prev, &mut curr);
        assert_close(&curr, &want, 1e-10, "csr loop");
    }

    /// With the coverage model's weights, every per-read EM (and the score-only EMs, which
    /// apply the same weights) reaches the reference estimates.
    #[test]
    fn coverage_weighted_ems_reach_the_reference() {
        let h = noodles_sam::Header::default();
        let ntx = MIXED_NTX;
        let mut txps = tinfo(&vec![1000; ntx]);
        let (store, weighted) = coverage_store(&h, &mut txps, 45);
        let emi = em_info(&store, &txps, 3000, 1e-10);
        let cl = ref_classes(&Flat::new(&emi), ntx);
        let want = |m| ref_em(&cl, &cl.uniform(), RefEm::new(3000, 1e-10, m)).0;
        let (w1, w50) = (want(1), want(50));
        check_invariants(&w1, &weighted, None, ntx, "reference");
        for (name, got, w) in [
            ("Local", em(&emi, 3, Impl::Local).0, &w1),
            ("Csr 1", em(&emi, 1, Impl::Csr).0, &w1),
            ("Csr 16", em(&emi, 16, Impl::Csr).0, &w1),
            ("em_par", crate::em::em_par(&emi, 4).0, &w1),
            ("em", crate::em::em(&emi, 1).0, &w50),
        ] {
            assert_close(&got, w, 1e-7, name);
            check_invariants(&got, &weighted, None, ntx, name);
        }
        for b in bootstrap(&emi, 2, 4) {
            let sum: f64 = b.iter().sum();
            assert!((sum - store.len() as f64).abs() < 1e-9 * store.len() as f64);
            assert!(b.iter().all(|&x| x >= 0.0) && b[62] == 0.0 && b[63] == 0.0);
        }
    }
}
