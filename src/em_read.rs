//! Per-read EM without atomics: experimental alternatives to `em::em_par`.
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

use std::time::Instant;

use itertools::izip;
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use tracing::{info, span, trace};

use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Impl {
    Local,
    Csr,
}

impl Impl {
    /// From `OARFISH_READ_EM` (`local` or `csr`); None keeps the atomic EM.
    pub fn from_env() -> Option<Impl> {
        match std::env::var("OARFISH_READ_EM").ok()?.as_str() {
            "local" => Some(Impl::Local),
            "csr" => Some(Impl::Csr),
            _ => None,
        }
    }
}

/// Alignments flattened in read order, with their fixed weights.
pub(crate) struct Flat {
    /// Read `r` owns alignments `off[r]..off[r + 1]`.
    pub off: Vec<usize>,
    pub target: Vec<u32>,
    pub w: Vec<f64>,
}

impl Flat {
    /// Built in [`read_order`], with capacities reserved up front, so the table exists once.
    fn new(em_info: &EMInfo) -> Self {
        let tinfo: &[TranscriptInfo] = em_info.txp_info;
        let store = &em_info.eq_map;
        let cov = store.filter_opts.model_coverage;
        let weight = |a: &AlnInfo, p: f32, cp: f64| {
            p as f64
                * if cov { cp } else { 1.0 }
                * match em_info.kde_model {
                    Some(ref k) => {
                        k[(
                            tinfo[a.ref_id as usize].lenf as usize,
                            a.alignment_span() as usize,
                        )]
                    }
                    None => 1.0,
                }
        };
        let order = read_order(store.len(), |r| {
            let (alns, probs, cps) = store.read(r);
            izip!(alns, probs, cps)
                .map(|(a, &p, &cp)| (weight(a, p, cp), a.ref_id))
                .max_by(|x, y| x.0.total_cmp(&y.0))
                .map_or(0, |x| x.1)
        });
        let n = store.alignments.len();
        let mut f = Flat {
            off: Vec::with_capacity(order.len() + 1),
            target: Vec::with_capacity(n),
            w: Vec::with_capacity(n),
        };
        f.off.push(0);
        for &r in &order {
            let (alns, probs, cps) = store.read(r as usize);
            for (a, &p, &cp) in izip!(alns, probs, cps) {
                f.target.push(a.ref_id);
                f.w.push(weight(a, p, cp));
            }
            f.off.push(f.target.len());
        }
        f
    }

    fn nreads(&self) -> usize {
        self.off.len() - 1
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
        for &t in &f.target {
            off[t as usize + 1] += 1;
        }
        for i in 0..ntx {
            off[i + 1] += off[i];
        }
        let mut fill = off.clone();
        let n = f.target.len();
        let mut read = vec![0u32; n];
        let mut w = vec![0f64; n];
        for r in 0..f.nreads() {
            for k in f.off[r]..f.off[r + 1] {
                let t = f.target[k] as usize;
                read[fill[t]] = r as u32;
                w[fill[t]] = f.w[k];
                fill[t] += 1;
            }
        }
        Transpose { off, read, w }
    }
}

/// Read chunks of roughly equal alignment counts, `per_thread` per worker.
fn chunks(f: &Flat, nthreads: usize, per_thread: usize) -> Vec<(usize, usize)> {
    let n = (nthreads * per_thread).max(1);
    let target = f.target.len().div_ceil(n).max(1);
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

pub fn em(em_info: &EMInfo, nthreads: usize, how: Impl) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em_read");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| run(em_info, nthreads, how))
}

fn run(em_info: &EMInfo, nthreads: usize, how: Impl) -> Vec<f64> {
    let t0 = Instant::now();
    let tinfo = em_info.txp_info;
    let ntx = tinfo.len();
    let f = Flat::new(em_info);
    let total_weight = em_info.eq_map.num_aligned_reads() as f64;
    let mut prev = match em_info.init_abundances {
        Some(ref v) => v.clone(),
        None => vec![total_weight / ntx as f64; ntx],
    };
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
        f.target.len().to_formatted_string(&Locale::en),
        t0.elapsed().as_secs_f64()
    );
    let t1 = Instant::now();
    let mut step = |prev: &[f64], curr: &mut [f64]| match how {
        Impl::Local => step_local(&f, &ch, &mut bufs, prev, curr),
        Impl::Csr => step_csr(&f, tr.as_ref().unwrap(), &mut inv, prev, curr),
    };
    let mut niter = 0u32;
    while niter < em_info.max_iter {
        step(&prev, &mut curr);
        let rel_diff = prev
            .par_iter()
            .zip(curr.par_iter())
            .filter(|(p, _)| **p > constants::MIN_READ_THRESH)
            .map(|(p, c)| (c - p) / p)
            .reduce(|| 0.0, f64::max);
        std::mem::swap(&mut prev, &mut curr);
        if rel_diff < em_info.convergence_thresh && niter > 1 {
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
    curr
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
                    denom += prev[f.target[k] as usize] * f.w[k];
                }
                if denom > constants::EM_DENOM_THRESH {
                    let inv = 1.0 / denom;
                    for k in lo..hi {
                        let t = f.target[k] as usize;
                        acc[t] += prev[t] * f.w[k] * inv;
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
            denom += prev[f.target[k] as usize] * f.w[k];
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
