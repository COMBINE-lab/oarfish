//! EM with a positional likelihood for 3'-anchored, possibly 5'-truncated long reads.
//!
//! Long-read protocols sequence one read per captured molecule, so read counts need no
//! fragmentation-style length correction. What does depend on the transcript is where a read
//! can land on it. Reads end near the transcript's 3' end (poly(A) priming, or direct RNA
//! sequencing from the 3' end), and they extend toward the 5' end until reverse transcription,
//! degradation or the pore stops them, at a roughly constant rate per base. For an alignment of a
//! read to transcript `t` of length `len` covering `[start, end)`, with `ext = len - start` the
//! read's extent from the transcript's 3' end:
//!
//! ```text
//! P(start | t) = [start < BW] (w + (1 - w) e^(-h len)) / BW     reached the 5' end
//!              + (1 - w) h e^(-h ext)                           stopped after ext bases
//! P(end | t)   = r [len - end < BW] / BW + (1 - r) g e^(-g (len - end))
//! ```
//!
//! `w` is the fraction of molecules that are full length regardless of length (full-length
//! selection in some protocols), `h` the per-base stop rate, `r` the fraction of reads ending at
//! the 3' end, and `g` the decay of the rest. For each transcript the start term is a
//! distribution over the starts it can produce, so a long transcript pays for its many possible
//! truncation points: a read that is full length for a short isoform but truncated for a longer
//! one favors the short one, by as much as the learned stop rate says truncation is unlikely.
//! The parametric (exponential) form keeps the model from explaining reads shared by two
//! isoforms as truncations that all stop at the same base, which a free-form truncation
//! histogram can learn. The four parameters are fitted in the EM together with the abundances
//! (censored-exponential M-steps); each alignment's likelihood is its alignment-score
//! probability times `P(start | t) P(end | t)`.

use std::sync::atomic::Ordering;

use atomic_float::AtomicF64;
use num_format::{Locale, ToFormattedString};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use tracing::{info, span, trace};

use crate::util::constants;
use crate::util::oarfish_types::{EMInfo, TranscriptInfo};

/// Width (nt) of the windows at the transcript ends that count as "at the end".
const BW: f64 = 50.0;

/// Alignments flattened with the positional quantities the model needs.
struct Flat {
    /// Read `r` owns alignments `off[r]..off[r + 1]`.
    off: Vec<usize>,
    target: Vec<u32>,
    prob: Vec<f64>,
    len: Vec<f64>,
    /// Extent from the alignment's start to the transcript's 3' end.
    ext: Vec<f64>,
    /// Distance from the alignment's end to the transcript's 3' end.
    d3: Vec<f64>,
}

impl Flat {
    fn push(&mut self, t: u32, p: f64, len: f64, start: f64, end: f64) {
        let start = start.clamp(0.0, (len - 1.0).max(0.0));
        let end = end.clamp(start, len);
        self.target.push(t);
        self.prob.push(p);
        self.len.push(len);
        self.ext.push(len - start);
        self.d3.push(len - end);
    }

    fn new(em_info: &EMInfo) -> Self {
        let tinfo = em_info.txp_info;
        let mut f = Flat::empty();
        for (alns, probs, _) in em_info.eq_map.iter() {
            for (a, &p) in alns.iter().zip(probs) {
                let len = tinfo[a.ref_id as usize].lenf;
                f.push(a.ref_id, p as f64, len, a.start as f64, a.end as f64);
            }
            f.off.push(f.target.len());
        }
        f
    }

    fn empty() -> Self {
        Flat {
            off: vec![0],
            target: vec![],
            prob: vec![],
            len: vec![],
            ext: vec![],
            d3: vec![],
        }
    }
}

/// The positional model's parameters.
#[derive(Clone, Copy, Debug)]
struct Params {
    /// Fraction of molecules full length regardless of length.
    w: f64,
    /// Per-base stop rate toward the 5' end.
    h: f64,
    /// Fraction of reads ending at the 3' end, and the decay of the others' distance from it.
    r: f64,
    g: f64,
}

/// Per-alignment terms of the positional likelihood.
struct Terms {
    /// Reached the 5' end because the molecule is full length by protocol (`w`), or by not
    /// stopping (`(1 - w) e^(-h len)`); stopped after `ext` bases.
    full_w: f64,
    full_surv: f64,
    stopped: f64,
    /// Ends at the 3' end; ends `d3` bases before it.
    at3: f64,
    off3: f64,
}

impl Terms {
    fn start(&self) -> f64 {
        self.full_w + self.full_surv + self.stopped
    }
    fn end(&self) -> f64 {
        self.at3 + self.off3
    }
}

impl Params {
    fn initial() -> Self {
        Params {
            w: 0.5,
            h: 1e-3,
            r: 0.5,
            g: 2e-3,
        }
    }

    #[inline]
    fn terms(&self, f: &Flat, k: usize) -> Terms {
        let full = f.ext[k] > f.len[k] - BW;
        Terms {
            full_w: if full { self.w / BW } else { 0.0 },
            full_surv: if full {
                (1.0 - self.w) * (-self.h * f.len[k]).exp() / BW
            } else {
                0.0
            },
            stopped: (1.0 - self.w) * self.h * (-self.h * f.ext[k]).exp(),
            at3: if f.d3[k] < BW { self.r / BW } else { 0.0 },
            off3: (1.0 - self.r) * self.g * (-self.g * f.d3[k]).exp(),
        }
    }
}

/// Expected sufficient statistics for the M-step, accumulated per thread.
#[derive(Default)]
struct Stats {
    /// Posterior mass, and its parts: full length by protocol, full length by not stopping,
    /// stopped; ending at the 3' end, ending before it.
    n: f64,
    n_w: f64,
    n_surv: f64,
    n_stop: f64,
    n_at3: f64,
    n_off3: f64,
    /// Exposure (bases) for the stop rate and for the 3' decay.
    expo_h: f64,
    expo_g: f64,
}

impl Stats {
    fn merge(mut self, o: Stats) -> Stats {
        self.n += o.n;
        self.n_w += o.n_w;
        self.n_surv += o.n_surv;
        self.n_stop += o.n_stop;
        self.n_at3 += o.n_at3;
        self.n_off3 += o.n_off3;
        self.expo_h += o.expo_h;
        self.expo_g += o.expo_g;
        self
    }

    fn add(&mut self, f: &Flat, k: usize, t: &Terms, post: f64) {
        let (s, e) = (t.start(), t.end());
        if s <= 0.0 || e <= 0.0 {
            return;
        }
        let (pw, psurv, pstop) = (
            post * t.full_w / s,
            post * t.full_surv / s,
            post * t.stopped / s,
        );
        let (pat, poff) = (post * t.at3 / e, post * t.off3 / e);
        self.n += post;
        self.n_w += pw;
        self.n_surv += psurv;
        self.n_stop += pstop;
        // A surviving read is right-censored at the transcript length; a stopped one is an
        // observed stop after `ext` bases.
        self.expo_h += psurv * f.len[k] + pstop * f.ext[k];
        self.n_at3 += pat;
        self.n_off3 += poff;
        self.expo_g += poff * f.d3[k];
    }

    fn params(&self, old: Params) -> Params {
        if self.n <= 0.0 {
            return old;
        }
        let clamp = |x: f64| x.clamp(1e-6, 1.0 - 1e-6);
        Params {
            w: clamp(self.n_w / self.n),
            h: if self.expo_h > 0.0 {
                (self.n_stop / self.expo_h).clamp(1e-7, 1.0)
            } else {
                old.h
            },
            r: clamp(self.n_at3 / self.n),
            g: if self.expo_g > 0.0 {
                (self.n_off3 / self.expo_g).clamp(1e-7, 1.0)
            } else {
                old.g
            },
        }
    }
}

/// One EM round: distribute each read over its alignments, accumulating abundances into `curr`
/// and the positional statistics into the returned `Stats`.
fn step(f: &Flat, par: &Params, prev: &[AtomicF64], curr: &[AtomicF64]) -> Stats {
    let nreads = f.off.len() - 1;
    (0..nreads)
        .into_par_iter()
        .fold(Stats::default, |mut st, r| {
            let (lo, hi) = (f.off[r], f.off[r + 1]);
            let mut denom = 0.0;
            for k in lo..hi {
                let t = par.terms(f, k);
                denom += prev[f.target[k] as usize].load(Ordering::Relaxed)
                    * f.prob[k]
                    * t.start()
                    * t.end();
            }
            if denom > constants::EM_DENOM_THRESH {
                for k in lo..hi {
                    let t = par.terms(f, k);
                    let tid = f.target[k] as usize;
                    let post =
                        prev[tid].load(Ordering::Relaxed) * f.prob[k] * t.start() * t.end() / denom;
                    curr[tid].fetch_add(post, Ordering::AcqRel);
                    st.add(f, k, &t, post);
                }
            }
            st
        })
        .reduce(Stats::default, Stats::merge)
}

/// Abundance estimates (expected reads per target) under the truncation-aware model.
pub fn em(em_info: &EMInfo, nthreads: usize) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em_truncation");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| run(em_info))
}

fn run(em_info: &EMInfo) -> Vec<f64> {
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let f = Flat::new(em_info);
    let total_weight = em_info.eq_map.num_aligned_reads() as f64;
    let init = match em_info.init_abundances {
        Some(ref v) => v.clone(),
        None => vec![total_weight / (tinfo.len() as f64); tinfo.len()],
    };
    let mut prev: Vec<AtomicF64> = init.into_iter().map(AtomicF64::new).collect();
    let mut curr: Vec<AtomicF64> = (0..tinfo.len()).map(|_| AtomicF64::new(0.0)).collect();
    let mut par = Params::initial();

    let mut niter = 0_u32;
    while niter < em_info.max_iter {
        par = step(&f, &par, &prev, &curr).params(par);
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
        if rel_diff < em_info.convergence_thresh && niter > 1 {
            break;
        }
        niter += 1;
        if niter.is_multiple_of(100) {
            info!(
                "iteration {}; rel diff {}; {:?}",
                niter.to_formatted_string(&Locale::en),
                rel_diff,
                par
            );
        } else if niter.is_multiple_of(10) {
            trace!(
                "iteration {}; rel diff {}",
                niter.to_formatted_string(&Locale::en),
                rel_diff
            );
        }
    }
    info!(
        "truncation model: full length by protocol {:.3}; stop rate {:.2e}/nt (mean extent {:.0} nt); \
         ending within {} nt of the 3' end {:.3}",
        par.w,
        par.h,
        1.0 / par.h,
        BW as u32,
        par.r
    );

    // Zero very small abundances, then one more round with the parameters fixed, as the
    // per-read EM does.
    for x in &prev {
        if x.load(Ordering::Relaxed) < constants::MIN_READ_THRESH {
            x.store(0.0, Ordering::Relaxed);
        }
    }
    step(&f, &par, &prev, &curr);
    curr.iter().map(|x| x.load(Ordering::Relaxed)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat(alns: &[Vec<(u32, f64, f64, f64)>], lens: &[f64]) -> Flat {
        let mut f = Flat::empty();
        for r in alns {
            for &(t, p, s, e) in r {
                f.push(t, p, lens[t as usize], s, e);
            }
            f.off.push(f.target.len());
        }
        f
    }

    fn fit(f: &Flat, ntx: usize, iters: usize) -> (Vec<f64>, Params) {
        let n = (f.off.len() - 1) as f64;
        let mut prev: Vec<AtomicF64> = (0..ntx).map(|_| AtomicF64::new(n / ntx as f64)).collect();
        let mut curr: Vec<AtomicF64> = (0..ntx).map(|_| AtomicF64::new(0.0)).collect();
        let mut par = Params::initial();
        for _ in 0..iters {
            par = step(f, &par, &prev, &curr).params(par);
            std::mem::swap(&mut prev, &mut curr);
            curr.iter().for_each(|x| x.store(0.0, Ordering::Relaxed));
        }
        (
            prev.iter().map(|x| x.load(Ordering::Relaxed)).collect(),
            par,
        )
    }

    /// A short isoform (1,000 nt) is the 3' half of a long one (2,000 nt). 300 reads are full
    /// length on the long isoform, 100 more are truncated on it, and 600 are full length on the
    /// short one, which are also compatible with the long one as truncated reads. The score-only
    /// EM can't tell the shared reads apart; the positional model recovers most of them.
    #[test]
    fn full_length_reads_favor_the_isoform_they_span() {
        let lens = [2000.0, 1000.0];
        let mut reads = Vec::new();
        for _ in 0..300 {
            reads.push(vec![(0, 1.0, 0.0, 2000.0)]);
        }
        for i in 0..100 {
            let s = 200.0 + 6.0 * i as f64; // truncated starts on the long isoform's 5' half
            reads.push(vec![(0, 1.0, s, 2000.0)]);
        }
        for _ in 0..600 {
            // Full length on the short isoform: transcript coordinates 0..1000 on it,
            // 1000..2000 on the long one.
            reads.push(vec![(0, 1.0, 1000.0, 2000.0), (1, 1.0, 0.0, 1000.0)]);
        }
        let f = flat(&reads, &lens);
        let (est, par) = fit(&f, 2, 300);
        assert!((est[0] + est[1] - 1000.0).abs() < 1e-6);
        assert!(est[1] > 500.0, "{est:?} {par:?}");
        assert!(par.r > 0.9, "{par:?}");
    }

    /// With every read full length and unambiguous, estimates are the read counts.
    #[test]
    fn unambiguous_reads_are_counted() {
        let lens = [1500.0, 800.0, 3000.0];
        let mut reads = Vec::new();
        for (t, n) in [(0u32, 50), (1, 120), (2, 30)] {
            for _ in 0..n {
                reads.push(vec![(t, 1.0, 0.0, lens[t as usize])]);
            }
        }
        let (est, _) = fit(&flat(&reads, &lens), 3, 50);
        for (e, want) in est.iter().zip([50.0, 120.0, 30.0]) {
            assert!((e - want).abs() < 1e-6, "{est:?}");
        }
    }
}
