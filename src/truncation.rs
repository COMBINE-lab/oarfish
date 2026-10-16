//! EM with a positional likelihood: a mixture of the ways long reads can sit on a transcript.
//!
//! Long-read protocols sequence one read per captured molecule, so read counts need no
//! fragmentation-style length correction. What does depend on the transcript is where a read
//! can land on it, and that differs by protocol: direct RNA and oligo(dT)-primed cDNA give reads
//! anchored at the 3' end and truncated toward the 5' end; some library preparations keep the
//! 5' end and lose the 3' end; full-length selection gives reads spanning the whole molecule;
//! and fragmentation (or a simulator that places reads uniformly) gives reads anywhere on it.
//! Which of these a sample contains is learned from the sample, as a mixture.
//!
//! For an alignment covering `[s, e)` of a transcript of length `L`, with read length
//! `l = e - s` and gaps `d5 = s`, `d3 = L - e`, let `F5 = [d5 < BW]`, `F3 = [d3 < BW]`. Every
//! component draws a read length from an exponential law, censored at the transcript length (a
//! read that would be longer spans the transcript), and places the read by its own rule:
//!
//! ```text
//! 3'-anchored  F3 / BW · ( h3 e^(-h3 l)  +  F5 e^(-h3 L) / BW )
//! 5'-anchored  F5 / BW · ( h5 e^(-h5 l)  +  F3 e^(-h5 L) / BW )
//! full length  F3 F5 / BW²
//! uniform      lam e^(-lam l) / (L - l + 1)  +  F3 F5 e^(-lam L) / BW²
//! P(s, e | t) = p3 · 3'-anchored + p5 · 5'-anchored + pf · full length + pu · uniform
//! ```
//!
//! Each component sums to one over the placements a transcript of length `L` allows, so a
//! transcript is never favored merely because one of its ends happens to coincide with a read's
//! end: under the uniform component every placement of a read of length `l` is equally likely,
//! and a candidate is judged only by how much room it gave the read. When the reads are not
//! anchored (`pu` near one) the model reduces to uniform placement; when they are (`p3` near one)
//! it is a 3'-anchored truncation model. The seven parameters are fitted in the EM together with
//! the abundances (mixture weights from the components' posterior shares, rates by
//! censored-exponential M-steps); each alignment's likelihood is its alignment-score probability
//! times `P(s, e | t)`.

use std::sync::atomic::Ordering;

use atomic_float::AtomicF64;
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use tracing::{info, span, trace};

use crate::presence::{Presence, PresenceOpts};
use crate::util::constants;
use crate::util::oarfish_types::{EMInfo, TranscriptInfo};

/// Width (nt) of the windows at the transcript ends that count as "at the end".
const BW: f64 = 50.0;

/// Alignments flattened with the positional quantities the model needs.
struct Flat {
    /// Read `r` owns alignments `off[r]..off[r + 1]`.
    off: Vec<usize>,
    target: Vec<u32>,
    prob: Vec<f32>,
    /// Distances in bases (integers, so that [`Exps`] can tabulate their exponentials): from
    /// the alignment's start to the transcript's 3' end, and from its end to the 3' end.
    ext: Vec<u32>,
    d3: Vec<u32>,
    /// Transcript lengths, per transcript.
    tlen: Vec<u32>,
}

impl Flat {
    fn push(&mut self, t: u32, p: f32, start: f64, end: f64) {
        let len = self.tlen[t as usize] as f64;
        let start = start.clamp(0.0, (len - 1.0).max(0.0));
        let end = end.clamp(start, len);
        self.target.push(t);
        self.prob.push(p);
        self.ext.push((len - start).round() as u32);
        self.d3.push((len - end).round() as u32);
    }

    #[inline]
    fn len(&self, k: usize) -> u32 {
        self.tlen[self.target[k] as usize]
    }

    /// Aligned length `l` and the gap `d5` from the transcript's 5' end.
    #[inline]
    fn l_d5(&self, k: usize) -> (u32, u32) {
        let (len, ext, d3) = (self.len(k), self.ext[k], self.d3[k]);
        (ext.saturating_sub(d3), len.saturating_sub(ext))
    }

    fn max_len(&self) -> u32 {
        self.tlen.iter().copied().max().unwrap_or(0)
    }

    /// Reads in the order of `crate::em_read::read_order`, by most probable target.
    fn new(em_info: &EMInfo) -> Self {
        let store = &em_info.eq_map;
        let order = crate::em_read::read_order(store.len(), |r| {
            let (alns, probs, _) = store.read(r);
            alns.iter()
                .zip(probs)
                .max_by(|x, y| x.1.total_cmp(y.1))
                .map_or(0, |x| x.0.ref_id)
        });
        let n = store.alignments.len();
        let mut f = Flat {
            off: Vec::with_capacity(order.len() + 1),
            target: Vec::with_capacity(n),
            prob: Vec::with_capacity(n),
            ext: Vec::with_capacity(n),
            d3: Vec::with_capacity(n),
            tlen: em_info
                .txp_info
                .iter()
                .map(|t| t.lenf.round() as u32)
                .collect(),
        };
        f.off.push(0);
        for &r in &order {
            let (alns, probs, _) = store.read(r as usize);
            for (a, &p) in alns.iter().zip(probs) {
                f.push(a.ref_id, p, a.start as f64, a.end as f64);
            }
            f.off.push(f.target.len());
        }
        f
    }

    #[cfg(test)]
    fn empty(lens: &[f64]) -> Self {
        Flat {
            off: vec![0],
            target: vec![],
            prob: vec![],
            ext: vec![],
            d3: vec![],
            tlen: lens.iter().map(|l| l.round() as u32).collect(),
        }
    }
}

/// The positional model's parameters.
#[derive(Clone, Copy, Debug)]
struct Params {
    /// Mixture weights: 3'-anchored, 5'-anchored, full length, uniformly placed.
    p3: f64,
    p5: f64,
    pf: f64,
    pu: f64,
    /// Per-base stop rates of the 3'-anchored and 5'-anchored components, and the rate of the
    /// uniform component's read-length law.
    h3: f64,
    h5: f64,
    lam: f64,
}

/// Per-alignment terms of the positional likelihood, each already weighted by its mixture
/// weight. `*_stop` is the part where the read's length was observed, `*_cens` the part where
/// it spans the transcript (length censored at `L`).
struct Terms {
    a3_stop: f64,
    a3_cens: f64,
    a5_stop: f64,
    a5_cens: f64,
    full: f64,
    u_stop: f64,
    u_cens: f64,
}

impl Terms {
    #[inline]
    fn total(&self) -> f64 {
        self.a3_stop
            + self.a3_cens
            + self.a5_stop
            + self.a5_cens
            + self.full
            + self.u_stop
            + self.u_cens
    }
}

impl Params {
    fn initial() -> Self {
        Params {
            p3: 0.25,
            p5: 0.25,
            pf: 0.25,
            pu: 0.25,
            h3: 1e-3,
            h5: 1e-3,
            lam: 1e-3,
        }
    }

    #[inline]
    fn terms(&self, f: &Flat, k: usize, e: &Exps) -> Terms {
        let len = f.len(k);
        let (l, d5) = f.l_d5(k);
        let f3 = (f.d3[k] as f64) < BW;
        let f5 = (d5 as f64) < BW;
        let both = f3 && f5;
        let (l_us, len_us) = (l as usize, len as usize);
        let bw2 = BW * BW;
        Terms {
            a3_stop: if f3 {
                self.p3 * self.h3 * e.h3[l_us] / BW
            } else {
                0.0
            },
            a3_cens: if both {
                self.p3 * e.h3[len_us] / bw2
            } else {
                0.0
            },
            a5_stop: if f5 {
                self.p5 * self.h5 * e.h5[l_us] / BW
            } else {
                0.0
            },
            a5_cens: if both {
                self.p5 * e.h5[len_us] / bw2
            } else {
                0.0
            },
            full: if both { self.pf / bw2 } else { 0.0 },
            u_stop: self.pu * self.lam * e.lam[l_us] / ((len - l) as f64 + 1.0),
            u_cens: if both {
                self.pu * e.lam[len_us] / bw2
            } else {
                0.0
            },
        }
    }
}

/// `exp(-rate x)` for each rate and every whole number of bases `x` up to the longest
/// transcript, built once per iteration so that no alignment computes an exponential.
struct Exps {
    h3: Vec<f64>,
    h5: Vec<f64>,
    lam: Vec<f64>,
}

impl Exps {
    fn new(par: &Params, max_len: u32) -> Self {
        let tab = |rate: f64| {
            (0..=max_len)
                .into_par_iter()
                .map(|x| (-rate * x as f64).exp())
                .collect()
        };
        Exps {
            h3: tab(par.h3),
            h5: tab(par.h5),
            lam: tab(par.lam),
        }
    }
}

/// Expected sufficient statistics for the M-step, accumulated per thread.
#[derive(Default)]
struct Stats {
    /// Posterior mass, and its share in each component.
    n: f64,
    n3: f64,
    n5: f64,
    nf: f64,
    nu: f64,
    /// Observed stops (events) and exposure (bases) for each rate.
    ev_h3: f64,
    ex_h3: f64,
    ev_h5: f64,
    ex_h5: f64,
    ev_lam: f64,
    ex_lam: f64,
}

impl Stats {
    fn merge(mut self, o: Stats) -> Stats {
        self.n += o.n;
        self.n3 += o.n3;
        self.n5 += o.n5;
        self.nf += o.nf;
        self.nu += o.nu;
        self.ev_h3 += o.ev_h3;
        self.ex_h3 += o.ex_h3;
        self.ev_h5 += o.ev_h5;
        self.ex_h5 += o.ex_h5;
        self.ev_lam += o.ev_lam;
        self.ex_lam += o.ex_lam;
        self
    }

    fn add(&mut self, f: &Flat, k: usize, t: &Terms, post: f64) {
        let tot = t.total();
        if tot <= 0.0 {
            return;
        }
        let q = post / tot;
        let (l, _) = f.l_d5(k);
        let (l, len) = (l as f64, f.len(k) as f64);
        let (a3s, a3c) = (q * t.a3_stop, q * t.a3_cens);
        let (a5s, a5c) = (q * t.a5_stop, q * t.a5_cens);
        let (us, uc) = (q * t.u_stop, q * t.u_cens);
        self.n += post;
        self.n3 += a3s + a3c;
        self.n5 += a5s + a5c;
        self.nf += q * t.full;
        self.nu += us + uc;
        // A read spanning the transcript is right-censored at its length; any other read is
        // an observed stop after `l` bases.
        self.ev_h3 += a3s;
        self.ex_h3 += a3s * l + a3c * len;
        self.ev_h5 += a5s;
        self.ex_h5 += a5s * l + a5c * len;
        self.ev_lam += us;
        self.ex_lam += us * l + uc * len;
    }

    fn params(&self, old: Params) -> Params {
        if self.n <= 0.0 {
            return old;
        }
        let frac = |x: f64| (x / self.n).clamp(1e-6, 1.0);
        let (p3, p5, pf, pu) = (frac(self.n3), frac(self.n5), frac(self.nf), frac(self.nu));
        let z = p3 + p5 + pf + pu;
        let rate = |ev: f64, ex: f64, old: f64| {
            if ex > 0.0 && ev > 0.0 {
                (ev / ex).clamp(1e-7, 1.0)
            } else {
                old
            }
        };
        Params {
            p3: p3 / z,
            p5: p5 / z,
            pf: pf / z,
            pu: pu / z,
            h3: rate(self.ev_h3, self.ex_h3, old.h3),
            h5: rate(self.ev_h5, self.ex_h5, old.h5),
            lam: rate(self.ev_lam, self.ex_lam, old.lam),
        }
    }
}

/// One EM round: distribute each read over its alignments, accumulating abundances into `curr`
/// and the positional statistics into the returned `Stats`. With `mult`, read `r` counts
/// `mult[r]` times (a bootstrap replicate).
fn step(
    f: &Flat,
    par: &Params,
    e: &Exps,
    prev: &[AtomicF64],
    curr: &[AtomicF64],
    mult: Option<&[u32]>,
) -> Stats {
    let nreads = f.off.len() - 1;
    (0..nreads)
        .into_par_iter()
        .fold(Stats::default, |mut st, r| {
            let m = mult.map_or(1.0, |m| m[r] as f64);
            if m == 0.0 {
                return st;
            }
            let (lo, hi) = (f.off[r], f.off[r + 1]);
            let mut denom = 0.0;
            for k in lo..hi {
                let t = par.terms(f, k, e);
                denom += prev[f.target[k] as usize].load(Ordering::Relaxed)
                    * f.prob[k] as f64
                    * t.total();
            }
            if denom > constants::EM_DENOM_THRESH {
                for k in lo..hi {
                    let t = par.terms(f, k, e);
                    let tid = f.target[k] as usize;
                    let post = m * prev[tid].load(Ordering::Relaxed) * f.prob[k] as f64 * t.total()
                        / denom;
                    curr[tid].fetch_add(post, Ordering::AcqRel);
                    st.add(f, k, &t, post);
                }
            }
            st
        })
        .reduce(Stats::default, Stats::merge)
}

/// Abundance estimates (expected reads per target) under the truncation-aware model.
pub fn em(em_info: &EMInfo, presence: Option<PresenceOpts>, nthreads: usize) -> Vec<f64> {
    let span = span!(tracing::Level::INFO, "em_truncation");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        let (counts, par) = run(em_info, &Flat::new(em_info), None, presence);
        info!(
            "truncation model: 3'-anchored {:.3} (stop rate {:.2e}/nt), 5'-anchored {:.3} \
             (stop rate {:.2e}/nt), full length {:.3}, uniformly placed {:.3} (mean length {:.0} nt)",
            par.p3,
            par.h3,
            par.p5,
            par.h5,
            par.pf,
            par.pu,
            1.0 / par.lam
        );
        counts
    })
}

/// Bootstrap replicates under the truncation-aware model. Each replicate resamples the reads
/// with replacement and re-fits both the abundances and the positional parameters, so the
/// parameters' uncertainty propagates into the replicates.
pub fn bootstrap(
    em_info: &EMInfo,
    presence: Option<PresenceOpts>,
    num_boot: u32,
    nthreads: usize,
) -> Vec<Vec<f64>> {
    let span = span!(tracing::Level::INFO, "bootstrap_truncation");
    let _guard = span.enter();
    info!("will collect {num_boot} bootstraps");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        let f = Flat::new(em_info);
        let nreads = f.off.len() - 1;
        let mut rng = rand::rng();
        (0..num_boot)
            .map(|i| {
                let mut mult = vec![0_u32; nreads];
                for r in crate::bootstrap::get_sample_inds(nreads, &mut rng) {
                    mult[r] += 1;
                }
                let (counts, par) = run(em_info, &f, Some(&mult), presence);
                info!("bootstrap replicate {i}: {par:?}");
                counts
            })
            .collect()
    })
}

/// Layout-specific state for the atomic-free variants of [`step`] (see `crate::em_read`).
enum Layout {
    Atomic,
    /// Read chunks and one abundance buffer per chunk.
    Local(Vec<(usize, usize)>, Vec<Vec<f64>>),
    /// Per-alignment posteriors in read order, and the transcript-major transpose: transcript
    /// `t` owns alignment indices `off[t]..off[t + 1]` of `aln`.
    Csr {
        post: Vec<f64>,
        off: Vec<usize>,
        aln: Vec<u32>,
    },
    /// `Csr` once the parameters are frozen, as in `crate::em_read`: transcript `t` owns
    /// slots `off[t]..off[t + 1]` of (`read`, `w`), each an alignment's read and fixed weight,
    /// and `inv` holds each read's multiplicity over its normalizer.
    CsrFrozen {
        off: Vec<usize>,
        read: Vec<u32>,
        w: Vec<f64>,
        inv: Vec<f64>,
    },
}

impl Layout {
    fn new(f: &Flat, ntx: usize) -> Self {
        Layout::with(f, ntx, crate::em_read::Impl::from_env())
    }

    fn with(f: &Flat, ntx: usize, how: Option<crate::em_read::Impl>) -> Self {
        let nthreads = rayon::current_num_threads();
        match how {
            None => Layout::Atomic,
            Some(crate::em_read::Impl::Local) => {
                let n = (nthreads * 2).max(1);
                let per = f.target.len().div_ceil(n).max(1);
                let nreads = f.off.len() - 1;
                let (mut ch, mut start) = (Vec::new(), 0);
                for r in 0..nreads {
                    if f.off[r + 1] - f.off[start] >= per {
                        ch.push((start, r + 1));
                        start = r + 1;
                    }
                }
                if start < nreads {
                    ch.push((start, nreads));
                }
                let bufs = (0..ch.len()).map(|_| vec![0.0; ntx]).collect();
                Layout::Local(ch, bufs)
            }
            Some(crate::em_read::Impl::Csr) => {
                let mut off = vec![0usize; ntx + 1];
                for &t in &f.target {
                    off[t as usize + 1] += 1;
                }
                for i in 0..ntx {
                    off[i + 1] += off[i];
                }
                let mut fill = off.clone();
                let mut aln = vec![0u32; f.target.len()];
                for (k, &t) in f.target.iter().enumerate() {
                    aln[fill[t as usize]] = k as u32;
                    fill[t as usize] += 1;
                }
                Layout::Csr {
                    post: vec![0.0; f.target.len()],
                    off,
                    aln,
                }
            }
        }
    }
}

impl Layout {
    /// The layout for frozen weights `ws`: `Csr` becomes `CsrFrozen`, the others stay. The
    /// per-alignment posteriors and indices are freed first and the transpose is refilled
    /// from read order, so the two layouts do not coexist.
    fn frozen(self, f: &Flat, ws: &[f64]) -> Self {
        match self {
            Layout::Csr { off, post, aln } => {
                drop((post, aln));
                let nreads = f.off.len() - 1;
                let mut fill = off.clone();
                let mut read = vec![0u32; f.target.len()];
                let mut w = vec![0f64; f.target.len()];
                for (r, win) in f.off.windows(2).enumerate() {
                    let (lo, hi) = (win[0], win[1]);
                    for (&t, &wk) in f.target[lo..hi].iter().zip(&ws[lo..hi]) {
                        let t = t as usize;
                        read[fill[t]] = r as u32;
                        w[fill[t]] = wk;
                        fill[t] += 1;
                    }
                }
                Layout::CsrFrozen {
                    off,
                    read,
                    w,
                    inv: vec![0.0; nreads],
                }
            }
            lay => lay,
        }
    }
}

/// One EM round in the chosen layout. The atomic-free layouts write each `curr` slot once.
fn step_with(
    lay: &mut Layout,
    f: &Flat,
    par: &Params,
    prev: &[AtomicF64],
    curr: &[AtomicF64],
    mult: Option<&[u32]>,
    e: &Exps,
) -> Stats {
    // Posteriors of read r's alignments, passed to `emit(k, t, post)`; nothing if the read
    // cannot be assigned. Each alignment's positional terms are computed once, into `scratch`,
    // and reused for the posterior.
    let read = |r: usize,
                st: &mut Stats,
                scratch: &mut Vec<(Terms, f64)>,
                emit: &mut dyn FnMut(usize, usize, f64)| {
        let m = mult.map_or(1.0, |m| m[r] as f64);
        if m == 0.0 {
            return;
        }
        let (lo, hi) = (f.off[r], f.off[r + 1]);
        scratch.clear();
        let mut denom = 0.0;
        for k in lo..hi {
            let tid = f.target[k] as usize;
            let t = par.terms(f, k, e);
            let w = prev[tid].load(Ordering::Relaxed) * f.prob[k] as f64 * t.total();
            denom += w;
            scratch.push((t, w));
        }
        if denom > constants::EM_DENOM_THRESH {
            for (i, k) in (lo..hi).enumerate() {
                let (t, w) = &scratch[i];
                let post = m * w / denom;
                emit(k, f.target[k] as usize, post);
                st.add(f, k, t, post);
            }
        }
    };
    match lay {
        Layout::Atomic => step(f, par, e, prev, curr, mult),
        Layout::CsrFrozen { .. } => unreachable!("only frozen weights use this layout"),
        Layout::Local(ch, bufs) => {
            let st = ch
                .par_iter()
                .zip(bufs.par_iter_mut())
                .map(|(&(r0, r1), acc)| {
                    acc.fill(0.0);
                    let (mut st, mut scratch) = (Stats::default(), Vec::new());
                    for r in r0..r1 {
                        read(r, &mut st, &mut scratch, &mut |_, t, p| acc[t] += p);
                    }
                    st
                })
                .reduce(Stats::default, Stats::merge);
            const BLOCK: usize = 4096;
            let bufs = &*bufs;
            curr.par_chunks(BLOCK).enumerate().for_each(|(b, out)| {
                let base = b * BLOCK;
                for (i, o) in out.iter().enumerate() {
                    o.store(
                        bufs.iter().map(|acc| acc[base + i]).sum(),
                        Ordering::Relaxed,
                    );
                }
            });
            st
        }
        Layout::Csr { post, off, aln } => {
            // Pass 1, over reads: posteriors into the read-ordered buffer (each read writes
            // only its own alignments' slots), and the sufficient statistics.
            let nreads = f.off.len() - 1;
            let post_ptr = post.as_mut_ptr() as usize;
            let st = (0..nreads)
                .into_par_iter()
                .fold(
                    || (Stats::default(), Vec::new()),
                    |(mut st, mut scratch), r| {
                        for k in f.off[r]..f.off[r + 1] {
                            // SAFETY: slot k belongs to read r alone.
                            unsafe { *(post_ptr as *mut f64).add(k) = 0.0 };
                        }
                        read(r, &mut st, &mut scratch, &mut |k, _, p| unsafe {
                            *(post_ptr as *mut f64).add(k) = p
                        });
                        (st, scratch)
                    },
                )
                .map(|(st, _)| st)
                .reduce(Stats::default, Stats::merge);
            // Pass 2, over transcripts: sum their alignments' posteriors.
            let post = &*post;
            curr.par_iter().enumerate().for_each(|(t, c)| {
                let s: f64 = aln[off[t]..off[t + 1]]
                    .iter()
                    .map(|&k| post[k as usize])
                    .sum();
                c.store(s, Ordering::Relaxed);
            });
            st
        }
    }
}

/// One EM round with frozen parameters: `ws[k]` is alignment k's fixed weight (score
/// probability times its positional terms), so no terms or statistics are computed.
fn step_static(
    lay: &mut Layout,
    f: &Flat,
    ws: &[f64],
    prev: &[AtomicF64],
    curr: &[AtomicF64],
    mult: Option<&[u32]>,
) {
    let read = |r: usize, emit: &mut dyn FnMut(usize, usize, f64)| {
        let m = mult.map_or(1.0, |m| m[r] as f64);
        if m == 0.0 {
            return;
        }
        let (lo, hi) = (f.off[r], f.off[r + 1]);
        let mut denom = 0.0;
        for (&t, &w) in f.target[lo..hi].iter().zip(&ws[lo..hi]) {
            denom += prev[t as usize].load(Ordering::Relaxed) * w;
        }
        if denom > constants::EM_DENOM_THRESH {
            let inv = m / denom;
            for (k, (&t, &w)) in (lo..hi).zip(f.target[lo..hi].iter().zip(&ws[lo..hi])) {
                let t = t as usize;
                emit(k, t, prev[t].load(Ordering::Relaxed) * w * inv);
            }
        }
    };
    match lay {
        Layout::Atomic => unreachable!("parameters are not frozen in the atomic layout"),
        Layout::CsrFrozen { off, read, w, inv } => {
            inv.par_iter_mut().enumerate().for_each(|(r, x)| {
                let m = mult.map_or(1.0, |m| m[r] as f64);
                let (lo, hi) = (f.off[r], f.off[r + 1]);
                let mut denom = 0.0;
                for (&t, &w) in f.target[lo..hi].iter().zip(&ws[lo..hi]) {
                    denom += prev[t as usize].load(Ordering::Relaxed) * w;
                }
                *x = if m > 0.0 && denom > constants::EM_DENOM_THRESH {
                    m / denom
                } else {
                    0.0
                };
            });
            let inv = &*inv;
            curr.par_iter().enumerate().for_each(|(t, c)| {
                let s: f64 = (off[t]..off[t + 1])
                    .map(|i| w[i] * inv[read[i] as usize])
                    .sum();
                c.store(prev[t].load(Ordering::Relaxed) * s, Ordering::Relaxed);
            });
        }
        Layout::Local(ch, bufs) => {
            ch.par_iter()
                .zip(bufs.par_iter_mut())
                .for_each(|(&(r0, r1), acc)| {
                    acc.fill(0.0);
                    for r in r0..r1 {
                        read(r, &mut |_, t, p| acc[t] += p);
                    }
                });
            const BLOCK: usize = 4096;
            let bufs = &*bufs;
            curr.par_chunks(BLOCK).enumerate().for_each(|(b, out)| {
                let base = b * BLOCK;
                for (i, o) in out.iter().enumerate() {
                    o.store(
                        bufs.iter().map(|acc| acc[base + i]).sum(),
                        Ordering::Relaxed,
                    );
                }
            });
        }
        Layout::Csr { post, off, aln } => {
            let nreads = f.off.len() - 1;
            let post_ptr = post.as_mut_ptr() as usize;
            (0..nreads).into_par_iter().for_each(|r| {
                for k in f.off[r]..f.off[r + 1] {
                    // SAFETY: slot k belongs to read r alone.
                    unsafe { *(post_ptr as *mut f64).add(k) = 0.0 };
                }
                read(r, &mut |k, _, p| unsafe {
                    *(post_ptr as *mut f64).add(k) = p
                });
            });
            let post = &*post;
            curr.par_iter().enumerate().for_each(|(t, c)| {
                let s: f64 = aln[off[t]..off[t + 1]]
                    .iter()
                    .map(|&k| post[k as usize])
                    .sum();
                c.store(s, Ordering::Relaxed);
            });
        }
    }
}

/// Parameters stop being re-estimated once none moves by more than this (see [`moved`]) for
/// `FREEZE_AFTER` consecutive iterations; the remaining iterations use fixed weights. Updating
/// the parameters less often is still a generalized EM; on K562 PacBio (truncation model)
/// freezing at 1e-6 moves no transcript estimate by more than 0.07% and cuts EM time 3.6x.
/// Only the atomic-free layouts freeze.
const FREEZE_TOL: f64 = 1e-6;
const FREEZE_AFTER: u32 = 3;

/// How far the parameters moved: absolute change for the mixture weights, relative change for
/// the rates. (A fraction decaying toward its floor changes by large relative amounts long after
/// its effect on the likelihood is negligible.)
fn moved(a: &Params, b: &Params) -> f64 {
    let rel = |x: f64, y: f64| (x - y).abs() / y.abs().max(1e-12);
    (a.p3 - b.p3)
        .abs()
        .max((a.p5 - b.p5).abs())
        .max((a.pf - b.pf).abs())
        .max((a.pu - b.pu).abs())
        .max(rel(a.h3, b.h3))
        .max(rel(a.h5, b.h5))
        .max(rel(a.lam, b.lam))
}

/// Refresh the presence posteriors if due at iteration `niter`, then fill `eff` with
/// `q * prev`. Each alignment's weight is its score probability times its positional
/// likelihood: the frozen weights `ws` if the parameters are frozen, else from `par`.
#[allow(clippy::too_many_arguments)]
fn refresh_presence(
    p: &mut Presence,
    niter: u32,
    f: &Flat,
    par: &Params,
    ws: Option<&[f64]>,
    prev: &[AtomicF64],
    eff: &[AtomicF64],
    mult: Option<&[u32]>,
    max_len: u32,
) {
    let counts: Vec<f64> = prev.iter().map(|x| x.load(Ordering::Relaxed)).collect();
    if p.due(niter) {
        let e = ws.is_none().then(|| Exps::new(par, max_len));
        p.update(&counts, f.off.len() - 1, |r, buf| {
            for k in f.off[r]..f.off[r + 1] {
                let w = match (ws, &e) {
                    (Some(ws), _) => ws[k],
                    (None, Some(e)) => f.prob[k] as f64 * par.terms(f, k, e).total(),
                    (None, None) => unreachable!(),
                };
                buf.push((f.target[k], w));
            }
            mult.map_or(1.0, |m| m[r] as f64)
        });
    }
    eff.par_iter()
        .zip(counts.par_iter().zip(&p.q))
        .for_each(|(e, (c, q))| e.store(q * c, Ordering::Relaxed));
}

fn run(
    em_info: &EMInfo,
    f: &Flat,
    mult: Option<&[u32]>,
    presence: Option<PresenceOpts>,
) -> (Vec<f64>, Params) {
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let total_weight = em_info.eq_map.num_aligned_reads() as f64;
    let init = match em_info.init_abundances {
        Some(ref v) => v.clone(),
        None => vec![total_weight / (tinfo.len() as f64); tinfo.len()],
    };
    let mut prev: Vec<AtomicF64> = init.into_iter().map(AtomicF64::new).collect();
    let mut curr: Vec<AtomicF64> = (0..tinfo.len()).map(|_| AtomicF64::new(0.0)).collect();
    let t_setup = std::time::Instant::now();
    let mut lay = Layout::new(f, tinfo.len());
    info!(
        "truncation EM layout set up in {:.2}s",
        t_setup.elapsed().as_secs_f64()
    );
    let mut par = Params::initial();

    let mut niter = 0_u32;
    let max_len = f.max_len();
    let can_freeze = !matches!(lay, Layout::Atomic);
    let (mut still, mut frozen): (u32, Option<Vec<f64>>) = (0, None);
    // With the presence model the EM sees `q * prev` (see `crate::presence`); `eff` holds it.
    let mut pres = presence.map(|o| Presence::new(tinfo.len(), o));
    let eff: Vec<AtomicF64> = if pres.is_some() {
        (0..tinfo.len()).map(|_| AtomicF64::new(0.0)).collect()
    } else {
        Vec::new()
    };
    while niter < em_info.max_iter {
        if let Some(p) = &mut pres {
            refresh_presence(
                p,
                niter,
                f,
                &par,
                frozen.as_deref(),
                &prev,
                &eff,
                mult,
                max_len,
            );
        }
        let src = if pres.is_some() { &eff } else { &prev };
        match &frozen {
            Some(ws) => step_static(&mut lay, f, ws, src, &curr, mult),
            None => {
                let e = Exps::new(&par, max_len);
                let next = step_with(&mut lay, f, &par, src, &curr, mult, &e).params(par);
                still = if moved(&next, &par) < FREEZE_TOL {
                    still + 1
                } else {
                    0
                };
                par = next;
                if can_freeze && still >= FREEZE_AFTER {
                    let e = Exps::new(&par, max_len);
                    let ws: Vec<f64> = (0..f.target.len())
                        .into_par_iter()
                        .map(|k| {
                            let t = par.terms(f, k, &e);
                            f.prob[k] as f64 * t.total()
                        })
                        .collect();
                    info!(
                        "truncation parameters frozen after {} iterations: {:?}",
                        niter + 1,
                        par
                    );
                    lay = std::mem::replace(&mut lay, Layout::Atomic).frozen(f, &ws);
                    frozen = Some(ws);
                }
            }
        }
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
    // Zero very small abundances, then one more round with the parameters fixed, as the
    // per-read EM does.
    for x in &prev {
        if x.load(Ordering::Relaxed) < constants::MIN_READ_THRESH {
            x.store(0.0, Ordering::Relaxed);
        }
    }
    if let Some(p) = &pres {
        if mult.is_none() {
            info!(
                "presence model: {} transcripts called absent (q < 0.5)",
                p.n_absent().to_formatted_string(&Locale::en)
            );
        }
        let counts: Vec<f64> = prev.iter().map(|x| x.load(Ordering::Relaxed)).collect();
        for (e, (c, q)) in eff.iter().zip(counts.iter().zip(&p.q)) {
            e.store(q * c, Ordering::Relaxed);
        }
    }
    let src = if pres.is_some() { &eff } else { &prev };
    match &frozen {
        Some(ws) => step_static(&mut lay, f, ws, src, &curr, mult),
        None => {
            step_with(
                &mut lay,
                f,
                &par,
                src,
                &curr,
                mult,
                &Exps::new(&par, max_len),
            );
        }
    }
    (
        curr.iter().map(|x| x.load(Ordering::Relaxed)).collect(),
        par,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat(alns: &[Vec<(u32, f64, f64, f64)>], lens: &[f64]) -> Flat {
        let mut f = Flat::empty(lens);
        for r in alns {
            for &(t, p, s, e) in r {
                f.push(t, p as f32, s, e);
            }
            f.off.push(f.target.len());
        }
        f
    }

    fn fit(f: &Flat, ntx: usize, iters: usize) -> (Vec<f64>, Params) {
        fit_mult(f, ntx, iters, None)
    }

    fn fit_mult(f: &Flat, ntx: usize, iters: usize, mult: Option<&[u32]>) -> (Vec<f64>, Params) {
        fit_layout(f, ntx, iters, mult, None)
    }

    fn fit_layout(
        f: &Flat,
        ntx: usize,
        iters: usize,
        mult: Option<&[u32]>,
        how: Option<crate::em_read::Impl>,
    ) -> (Vec<f64>, Params) {
        let mut lay = Layout::with(f, ntx, how);
        let n = mult.map_or((f.off.len() - 1) as f64, |m| m.iter().sum::<u32>() as f64);
        let mut prev: Vec<AtomicF64> = (0..ntx).map(|_| AtomicF64::new(n / ntx as f64)).collect();
        let mut curr: Vec<AtomicF64> = (0..ntx).map(|_| AtomicF64::new(0.0)).collect();
        let mut par = Params::initial();
        for _ in 0..iters {
            let e = Exps::new(&par, f.max_len());
            par = step_with(&mut lay, f, &par, &prev, &curr, mult, &e).params(par);
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
        // With exactly placed ends the mixture weights are not identifiable here (uniform
        // placement explains a full-length read as well as anchoring does), so only the
        // estimates are checked.
        assert!(est[1] > 500.0, "{est:?} {par:?}");
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

    /// A bootstrap multiplicity is the same as repeating the read: counting reads twice or
    /// dropping them through `mult` matches the fit on the resampled reads written out.
    #[test]
    fn multiplicities_match_repeated_reads() {
        let lens = [2000.0, 1000.0];
        let reads: Vec<Vec<(u32, f64, f64, f64)>> = (0..60)
            .map(|i| match i % 3 {
                0 => vec![(0, 1.0, 0.0, 2000.0)],
                1 => vec![(0, 1.0, 300.0 + 10.0 * i as f64, 1990.0)],
                _ => vec![(0, 1.0, 1000.0, 2000.0), (1, 0.9, 0.0, 1000.0)],
            })
            .collect();
        let mult: Vec<u32> = (0..reads.len()).map(|i| (i % 4) as u32).collect();
        let repeated: Vec<_> = reads
            .iter()
            .zip(&mult)
            .flat_map(|(r, &m)| std::iter::repeat_n(r.clone(), m as usize))
            .collect();
        let (a, pa) = fit_mult(&flat(&reads, &lens), 2, 100, Some(&mult));
        let (b, pb) = fit(&flat(&repeated, &lens), 2, 100);
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() < 1e-6 * y.max(1.0), "{a:?} {b:?}");
        }
        assert!(
            (pa.h3 - pb.h3).abs() < 1e-9 && (pa.p3 - pb.p3).abs() < 1e-9,
            "{pa:?} {pb:?}"
        );
    }

    /// The atomic-free layouts give the atomic EM's estimates and parameters.
    #[test]
    fn layouts_agree() {
        use crate::em_read::Impl;
        let lens = [2000.0, 1000.0];
        let mut reads = Vec::new();
        for i in 0..300 {
            reads.push(vec![(0, 1.0, 0.0, 2000.0)]);
            reads.push(vec![(0, 1.0, 200.0 + (i % 50) as f64 * 6.0, 2000.0)]);
            reads.push(vec![(0, 1.0, 1000.0, 2000.0), (1, 0.9, 0.0, 1000.0)]);
        }
        let f = flat(&reads, &lens);
        let mult: Vec<u32> = (0..reads.len()).map(|i| (i % 3) as u32).collect();
        for m in [None, Some(&mult[..])] {
            let (a, pa) = fit_layout(&f, 2, 60, m, None);
            for how in [Impl::Local, Impl::Csr] {
                let (b, pb) = fit_layout(&f, 2, 60, m, Some(how));
                for (x, y) in a.iter().zip(&b) {
                    assert!((x - y).abs() < 1e-9 * x.max(1.0), "{how:?} {a:?} {b:?}");
                }
                assert!(
                    (pa.h3 - pb.h3).abs() < 1e-12 && (pa.p3 - pb.p3).abs() < 1e-12,
                    "{how:?}"
                );
            }
        }
    }

    /// With frozen weights, the thread-local and transposed layouts give the same round, with
    /// and without bootstrap multiplicities; reordering the reads does not change the fit.
    #[test]
    fn frozen_layouts_and_read_order_agree() {
        use crate::em_read::Impl;
        let lens = [2000.0, 1000.0, 1500.0];
        let mut reads = Vec::new();
        for i in 0..200 {
            reads.push(vec![(2, 1.0, 0.0, 1500.0)]);
            reads.push(vec![(0, 1.0, 200.0 + (i % 50) as f64 * 6.0, 2000.0)]);
            reads.push(vec![(0, 1.0, 1000.0, 2000.0), (1, 0.9, 0.0, 1000.0)]);
            reads.push(vec![(1, 1.0, 400.0, 1000.0), (2, 0.8, 900.0, 1500.0)]);
        }
        let f = flat(&reads, &lens);
        let par = Params::initial();
        let e = Exps::new(&par, f.max_len());
        let ws: Vec<f64> = (0..f.target.len())
            .map(|k| {
                let t = par.terms(&f, k, &e);
                f.prob[k] as f64 * t.total()
            })
            .collect();
        let prev: Vec<AtomicF64> = [3.0, 1.0, 2.0].into_iter().map(AtomicF64::new).collect();
        let mult: Vec<u32> = (0..reads.len()).map(|i| (i % 3) as u32).collect();
        for m in [None, Some(&mult[..])] {
            let round = |how| {
                let mut lay = Layout::with(&f, 3, Some(how)).frozen(&f, &ws);
                let curr: Vec<AtomicF64> = (0..3).map(|_| AtomicF64::new(0.0)).collect();
                step_static(&mut lay, &f, &ws, &prev, &curr, m);
                curr.iter()
                    .map(|x| x.load(Ordering::Relaxed))
                    .collect::<Vec<_>>()
            };
            let (a, b) = (round(Impl::Local), round(Impl::Csr));
            for (x, y) in a.iter().zip(&b) {
                assert!((x - y).abs() < 1e-9 * x.max(1.0), "{a:?} {b:?}");
            }
        }
        let (a, pa) = fit(&f, 3, 50);
        let rev: Vec<_> = reads.iter().rev().cloned().collect();
        let (b, pb) = fit(&flat(&rev, &lens), 3, 50);
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() < 1e-9 * x.max(1.0), "{a:?} {b:?}");
        }
        assert!((pa.h3 - pb.h3).abs() < 1e-12 && (pa.lam - pb.lam).abs() < 1e-12);
    }
}
