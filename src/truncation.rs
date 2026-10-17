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
//! `l = e - s` and gaps `d5 = s`, `d3 = L - e`: each end has a kernel `K3(d3)`, `K5(d5)`, the
//! distribution of a read's gap to that end when the read sits at it (a mixture of a short- and
//! a long-scale geometric law, learned: it absorbs both how tightly molecules are anchored and
//! how many bases the aligner trims). Every component draws a read length from its own law,
//! censored at the length available to it (a read that would be longer runs to the far
//! end), and places the read by its own rule. With `f_c`, `S_c` the density and survival of
//! component `c`'s read-length law (a histogram, see [`LenLaw`]):
//!
//! ```text
//! 3'-anchored  K3(d3) ( f_3(l)  +  S_3(L - d3) K5(d5) )
//! 5'-anchored  K5(d5) ( f_5(l)  +  S_5(L - d5) K3(d3) )
//! full length  K3(d3) K5(d5)
//! uniform      f_u(l) / (L - l + 1)  +  S_u(L) K3(d3) K5(d5)
//! P(s, e | t) = p3 · 3'-anchored + p5 · 5'-anchored + pf · full length + pu · uniform
//! ```
//!
//! Each component sums to (about) one over the placements a transcript of length `L` allows, so
//! a transcript is never favored merely because one of its ends happens to coincide with a
//! read's end: under the uniform component every placement of a read of length `l` is equally
//! likely, and a candidate is judged only by how much room it gave the read. When the reads are
//! not anchored (`pu` near one) the model reduces to uniform placement; when they are (`p3` near
//! one) it is a 3'-anchored truncation model. The parameters are fitted in the EM together with
//! the abundances (mixture weights from the components' posterior shares, rates by
//! censored-histogram M-steps, kernels by geometric-mixture M-steps); each alignment's
//! likelihood is its alignment-score probability times `P(s, e | t)`.

use std::sync::atomic::Ordering;

use atomic_float::AtomicF64;
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use tracing::{info, span, trace};

use crate::presence::{Presence, PresenceOpts};
use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use std::borrow::Cow;

/// The alignment store as the parameter-estimating EM sees it: reads in store order, read `r`
/// owning alignments `off[r]..off[r + 1]` (positions in the store). Nothing per alignment is
/// copied: targets, score probabilities and spans are read from the store, which those passes
/// stream through in order.
struct Flat<'a> {
    off: Vec<usize>,
    alns: Cow<'a, [AlnInfo]>,
    probs: Cow<'a, [f32]>,
    /// Transcript lengths, per transcript.
    tlen: Vec<u32>,
}

impl<'a> Flat<'a> {
    #[cfg(test)]
    fn push(&mut self, t: u32, p: f32, start: f64, end: f64) {
        self.alns.to_mut().push(AlnInfo {
            ref_id: t,
            start: start.max(0.0).round() as u32,
            end: end.max(0.0).round() as u32,
            strand: bio_types::strand::Strand::Forward,
        });
        self.probs.to_mut().push(p);
    }

    /// Number of alignments.
    #[inline]
    fn n(&self) -> usize {
        self.alns.len()
    }

    #[inline]
    fn target(&self, k: usize) -> u32 {
        self.alns[k].ref_id
    }

    /// The geometry of alignment `k`.
    #[inline]
    fn geo(&self, k: usize) -> Geo {
        let a = &self.alns[k];
        let len = self.tlen[a.ref_id as usize];
        let start = a.start.min(len.saturating_sub(1));
        let end = a.end.clamp(start, len);
        Geo {
            t: a.ref_id,
            p: self.probs[k],
            len,
            l: end - start,
            d5: start,
            d3: len - end,
        }
    }

    /// The geometries of read `r`'s alignments.
    #[inline]
    fn read_geos(&self, r: usize) -> impl Iterator<Item = Geo> + '_ {
        (self.off[r]..self.off[r + 1]).map(move |k| self.geo(k))
    }

    fn max_len(&self) -> u32 {
        self.tlen.iter().copied().max().unwrap_or(0)
    }

    fn new(em_info: &'a EMInfo) -> Self {
        let store = &em_info.eq_map;
        let n = store.alignments.len();
        assert!(
            n <= u32::MAX as usize,
            "the truncation model indexes alignments with 32 bits; {n} is too many"
        );
        let mut off = Vec::with_capacity(store.len() + 1);
        off.push(0);
        off.extend((0..store.len()).map(|r| store.read_range(r).end));
        Flat {
            off,
            alns: Cow::Borrowed(&store.alignments),
            probs: Cow::Borrowed(&store.as_probabilities),
            tlen: em_info
                .txp_info
                .iter()
                .map(|t| t.lenf.round() as u32)
                .collect(),
        }
    }

    /// Every `k`-th read (from the first), as a table of its own. Reads are in store order, by
    /// likeliest target, so the sample spreads over the transcripts in proportion to their reads.
    fn every(&self, k: usize) -> Flat<'static> {
        let nreads = self.off.len() - 1;
        let picked = (0..nreads).step_by(k.max(1));
        let (mut off, mut alns, mut probs) = (vec![0usize], Vec::new(), Vec::new());
        for r in picked {
            let (lo, hi) = (self.off[r], self.off[r + 1]);
            alns.extend_from_slice(&self.alns[lo..hi]);
            probs.extend_from_slice(&self.probs[lo..hi]);
            off.push(alns.len());
        }
        Flat {
            off,
            alns: Cow::Owned(alns),
            probs: Cow::Owned(probs),
            tlen: self.tlen.clone(),
        }
    }

    #[cfg(test)]
    fn empty(lens: &[f64]) -> Self {
        Flat {
            off: vec![0],
            alns: Cow::Owned(vec![]),
            probs: Cow::Owned(vec![]),
            tlen: lens.iter().map(|l| l.round() as u32).collect(),
        }
    }
}

/// What the positional model needs of one alignment: its target and score probability, the
/// target's length, the aligned length, and the gaps to the transcript's 5' and 3' ends (the
/// span clamped to the transcript).
#[derive(Clone, Copy)]
struct Geo {
    t: u32,
    p: f32,
    len: u32,
    l: u32,
    d5: u32,
    d3: u32,
}

/// The frozen-parameter EM's inputs: each alignment's target and fixed weight (score
/// probability times positional likelihood), in [`Flat`]'s order, whose read boundaries they
/// share: entry `k` is alignment `k`, and read `r` owns entries `off[r]..off[r + 1]`. That is
/// store order, reads by likeliest target (`bulk` orders the store so), which the
/// transcript-major pass reads from nearby memory.
struct Fixed<'f> {
    off: &'f [usize],
    t: Vec<u32>,
    w: Vec<f64>,
}

impl<'f> Fixed<'f> {
    fn new(f: &'f Flat, par: &Params) -> Self {
        let e = Exps::new(par, f.max_len());
        let (t, w) = (0..f.n())
            .into_par_iter()
            .map(|k| {
                let g = f.geo(k);
                (g.t, g.p as f64 * par.terms(&g, &e).total())
            })
            .unzip();
        Fixed { off: &f.off, t, w }
    }

    fn nreads(&self) -> usize {
        self.off.len() - 1
    }

    /// Multiplicity of read `r`.
    #[inline]
    fn mult(&self, mult: Option<&[u32]>, r: usize) -> f64 {
        mult.map_or(1.0, |m| m[r] as f64)
    }
}

/// An end kernel: the distribution of the gap `d` (nt) between a read's end and the transcript's
/// end, for reads placed at that end. Nonparametric, on log-spaced bins `{0}`, `{1}`, `[2, 3]`,
/// `[4, 7]`, ... (uniform within a bin), so it can represent both reads flush with the end and
/// ends spread hundreds of bases from it (alternative polyadenylation, annotation mismatch).
/// For a transcript of length `L` it is normalized over the gaps the transcript allows,
/// `0..L`, and fitted by the truncated-data EM.
#[derive(Clone)]
struct GapLaw {
    h: Vec<f64>,
}

impl std::fmt::Debug for GapLaw {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.summary();
        write!(
            f,
            "GapLaw(0: {:.2}, 1-15: {:.2}, 16-255: {:.2}, 256+: {:.2})",
            s[0], s[1], s[2], s[3]
        )
    }
}

impl GapLaw {
    #[inline]
    fn bin(d: u32) -> usize {
        (32 - d.leading_zeros()) as usize
    }

    #[inline]
    fn width(b: usize) -> f64 {
        if b == 0 {
            1.0
        } else {
            (1u64 << (b - 1)) as f64
        }
    }

    /// Discretized exponential law of mean 100 nt over gaps up to `max_len`.
    fn initial(max_len: u32) -> Self {
        let nb = Self::bin(max_len) + 1;
        let cdf = |x: f64| 1.0 - (-x / 100.0).exp();
        let h: Vec<f64> = (0..nb)
            .map(|b| {
                let lo = if b == 0 { 0.0 } else { Self::width(b) };
                cdf(lo + Self::width(b)) - cdf(lo)
            })
            .collect();
        let s: f64 = h.iter().sum();
        GapLaw {
            h: h.into_iter().map(|x| x / s).collect(),
        }
    }

    /// Per-nt density `dens[d]` and distribution function `cdf[x] = P(gap <= x)` for gaps up
    /// to `max_len`.
    fn tables(&self, max_len: u32) -> (Vec<f64>, Vec<f64>) {
        let dens: Vec<f64> = (0..=max_len)
            .map(|d| {
                let b = Self::bin(d);
                self.h.get(b).map_or(0.0, |p| p / Self::width(b))
            })
            .collect();
        let mut acc = 0.0;
        let cdf = dens
            .iter()
            .map(|x| {
                acc += x;
                acc
            })
            .collect();
        (dens, cdf)
    }

    /// Mass at gap 0, in 1-15, 16-255 and 256 or more.
    fn summary(&self) -> [f64; 4] {
        let mut s = [0.0; 4];
        for (b, p) in self.h.iter().enumerate() {
            s[match b {
                0 => 0,
                1..=4 => 1,
                5..=8 => 2,
                _ => 3,
            }] += p;
        }
        s
    }

    fn moved(&self, o: &GapLaw) -> f64 {
        0.5 * self
            .h
            .iter()
            .zip(&o.h)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
    }

    /// Truncated-data EM step from `old`: `s.stop[b]` is the observed weight in bin `b`,
    /// `s.cens[c]` the weight of observations truncated at a gap in bin `c` (the transcript
    /// allowed no larger gap); each implies unobserved draws beyond its truncation point,
    /// spread over the bins above `c` in proportion to `old.h`.
    fn law(s: &LStats, old: &GapLaw) -> GapLaw {
        let nb = old.h.len();
        let mut n = vec![0.0; nb];
        let (mut f, mut acc) = (0.0, 0.0);
        for (b, (nb_, &hb)) in n.iter_mut().zip(&old.h).enumerate() {
            *nb_ = s.stop.get(b).copied().unwrap_or(0.0) + acc * hb;
            f += hb;
            let c = s.cens.get(b).copied().unwrap_or(0.0);
            if c > 0.0 && f > 0.0 {
                // (1 - F) / F unobserved draws per truncated one, spread as h_b / (1 - F)
                acc += c / f;
            }
        }
        let tot: f64 = n.iter().sum();
        if tot <= 0.0 {
            return old.clone();
        }
        let eps = 1e-9 * tot;
        let z = tot + eps * nb as f64;
        GapLaw {
            h: n.into_iter().map(|x| (x + eps) / z).collect(),
        }
    }
}

/// Bins of the relative part of an end kernel (gap as a fraction of transcript length).
const RB: usize = 20;

/// An end kernel: how far a read placed at an end sits from it. A mixture of an absolute part
/// (a [`GapLaw`] in nucleotides, for reads flush with the end or a fixed distance from it) and a
/// relative part (a histogram over the gap as a fraction of the transcript's length, for ends
/// that fall at a relative position, as poly(A) sites in long 3' UTRs do), each a distribution
/// over the gaps `0..L` a transcript of length `L` allows: `K(d | L) = (1 - rho) K_abs(d) /
/// C_abs(L - 1) + rho g(d / L) / L`.
#[derive(Clone)]
struct EndLaw {
    abs: GapLaw,
    rel: Vec<f64>,
    rho: f64,
}

impl std::fmt::Debug for EndLaw {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let near: f64 = self.rel[..RB / 10].iter().sum();
        let far: f64 = self.rel[RB / 2..].iter().sum();
        write!(
            f,
            "{:?} + relative {:.2} (gap/L < 0.1: {:.2}, >= 0.5: {:.2})",
            self.abs, self.rho, near, far
        )
    }
}

impl EndLaw {
    fn initial(max_len: u32) -> Self {
        EndLaw {
            abs: GapLaw::initial(max_len),
            rel: vec![1.0 / RB as f64; RB],
            rho: 0.5,
        }
    }

    #[inline]
    fn rbin(d: u32, len: u32) -> usize {
        ((d as u64 * RB as u64) / len.max(1) as u64).min(RB as u64 - 1) as usize
    }

    /// The absolute and relative parts of `K(d | len)`, given the absolute part's tables.
    #[inline]
    fn parts(&self, d: u32, len: u32, dens: &[f64], cdf: &[f64]) -> (f64, f64) {
        let last = len.saturating_sub(1) as usize;
        let a = (1.0 - self.rho) * dens[d as usize] / cdf[last].max(1e-300);
        let r = self.rho * self.rel[Self::rbin(d, len)] * RB as f64 / len.max(1) as f64;
        (a, r)
    }

    fn moved(&self, o: &EndLaw) -> f64 {
        let tv: f64 = 0.5
            * self
                .rel
                .iter()
                .zip(&o.rel)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>();
        self.abs.moved(&o.abs).max(tv).max((self.rho - o.rho).abs())
    }

    fn law(s: &EStats, old: &EndLaw) -> EndLaw {
        let n = s.n_abs + s.n_rel;
        if n <= 0.0 {
            return old.clone();
        }
        let tot: f64 = s.rel.iter().sum();
        let eps = 1e-9 * tot.max(1e-300);
        let z = tot + eps * RB as f64;
        EndLaw {
            abs: GapLaw::law(&s.abs, &old.abs),
            rel: if tot > 0.0 {
                (0..RB)
                    .map(|b| (s.rel.get(b).copied().unwrap_or(0.0) + eps) / z)
                    .collect()
            } else {
                old.rel.clone()
            },
            rho: (s.n_rel / n).clamp(1e-6, 1.0 - 1e-6),
        }
    }
}

/// Expected sufficient statistics of an end kernel.
#[derive(Default, Clone)]
struct EStats {
    abs: LStats,
    rel: Vec<f64>,
    n_abs: f64,
    n_rel: f64,
}

impl EStats {
    fn merge(mut self, o: EStats) -> EStats {
        self.abs = self.abs.merge(o.abs);
        if self.rel.len() < o.rel.len() {
            self.rel.resize(o.rel.len(), 0.0);
        }
        self.rel.iter_mut().zip(o.rel).for_each(|(a, b)| *a += b);
        self.n_abs += o.n_abs;
        self.n_rel += o.n_rel;
        self
    }

    /// Weight `w` of a gap `d` on a transcript of length `len`, split between the kernel's
    /// parts in proportion to their values `a`, `r` there.
    fn add(&mut self, w: f64, d: u32, len: u32, a: f64, r: f64) {
        let z = a + r;
        if w <= 0.0 || z <= 0.0 {
            return;
        }
        let (wa, wr) = (w * a / z, w * r / z);
        LStats::add(&mut self.abs.stop, GapLaw::bin(d), wa);
        LStats::add(&mut self.abs.cens, GapLaw::bin(len.saturating_sub(1)), wa);
        LStats::add(&mut self.rel, EndLaw::rbin(d, len), wr);
        self.n_abs += wa;
        self.n_rel += wr;
    }
}

/// A read-length law: a histogram over lengths in bins of [`LenLaw::W`] nt, uniform within a
/// bin. Lengths are right-censored at the length a transcript makes available (a read that
/// would be longer runs to the transcript's end), so the law is fitted by the censored-data EM
/// (each censored read is spread over the bins at or above its censoring length).
#[derive(Clone)]
struct LenLaw {
    h: Vec<f64>,
}

impl std::fmt::Debug for LenLaw {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LenLaw(mean {:.0} nt)", self.mean())
    }
}

impl LenLaw {
    const W: f64 = 50.0;

    /// Discretized exponential law with the given mean, over lengths up to `max_len`.
    fn exponential(mean: f64, max_len: u32) -> Self {
        let nb = (max_len as f64 / Self::W) as usize + 1;
        let h: Vec<f64> = (0..nb)
            .map(|b| {
                let (a, z) = (b as f64 * Self::W, (b + 1) as f64 * Self::W);
                (-a / mean).exp() - (-z / mean).exp()
            })
            .collect();
        let s: f64 = h.iter().sum();
        LenLaw {
            h: h.into_iter().map(|x| x / s).collect(),
        }
    }

    #[inline]
    fn bin(x: f64) -> usize {
        (x / Self::W) as usize
    }

    /// [`LenLaw::bin`] of a whole number of bases, by integer division (the same bin: `x / W`
    /// is at least `1 / W` from the next whole number, far beyond rounding).
    #[inline]
    fn bin_nt(x: u32) -> usize {
        const W_NT: u32 = LenLaw::W as u32;
        (x / W_NT) as usize
    }

    fn mean(&self) -> f64 {
        self.h
            .iter()
            .enumerate()
            .map(|(b, p)| p * (b as f64 + 0.5) * Self::W)
            .sum()
    }

    /// Per-nt density `dens[x]` and survival `surv[x] = P(length >= x)` for `x` up to `max_len`.
    fn tables(&self, max_len: u32) -> (Vec<f64>, Vec<f64>) {
        let nb = self.h.len();
        let mut tail = vec![0.0; nb + 1];
        for b in (0..nb).rev() {
            tail[b] = tail[b + 1] + self.h[b];
        }
        let dens = (0..=max_len)
            .map(|x| self.h.get(Self::bin(x as f64)).map_or(0.0, |p| p / Self::W))
            .collect();
        let surv = (0..=max_len)
            .map(|x| {
                let b = Self::bin(x as f64);
                if b >= nb {
                    return 0.0;
                }
                let within = (x as f64 - b as f64 * Self::W) / Self::W;
                (tail[b + 1] + self.h[b] * (1.0 - within)).max(0.0)
            })
            .collect();
        (dens, surv)
    }

    /// Total-variation distance between two laws.
    fn moved(&self, o: &LenLaw) -> f64 {
        0.5 * self
            .h
            .iter()
            .zip(&o.h)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
    }
}

/// Expected sufficient statistics of a length law: observed lengths and censoring lengths,
/// binned.
#[derive(Default, Clone)]
struct LStats {
    stop: Vec<f64>,
    cens: Vec<f64>,
}

impl LStats {
    #[inline]
    fn add(v: &mut Vec<f64>, b: usize, w: f64) {
        if w <= 0.0 {
            return;
        }
        if v.len() <= b {
            v.resize(b + 1, 0.0);
        }
        v[b] += w;
    }

    fn merge(mut self, o: LStats) -> LStats {
        for (dst, src) in [(&mut self.stop, o.stop), (&mut self.cens, o.cens)] {
            if dst.len() < src.len() {
                dst.resize(src.len(), 0.0);
            }
            dst.iter_mut().zip(src).for_each(|(a, b)| *a += b);
        }
        self
    }

    /// Censored-histogram EM step from `old`: an observed length counts in its bin; a read
    /// censored in bin `c` is spread over the bins `b >= c` in proportion to `old.h[b]`.
    fn law(&self, old: &LenLaw) -> LenLaw {
        let nb = old.h.len();
        let mut tail = vec![0.0; nb + 1];
        for b in (0..nb).rev() {
            tail[b] = tail[b + 1] + old.h[b];
        }
        let mut n = vec![0.0; nb];
        // running sum over c <= b of cens[c] / P(length in bins >= c)
        let mut acc = 0.0;
        for b in 0..nb {
            let c = self.cens.get(b).copied().unwrap_or(0.0);
            if c > 0.0 && tail[b] > 0.0 {
                acc += c / tail[b];
            }
            n[b] = self.stop.get(b).copied().unwrap_or(0.0) + acc * old.h[b];
        }
        let tot: f64 = n.iter().sum();
        if tot <= 0.0 {
            return old.clone();
        }
        // A small floor keeps every length possible.
        let eps = 1e-9 * tot;
        let z = tot + eps * nb as f64;
        LenLaw {
            h: n.into_iter().map(|x| (x + eps) / z).collect(),
        }
    }
}

/// The positional model's parameters.
#[derive(Clone, Debug)]
struct Params {
    /// Mixture weights: 3'-anchored, 5'-anchored, full length, uniformly placed.
    p3: f64,
    p5: f64,
    pf: f64,
    pu: f64,
    /// Read-length laws of the 3'-anchored, 5'-anchored and uniform components.
    len3: LenLaw,
    len5: LenLaw,
    lenu: LenLaw,
    /// End kernels at the 3' and 5' ends.
    k3: EndLaw,
    k5: EndLaw,
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
    /// Absolute and relative parts of the end kernels at this alignment's gaps.
    k3a: f64,
    k3r: f64,
    k5a: f64,
    k5r: f64,
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
    fn initial(max_len: u32) -> Self {
        Params {
            p3: 0.25,
            p5: 0.25,
            pf: 0.25,
            pu: 0.25,
            len3: LenLaw::exponential(1000.0, max_len),
            len5: LenLaw::exponential(1000.0, max_len),
            lenu: LenLaw::exponential(1000.0, max_len),
            k3: EndLaw::initial(max_len),
            k5: EndLaw::initial(max_len),
        }
    }

    #[inline]
    fn terms(&self, g: &Geo, e: &Exps) -> Terms {
        let (len, l, d5, d3) = (g.len, g.l, g.d5, g.d3);
        let (l_us, len_us) = (l as usize, len as usize);
        // end kernels, each a distribution over the gaps this transcript allows (0..len)
        let (k3a, k3r) = self.k3.parts(d3, len, &e.g3d, &e.g3c);
        let (k5a, k5r) = self.k5.parts(d5, len, &e.g5d, &e.g5c);
        let (kk3, kk5) = (k3a + k3r, k5a + k5r);
        Terms {
            a3_stop: self.p3 * kk3 * e.d3[l_us],
            a3_cens: self.p3 * kk3 * e.s3[(len - d3) as usize] * kk5,
            a5_stop: self.p5 * kk5 * e.d5[l_us],
            a5_cens: self.p5 * kk5 * e.s5[(len - d5) as usize] * kk3,
            full: self.pf * kk3 * kk5,
            u_stop: self.pu * e.du[l_us] / ((len - l) as f64 + 1.0),
            u_cens: self.pu * e.su[len_us] * kk3 * kk5,
            k3a,
            k3r,
            k5a,
            k5r,
        }
    }
}

/// `exp(-rate x)` for each rate, and `(1 - rho)^x` for each kernel scale, for every whole number
/// of bases `x` up to the longest transcript, built once per iteration so that no alignment
/// computes an exponential.
struct Exps {
    /// Length-law densities and survivals of the three components.
    d3: Vec<f64>,
    s3: Vec<f64>,
    d5: Vec<f64>,
    s5: Vec<f64>,
    du: Vec<f64>,
    su: Vec<f64>,
    /// End-kernel densities and distribution functions.
    g3d: Vec<f64>,
    g3c: Vec<f64>,
    g5d: Vec<f64>,
    g5c: Vec<f64>,
}

impl Exps {
    fn new(par: &Params, max_len: u32) -> Self {
        let (d3, s3) = par.len3.tables(max_len);
        let (d5, s5) = par.len5.tables(max_len);
        let (du, su) = par.lenu.tables(max_len);
        let (g3d, g3c) = par.k3.abs.tables(max_len);
        let (g5d, g5c) = par.k5.abs.tables(max_len);
        Exps {
            d3,
            s3,
            d5,
            s5,
            du,
            su,
            g3d,
            g3c,
            g5d,
            g5c,
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
    /// Length laws.
    l3: LStats,
    l5: LStats,
    lu: LStats,
    /// End kernels: observed gaps (`stop`) and truncation points (`cens`), binned.
    k3: EStats,
    k5: EStats,
    /// Log-likelihood of the reads at the step's input (abundances unnormalized).
    ll: f64,
}

impl Stats {
    fn merge(mut self, o: Stats) -> Stats {
        self.ll += o.ll;
        self.n += o.n;
        self.n3 += o.n3;
        self.n5 += o.n5;
        self.nf += o.nf;
        self.nu += o.nu;
        self.l3 = self.l3.merge(o.l3);
        self.l5 = self.l5.merge(o.l5);
        self.lu = self.lu.merge(o.lu);
        self.k3 = self.k3.merge(o.k3);
        self.k5 = self.k5.merge(o.k5);
        self
    }

    fn add(&mut self, g: &Geo, t: &Terms, post: f64) {
        let tot = t.total();
        if tot <= 0.0 {
            return;
        }
        let q = post / tot;
        let (a3s, a3c) = (q * t.a3_stop, q * t.a3_cens);
        let (a5s, a5c) = (q * t.a5_stop, q * t.a5_cens);
        let (us, uc) = (q * t.u_stop, q * t.u_cens);
        let fl = q * t.full;
        self.n += post;
        self.n3 += a3s + a3c;
        self.n5 += a5s + a5c;
        self.nf += fl;
        self.nu += us + uc;
        // A read reaching the far end is right-censored at the length available to it; any
        // other read is an observed stop after `l` bases.
        LStats::add(&mut self.l3.stop, LenLaw::bin_nt(g.l), a3s);
        LStats::add(&mut self.l3.cens, LenLaw::bin_nt(g.len - g.d3), a3c);
        LStats::add(&mut self.l5.stop, LenLaw::bin_nt(g.l), a5s);
        LStats::add(&mut self.l5.cens, LenLaw::bin_nt(g.len - g.d5), a5c);
        LStats::add(&mut self.lu.stop, LenLaw::bin_nt(g.l), us);
        LStats::add(&mut self.lu.cens, LenLaw::bin_nt(g.len), uc);
        // Every part that places the read at an end through that end's kernel observes its gap.
        let (w3, w5) = (a3s + a3c + a5c + fl + uc, a5s + a5c + a3c + fl + uc);
        let lu = g.len;
        self.k3.add(w3, g.d3, lu, t.k3a, t.k3r);
        self.k5.add(w5, g.d5, lu, t.k5a, t.k5r);
    }

    fn params(&self, old: &Params) -> Params {
        if self.n <= 0.0 {
            return old.clone();
        }
        let frac = |x: f64| (x / self.n).clamp(1e-6, 1.0);
        let (p3, p5, pf, pu) = (frac(self.n3), frac(self.n5), frac(self.nf), frac(self.nu));
        let z = p3 + p5 + pf + pu;
        Params {
            p3: p3 / z,
            p5: p5 / z,
            pf: pf / z,
            pu: pu / z,
            len3: self.l3.law(&old.len3),
            len5: self.l5.law(&old.len5),
            lenu: self.lu.law(&old.lenu),
            k3: EndLaw::law(&self.k3, &old.k3),
            k5: EndLaw::law(&self.k5, &old.k5),
        }
    }
}

/// Items per block of [`fold_blocks`].
const FOLD_BLOCK: usize = 1 << 12;

/// Folds `f` over `0..n` in fixed blocks of [`FOLD_BLOCK`] items: the blocks run in parallel,
/// each folding its items in order from `init()` (with a reusable `scratch`), and their results
/// are merged in block order. Floating-point sums so come out the same on every run and for
/// any number of threads, unlike rayon's `fold`/`reduce`, whose grouping follows the work
/// stealing.
fn fold_blocks<T: Send, S>(
    n: usize,
    init: impl Fn() -> T + Sync + Send,
    scratch: impl Fn() -> S + Sync + Send,
    f: impl Fn(&mut T, &mut S, usize) + Sync + Send,
    merge: impl Fn(T, T) -> T,
) -> T {
    let parts: Vec<T> = (0..n.div_ceil(FOLD_BLOCK))
        .into_par_iter()
        .map_init(scratch, |s, b| {
            let mut acc = init();
            for i in b * FOLD_BLOCK..((b + 1) * FOLD_BLOCK).min(n) {
                f(&mut acc, s, i);
            }
            acc
        })
        .collect();
    parts.into_iter().reduce(merge).unwrap_or_else(init)
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
            let mut denom = 0.0;
            for g in f.read_geos(r) {
                let t = par.terms(&g, e);
                denom += prev[g.t as usize].load(Ordering::Relaxed) * g.p as f64 * t.total();
            }
            if denom > constants::EM_DENOM_THRESH {
                for g in f.read_geos(r) {
                    let t = par.terms(&g, e);
                    let tid = g.t as usize;
                    let post =
                        m * prev[tid].load(Ordering::Relaxed) * g.p as f64 * t.total() / denom;
                    curr[tid].fetch_add(post, Ordering::AcqRel);
                    st.add(&g, &t, post);
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
        let f = Flat::new(em_info);
        let nreads = f.off.len() - 1;
        let sample = sample_reads();
        let start = if sample > 0 && nreads > sample && Layout::can_freeze_default() {
            let k = nreads.div_ceil(sample);
            let (a, par) = {
                let sub = f.every(k);
                info!(
                    "estimating the positional parameters on one read in {k} ({} reads)",
                    (sub.off.len() - 1).to_formatted_string(&Locale::en)
                );
                let share = (sub.off.len() - 1) as f64 / nreads as f64;
                run_from(em_info, &sub, None, None, Start::EstimateOnly(share))
            };
            // The sample's tables are freed; return them to the OS before the full-data ones
            // are built (the allocator would otherwise keep them resident for a while).
            crate::util::oarfish_types::release_freed_memory();
            Start::Frozen(Box::new(par), warm_start(&a, nreads as f64))
        } else {
            Start::Estimate
        };
        let (counts, par) = run_from(em_info, &f, None, presence, start);
        info!(
            "truncation model: 3'-anchored {:.3} (mean length {:.0} nt), 5'-anchored {:.3} \
             (mean length {:.0} nt), full length {:.3}, uniformly placed {:.3} (mean length {:.0} nt); \
             3' end gap {:?}, 5' end gap {:?}",
            par.p3,
            par.len3.mean(),
            par.p5,
            par.len5.mean(),
            par.pf,
            par.pu,
            par.lenu.mean(),
            par.k3,
            par.k5
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
                let per = f.n().div_ceil(n).max(1);
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
                for t in (0..f.n()).map(|k| f.target(k)) {
                    off[t as usize + 1] += 1;
                }
                for i in 0..ntx {
                    off[i + 1] += off[i];
                }
                let mut fill = off.clone();
                let mut aln = vec![0u32; f.n()];
                for (k, t) in (0..f.n()).map(|k| (k, f.target(k))) {
                    aln[fill[t as usize]] = k as u32;
                    fill[t as usize] += 1;
                }
                Layout::Csr {
                    post: vec![0.0; f.n()],
                    off,
                    aln,
                }
            }
        }
    }
}

impl Layout {
    /// Whether the default layout freezes its parameters (all but the atomic one do).
    fn can_freeze_default() -> bool {
        crate::em_read::Impl::from_env().is_some()
    }

    /// The default layout for frozen weights `fx`, built directly (as [`Layout::frozen`]
    /// would make it from the estimating layout).
    fn frozen_new(f: &Flat, fx: &Fixed, ntx: usize) -> Self {
        match crate::em_read::Impl::from_env() {
            Some(crate::em_read::Impl::Csr) => {
                let mut off = vec![0usize; ntx + 1];
                for &t in &fx.t {
                    off[t as usize + 1] += 1;
                }
                for i in 0..ntx {
                    off[i + 1] += off[i];
                }
                Layout::Csr {
                    post: Vec::new(),
                    off,
                    aln: Vec::new(),
                }
                .frozen(fx)
            }
            how => Layout::with(f, ntx, how),
        }
    }

    /// The layout for frozen weights `ws`: `Csr` becomes `CsrFrozen`, the others stay. The
    /// per-alignment posteriors and indices are freed first and the transpose is refilled
    /// from read order, so the two layouts do not coexist.
    fn frozen(self, fx: &Fixed) -> Self {
        match self {
            Layout::Csr { off, post, aln } => {
                drop((post, aln));
                let nreads = fx.nreads();
                let mut fill = off.clone();
                let mut read = vec![0u32; fx.t.len()];
                let mut w = vec![0f64; fx.t.len()];
                for (r, win) in fx.off.windows(2).enumerate() {
                    let (lo, hi) = (win[0], win[1]);
                    for (&t, &wk) in fx.t[lo..hi].iter().zip(&fx.w[lo..hi]) {
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
                scratch: &mut Vec<(Geo, Terms, f64)>,
                emit: &mut dyn FnMut(usize, usize, f64)| {
        let m = mult.map_or(1.0, |m| m[r] as f64);
        if m == 0.0 {
            return;
        }
        let (lo, hi) = (f.off[r], f.off[r + 1]);
        scratch.clear();
        let mut denom = 0.0;
        for g in f.read_geos(r) {
            let t = par.terms(&g, e);
            let w = prev[g.t as usize].load(Ordering::Relaxed) * g.p as f64 * t.total();
            denom += w;
            scratch.push((g, t, w));
        }
        if denom > constants::EM_DENOM_THRESH {
            st.ll += m * denom.ln();
            for (k, (g, t, w)) in (lo..hi).zip(scratch.iter()) {
                let post = m * w / denom;
                emit(k, g.t as usize, post);
                st.add(g, t, post);
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
                .collect::<Vec<_>>()
                .into_iter()
                .fold(Stats::default(), Stats::merge);
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
            let st = fold_blocks(
                nreads,
                Stats::default,
                Vec::new,
                |st, scratch, r| {
                    for k in f.off[r]..f.off[r + 1] {
                        // SAFETY: slot k belongs to read r alone.
                        unsafe { *(post_ptr as *mut f64).add(k) = 0.0 };
                    }
                    read(r, st, scratch, &mut |k, _, p| unsafe {
                        *(post_ptr as *mut f64).add(k) = p
                    });
                },
                Stats::merge,
            );
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
    fx: &Fixed,
    prev: &[AtomicF64],
    curr: &[AtomicF64],
    mult: Option<&[u32]>,
) {
    let read = |r: usize, emit: &mut dyn FnMut(usize, usize, f64)| {
        let m = fx.mult(mult, r);
        if m == 0.0 {
            return;
        }
        let (lo, hi) = (fx.off[r], fx.off[r + 1]);
        let mut denom = 0.0;
        for (&t, &w) in fx.t[lo..hi].iter().zip(&fx.w[lo..hi]) {
            denom += prev[t as usize].load(Ordering::Relaxed) * w;
        }
        if denom > constants::EM_DENOM_THRESH {
            let inv = m / denom;
            for (k, (&t, &w)) in (lo..hi).zip(fx.t[lo..hi].iter().zip(&fx.w[lo..hi])) {
                let t = t as usize;
                emit(k, t, prev[t].load(Ordering::Relaxed) * w * inv);
            }
        }
    };
    match lay {
        Layout::Atomic => unreachable!("parameters are not frozen in the atomic layout"),
        Layout::CsrFrozen { off, read, w, inv } => {
            inv.par_iter_mut().enumerate().for_each(|(r, x)| {
                let m = fx.mult(mult, r);
                let (lo, hi) = (fx.off[r], fx.off[r + 1]);
                let mut denom = 0.0;
                for (&t, &w) in fx.t[lo..hi].iter().zip(&fx.w[lo..hi]) {
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
        Layout::Csr { .. } => unreachable!("frozen weights use the CsrFrozen layout"),
    }
}

/// Parameters stop being re-estimated once none moves by more than this (see [`moved`]) for
/// `FREEZE_AFTER` consecutive iterations; the remaining iterations use fixed weights. Updating
/// the parameters less often is still a generalized EM; on K562 PacBio (truncation model)
/// freezing at 1e-6 moves no transcript estimate by more than 0.07% and cuts EM time 3.6x.
/// Only the atomic-free layouts freeze.
const FREEZE_TOL: f64 = 1e-6;
const FREEZE_AFTER: u32 = 3;
/// Parameters also stop being re-estimated once an iteration raises the log-likelihood by less
/// than this many nats, for `FREEZE_AFTER` consecutive iterations. The parameters can drift
/// along directions that barely change any alignment's weight (mass trading between the
/// full-length component and an anchored one whose read lengths mostly exceed the transcript)
/// long after their distance test would pass; the likelihood bounds what such drift is worth.
/// One nat over all the reads is far below what the data could distinguish. On a sample, whose
/// gains are smaller by its share of the reads, the threshold shrinks by the same share, so the
/// sample freezes where a fit to all the reads would.
const FREEZE_LL_GAIN: f64 = 1.0;

/// How far the parameters moved: absolute change for the mixture weights; for each length law,
/// its change scaled by the weight of the one component that uses it, as that bounds its effect
/// on the likelihood (the law of a component holding a negligible share of the reads drifts for
/// as long as the EM runs, and must not hold the parameters open). The end kernels enter
/// several components and count in full.
fn moved(a: &Params, b: &Params) -> f64 {
    (a.p3 - b.p3)
        .abs()
        .max((a.p5 - b.p5).abs())
        .max((a.pf - b.pf).abs())
        .max((a.pu - b.pu).abs())
        .max(a.p3.max(b.p3) * a.len3.moved(&b.len3))
        .max(a.p5.max(b.p5) * a.len5.moved(&b.len5))
        .max(a.pu.max(b.pu) * a.lenu.moved(&b.lenu))
        .max(a.k3.moved(&b.k3))
        .max(a.k5.moved(&b.k5))
}

/// State of the SQUAREM acceleration (Varadhan & Roland 2008, scheme S3) of the frozen-weight
/// EM: buffers for the iterates of a cycle, and the bound on the extrapolation step length,
/// which grows while long steps pay off and shrinks when one is rejected.
struct Squarem {
    t0: Vec<f64>,
    t1: Vec<f64>,
    t2: Vec<f64>,
    step_max: f64,
    /// Cycles that extrapolated and kept the result, rejected it, or did not extrapolate.
    acc: u32,
    rej: u32,
    plain: u32,
}

impl Squarem {
    const STEP_GROW: f64 = 4.0;

    fn new(n: usize) -> Self {
        Squarem {
            t0: vec![0.0; n],
            t1: vec![0.0; n],
            t2: vec![0.0; n],
            step_max: 1.0,
            acc: 0,
            rej: 0,
            plain: 0,
        }
    }
}

fn load_into(dst: &mut [f64], src: &[AtomicF64]) {
    dst.par_iter_mut()
        .zip(src.par_iter())
        .for_each(|(d, s)| *d = s.load(Ordering::Relaxed));
}

fn store_from(dst: &[AtomicF64], src: &[f64]) {
    dst.par_iter()
        .zip(src.par_iter())
        .for_each(|(d, &s)| d.store(s, Ordering::Relaxed));
}

/// Log-likelihood of abundances `theta` (normalized to proportions) under the fixed weights
/// `ws`, up to a constant: the objective the frozen-weight EM ascends.
fn loglik_static(fx: &Fixed, theta: &[f64], mult: Option<&[u32]>) -> f64 {
    let total: f64 = theta.iter().sum();
    let (ll, n) = fold_blocks(
        fx.nreads(),
        || (0.0, 0.0),
        || (),
        |acc, _, r| {
            let m = fx.mult(mult, r);
            let (lo, hi) = (fx.off[r], fx.off[r + 1]);
            let d: f64 = fx.t[lo..hi]
                .iter()
                .zip(&fx.w[lo..hi])
                .map(|(&t, &w)| theta[t as usize] * w)
                .sum();
            if m > 0.0 && d > constants::EM_DENOM_THRESH {
                acc.0 += m * d.ln();
                acc.1 += m;
            }
        },
        |a, b| (a.0 + b.0, a.1 + b.1),
    );
    ll - n * total.ln()
}

/// One SQUAREM cycle from the abundances in `prev`: two EM steps (θ1, θ2), a step along the
/// extrapolated path from θ0, and one EM step from there to stabilize it. The result is kept
/// only if its likelihood is at least that of θ2, the result of plain EM; otherwise θ2 is kept
/// and the step bound shrinks. On return `curr` holds the new abundances and `prev` the input
/// of the cycle's last EM step, so their difference is a one-step change as in plain EM.
fn squarem_cycle(
    lay: &mut Layout,
    fx: &Fixed,
    prev: &[AtomicF64],
    curr: &[AtomicF64],
    mult: Option<&[u32]>,
    sq: &mut Squarem,
) {
    load_into(&mut sq.t0, prev);
    step_static(lay, fx, prev, curr, mult);
    load_into(&mut sq.t1, curr);
    store_from(prev, &sq.t1);
    step_static(lay, fx, prev, curr, mult);
    load_into(&mut sq.t2, curr);
    let (rr, vv) = fold_blocks(
        sq.t0.len(),
        || (0.0, 0.0),
        || (),
        |acc, _, i| {
            let (a, b, c) = (sq.t0[i], sq.t1[i], sq.t2[i]);
            let (r, v) = (b - a, c - 2.0 * b + a);
            acc.0 += r * r;
            acc.1 += v * v;
        },
        |x, y| (x.0 + y.0, x.1 + y.1),
    );
    // Step length -alpha in [1, step_max]; 1 reproduces θ2 exactly.
    let mut s = if vv > 0.0 { (rr / vv).sqrt() } else { 1.0 };
    if !s.is_finite() || s < 1.0 {
        s = 1.0;
    }
    let hit = s > 1.0 && s >= sq.step_max;
    s = s.min(sq.step_max);
    if s == 1.0 {
        // Plain EM: θ2 stands; report its last step. A step held back by the bound lets the
        // bound grow, as after an accepted long step.
        if hit {
            sq.step_max *= Squarem::STEP_GROW;
        }
        sq.plain += 1;
        store_from(prev, &sq.t1);
        return;
    }
    let ext: Vec<f64> = sq
        .t0
        .par_iter()
        .zip(sq.t1.par_iter().zip(sq.t2.par_iter()))
        .map(|(&a, (&b, &c))| {
            let (r, v) = (b - a, c - 2.0 * b + a);
            (a + 2.0 * s * r + s * s * v).max(0.0)
        })
        .collect();
    store_from(prev, &ext);
    step_static(lay, fx, prev, curr, mult);
    let t3: Vec<f64> = curr.par_iter().map(|x| x.load(Ordering::Relaxed)).collect();
    let (l3, l2) = (
        loglik_static(fx, &t3, mult),
        loglik_static(fx, &sq.t2, mult),
    );
    if l3.is_finite() && l3 >= l2 {
        sq.acc += 1;
        if hit {
            sq.step_max *= Squarem::STEP_GROW;
        }
    } else {
        sq.rej += 1;
        sq.step_max = (sq.step_max / Squarem::STEP_GROW).max(1.0);
        store_from(prev, &sq.t1);
        store_from(curr, &sq.t2);
    }
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
        let reads = crate::presence::ReadAlns {
            off: &f.off,
            target: &|k| f.target(k),
        };
        p.refresh(
            &counts,
            reads,
            |k| match (ws, &e) {
                (Some(ws), _) => ws[k],
                (None, Some(e)) => {
                    let g = f.geo(k);
                    g.p as f64 * par.terms(&g, e).total()
                }
                (None, None) => unreachable!(),
            },
            |r| mult.map_or(1.0, |m| m[r] as f64),
        );
    }
    eff.par_iter()
        .zip(counts.par_iter().zip(&p.q))
        .for_each(|(e, (c, q))| e.store(q * c, Ordering::Relaxed));
}

/// How [`run_from`] begins.
enum Start {
    /// Estimate the parameters, freeze them once they settle, and continue with them fixed.
    Estimate,
    /// As `Estimate` on a sample holding this share of the reads, but return as soon as the
    /// parameters freeze, with the abundances then.
    EstimateOnly(f64),
    /// With these parameters fixed from the start, and these initial abundances.
    Frozen(Box<Params>, Vec<f64>),
}

/// Initial abundances for all reads from the abundances `a` estimated on a sample: `a` scaled to
/// `total` reads, with a share spread evenly over the transcripts, so that those the sample
/// missed can still gain reads (the EM cannot revive a zero).
fn warm_start(a: &[f64], total: f64) -> Vec<f64> {
    let lambda = WARM_EVEN;
    let sum: f64 = a.iter().sum();
    let even = total / a.len().max(1) as f64;
    if sum <= 0.0 {
        return vec![even; a.len()];
    }
    a.iter()
        .map(|&x| (1.0 - lambda) * x * total / sum + lambda * even)
        .collect()
}

/// With more reads than this, the positional parameters are estimated on a sample of about
/// this many (every k-th read), and the EM over all reads starts with them fixed
/// ([`Start::Frozen`]). They are a few global laws, which a million reads determine about as
/// well as all of them: on simulated cDNA, direct RNA and HiFi samples the mixture weights and
/// the laws of the components holding reads match the full-data fit, and the estimates agree
/// once both have converged, while each estimating iteration costs several frozen ones.
const SAMPLE_READS: usize = 1_000_000;
/// Share of the reads spread evenly over the transcripts in the initial abundances after a
/// sample ([`warm_start`]).
const WARM_EVEN: f64 = 0.1;
/// Plain EM steps before SQUAREM, when the EM starts with the parameters fixed. Without them the
/// estimates of transcripts with few or no reads have not settled within the default iteration
/// limit (as they have after the estimating iterations of a full-data fit).
const PLAIN_FIRST: u32 = 100;

/// [`SAMPLE_READS`], or `OARFISH_TRUNC_SAMPLE` if set (0: estimate on all reads).
fn sample_reads() -> usize {
    std::env::var("OARFISH_TRUNC_SAMPLE")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(SAMPLE_READS)
}

fn run(
    em_info: &EMInfo,
    f: &Flat,
    mult: Option<&[u32]>,
    presence: Option<PresenceOpts>,
) -> (Vec<f64>, Params) {
    run_from(em_info, f, mult, presence, Start::Estimate)
}

fn run_from(
    em_info: &EMInfo,
    f: &Flat,
    mult: Option<&[u32]>,
    presence: Option<PresenceOpts>,
    start: Start,
) -> (Vec<f64>, Params) {
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let total_weight = em_info.eq_map.num_aligned_reads() as f64;
    let init = match (&start, &em_info.init_abundances) {
        (Start::Frozen(_, a), _) => a.clone(),
        (_, Some(v)) => v.clone(),
        (_, None) => vec![total_weight / (tinfo.len() as f64); tinfo.len()],
    };
    let mut prev: Vec<AtomicF64> = init.into_iter().map(AtomicF64::new).collect();
    let mut curr: Vec<AtomicF64> = (0..tinfo.len()).map(|_| AtomicF64::new(0.0)).collect();
    let t_setup = std::time::Instant::now();
    let max_len = f.max_len();
    let started_frozen = matches!(start, Start::Frozen(..));
    let freeze_ll_gain = match start {
        Start::EstimateOnly(share) => FREEZE_LL_GAIN * share,
        _ => FREEZE_LL_GAIN,
    };
    let (mut par, mut lay, mut frozen, estimate_only) = match start {
        Start::Frozen(par, _) => {
            let fx = Fixed::new(f, &par);
            let lay = Layout::frozen_new(f, &fx, tinfo.len());
            (*par, lay, Some(fx), false)
        }
        s => (
            Params::initial(max_len),
            Layout::new(f, tinfo.len()),
            None,
            matches!(s, Start::EstimateOnly(_)),
        ),
    };
    info!(
        "truncation EM layout set up in {:.2}s",
        t_setup.elapsed().as_secs_f64()
    );

    let mut niter = 0_u32;
    let can_freeze = !matches!(lay, Layout::Atomic);
    let mut still: u32 = 0;
    // Log-likelihood gain of the last parameter-updating step.
    let (mut last_ll, mut ll_gain) = (f64::NEG_INFINITY, f64::INFINITY);
    // With the presence model the EM sees `q * prev` (see `crate::presence`); `eff` holds it.
    let mut pres = presence.map(|o| Presence::new(tinfo.len(), o));
    let eff: Vec<AtomicF64> = if pres.is_some() {
        (0..tinfo.len()).map(|_| AtomicF64::new(0.0)).collect()
    } else {
        Vec::new()
    };
    let mut sq = Squarem::new(tinfo.len());
    let plain_first = if started_frozen { PLAIN_FIRST } else { 0 };
    while niter < em_info.max_iter {
        if let Some(p) = &mut pres {
            refresh_presence(
                p,
                niter,
                f,
                &par,
                frozen.as_ref().map(|x| x.w.as_slice()),
                &prev,
                &eff,
                mult,
                max_len,
            );
        }
        let src = if pres.is_some() { &eff } else { &prev };
        match &frozen {
            // With the weights fixed, the abundance update is a fixed-point map; accelerate
            // it (the presence model rescales the abundances between steps, so it does not).
            Some(ws) if pres.is_none() && niter >= plain_first => {
                squarem_cycle(&mut lay, ws, &prev, &curr, mult, &mut sq);
                niter += 2;
            }
            Some(ws) => step_static(&mut lay, ws, src, &curr, mult),
            None => {
                let e = Exps::new(&par, max_len);
                let st = step_with(&mut lay, f, &par, src, &curr, mult, &e);
                ll_gain = st.ll - last_ll;
                last_ll = st.ll;
                let next = st.params(&par);
                still = if moved(&next, &par) < FREEZE_TOL
                    || (0.0..freeze_ll_gain).contains(&ll_gain)
                {
                    still + 1
                } else {
                    0
                };
                par = next;
                if can_freeze && still >= FREEZE_AFTER {
                    let fx = Fixed::new(f, &par);
                    info!(
                        "truncation parameters frozen after {} iterations: {:?}",
                        niter + 1,
                        par
                    );
                    if estimate_only {
                        let a = curr.iter().map(|x| x.load(Ordering::Relaxed)).collect();
                        return (a, par);
                    }
                    lay = std::mem::replace(&mut lay, Layout::Atomic).frozen(&fx);
                    frozen = Some(fx);
                }
            }
        }
        // Largest relative change in either direction. (The per-read EM counts increases only;
        // after a SQUAREM extrapolation most of the motion is in decreases, which that test
        // ignores, so it would stop at an arbitrary point.)
        let mut rel_diff = 0.0_f64;
        for (p, c) in prev.iter().zip(&curr) {
            let pc = p.load(Ordering::Relaxed);
            if pc > constants::MIN_READ_THRESH {
                rel_diff = rel_diff.max((c.load(Ordering::Relaxed) - pc).abs() / pc);
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        curr.par_iter()
            .for_each(|x| x.store(0.0, Ordering::Relaxed));
        if rel_diff < em_info.convergence_thresh && niter > 1 {
            break;
        }
        niter += 1;
        if niter.is_multiple_of(100) && frozen.is_some() {
            // The parameters were logged when they froze.
            info!(
                "iteration {}; rel diff {}",
                niter.to_formatted_string(&Locale::en),
                rel_diff
            );
        } else if niter.is_multiple_of(100) {
            info!(
                "iteration {}; rel diff {}; log-likelihood gain {:.3}; {:?}",
                niter.to_formatted_string(&Locale::en),
                rel_diff,
                ll_gain,
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
    if sq.acc + sq.rej + sq.plain > 0 {
        info!(
            "truncation EM: {} iterations; SQUAREM cycles {} extrapolated, {} rejected, {} plain (step bound {})",
            niter, sq.acc, sq.rej, sq.plain, sq.step_max
        );
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
                "presence model: {} transcripts called absent (q < 0.5){}",
                p.n_absent().to_formatted_string(&Locale::en),
                p.base_rate()
                    .map_or(String::new(), |pi| format!("; base rate {pi:.3}"))
            );
        }
        let counts: Vec<f64> = prev.iter().map(|x| x.load(Ordering::Relaxed)).collect();
        for (e, (c, q)) in eff.iter().zip(counts.iter().zip(&p.q)) {
            e.store(q * c, Ordering::Relaxed);
        }
    }
    let src = if pres.is_some() { &eff } else { &prev };
    match &frozen {
        Some(ws) => step_static(&mut lay, ws, src, &curr, mult),
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

    fn flat(alns: &[Vec<(u32, f64, f64, f64)>], lens: &[f64]) -> Flat<'static> {
        let mut f = Flat::empty(lens);
        for r in alns {
            for &(t, p, s, e) in r {
                f.push(t, p as f32, s, e);
            }
            f.off.push(f.n());
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
        let mut par = Params::initial(f.max_len());
        for _ in 0..iters {
            let e = Exps::new(&par, f.max_len());
            par = step_with(&mut lay, f, &par, &prev, &curr, mult, &e).params(&par);
            std::mem::swap(&mut prev, &mut curr);
            curr.iter().for_each(|x| x.store(0.0, Ordering::Relaxed));
        }
        (
            prev.iter().map(|x| x.load(Ordering::Relaxed)).collect(),
            par,
        )
    }

    /// A sample keeps every k-th read whole, and the warm start spreads the full read count
    /// over all transcripts, none of them at zero.
    #[test]
    fn sample_keeps_whole_reads_and_warm_start_covers_every_transcript() {
        let lens = [2000.0, 1000.0, 1500.0];
        let reads: Vec<Vec<_>> = (0..10)
            .map(|i| {
                (0..i % 3 + 1)
                    .map(|j| (j as u32, 1.0, 10.0 * i as f64, 900.0))
                    .collect()
            })
            .collect();
        let f = flat(&reads, &lens);
        let s = f.every(3);
        assert_eq!(s.off.len() - 1, 4);
        for (i, r) in [0, 3, 6, 9].into_iter().enumerate() {
            assert_eq!(
                &s.alns[s.off[i]..s.off[i + 1]],
                &f.alns[f.off[r]..f.off[r + 1]]
            );
        }
        let a = warm_start(&[4.0, 0.0, 1.0], 100.0);
        assert!((a.iter().sum::<f64>() - 100.0).abs() < 1e-9);
        assert!(a.iter().all(|&x| x > 0.0) && a[0] > a[2]);
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
            pa.len3.moved(&pb.len3) < 1e-9 && (pa.p3 - pb.p3).abs() < 1e-9,
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
                    pa.len3.moved(&pb.len3) < 1e-12 && (pa.p3 - pb.p3).abs() < 1e-12,
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
        let par = Params::initial(f.max_len());
        let ws = Fixed::new(&f, &par);
        let prev: Vec<AtomicF64> = [3.0, 1.0, 2.0].into_iter().map(AtomicF64::new).collect();
        let mult: Vec<u32> = (0..reads.len()).map(|i| (i % 3) as u32).collect();
        for m in [None, Some(&mult[..])] {
            let round = |how| {
                let mut lay = Layout::with(&f, 3, Some(how)).frozen(&ws);
                let curr: Vec<AtomicF64> = (0..3).map(|_| AtomicF64::new(0.0)).collect();
                step_static(&mut lay, &ws, &prev, &curr, m);
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
        assert!(pa.len3.moved(&pb.len3) < 1e-12 && pa.lenu.moved(&pb.lenu) < 1e-12);
    }

    /// SQUAREM reaches the fixed point of the frozen-weight EM, in fewer EM evaluations, and
    /// never lowers the likelihood from one cycle to the next.
    #[test]
    fn squarem_reaches_the_em_fixed_point() {
        use crate::em_read::Impl;
        let lens = [2000.0, 1800.0, 1500.0, 900.0];
        let mut reads = Vec::new();
        for i in 0..60 {
            let s = (i % 60) as f64 * 10.0;
            // Nested, nearly indistinguishable targets: the slow case for plain EM.
            reads.push(vec![
                (0, 1.0, s, s + 700.0),
                (1, 1.0, s, s + 700.0),
                (2, 0.9, s, s + 700.0),
            ]);
            if i % 3 == 0 {
                reads.push(vec![(0, 1.0, 1300.0 + s, 2000.0)]);
            }
            if i % 5 == 0 {
                reads.push(vec![(3, 1.0, 0.0, 900.0), (2, 0.5, 600.0, 1500.0)]);
            }
        }
        let f = flat(&reads, &lens);
        let par = Params::initial(f.max_len());
        let ws = Fixed::new(&f, &par);
        let init = || -> Vec<AtomicF64> { (0..4).map(|_| AtomicF64::new(100.0)).collect() };
        let fresh = || -> Vec<AtomicF64> { (0..4).map(|_| AtomicF64::new(0.0)).collect() };
        let get = |v: &[AtomicF64]| {
            v.iter()
                .map(|x| x.load(Ordering::Relaxed))
                .collect::<Vec<_>>()
        };

        let mut lay = Layout::with(&f, 4, Some(Impl::Csr)).frozen(&ws);
        let (prev, curr) = (init(), fresh());
        for _ in 0..5_000 {
            step_static(&mut lay, &ws, &prev, &curr, None);
            store_from(&prev, &get(&curr));
        }
        let plain = get(&prev);

        let mut lay = Layout::with(&f, 4, Some(Impl::Csr)).frozen(&ws);
        let (prev, curr) = (init(), fresh());
        let mut sq = Squarem::new(4);
        let mut last = f64::NEG_INFINITY;
        for _ in 0..100 {
            squarem_cycle(&mut lay, &ws, &prev, &curr, None, &mut sq);
            let now = get(&curr);
            let ll = loglik_static(&ws, &now, None);
            assert!(ll >= last - 1e-9, "likelihood fell: {last} -> {ll}");
            last = ll;
            store_from(&prev, &now);
        }
        let fast = get(&prev);
        assert!(sq.acc > 0, "never extrapolated");
        for (x, y) in plain.iter().zip(&fast) {
            assert!((x - y).abs() < 1e-3 * x.max(1.0), "{plain:?} {fast:?}");
        }
    }
}
