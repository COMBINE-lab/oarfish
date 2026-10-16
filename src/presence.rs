//! Spike-and-slab presence model (`--presence-model spike-slab`), shared by the EM paths.
//!
//! Each transcript carries a Bernoulli presence posterior `q_t`. The EM sees the *effective*
//! abundance `q_t * counts[t]`, which is the zero-inflated EM update: a transcript with low
//! presence probability competes for shared reads at a discount and its ambiguous mass drains to
//! its supported competitors, while reads unique to a transcript keep it (their evidence drives
//! `q_t -> 1`). `q` is refreshed every `period` iterations after `warmup` from a leave-one-out
//! evidence pass:
//!
//! ```text
//! delta_t = sum over reads touching t of  m_r * ln( (rest_r + slab_t * w) / rest_r )
//! q_t     = sigmoid( delta_t + logit(rho) )
//! ```
//!
//! where `w` is the read's (summed) weight on `t`, `rest_r` the read's normalizer without `t`,
//! and `slab_t = counts[t] / max(q_t, 1e-3)` the abundance conditional on presence. `delta_t` is
//! a leave-one-out log-likelihood gain and so is never negative; the prior `rho` is therefore a
//! fixed penalty (an empirical-Bayes `rho = mean(q)` only ratchets upward). For a transcript
//! holding small posterior shares, `delta_t` is about its assigned read count, so the default
//! `rho = 0.05` (`logit ~ -2.94`) asks for roughly three effective exclusive reads for `q > 0.5`.
//!
//! Ported from `oarfish2-dev` (`em.rs`, `PresenceState`) without its optional evidence terms
//! (endpoint, 5'-flush, external prior).
//!
//! # Bayes-factor presence (`--presence-model bayes`)
//!
//! The leave-one-out gain above is a *maximized* likelihood ratio: never negative, blind to how
//! sharply the data pin a transcript's abundance, so its prior must be a fixed penalty. That
//! cannot tell a transcript the data support from one along a direction the data cannot
//! resolve (a nested isoform whose reads fit its longer sibling equally well), where noise is
//! amplified into hundreds of reads. The Bayes-factor model replaces it with the marginal
//! likelihood of each transcript's presence. Moving a share `a` of the reads to transcript `t`
//! from the rest of the current fit changes the log-likelihood by
//!
//! ```text
//! l_t(a) = sum over reads r touching t of m_r ln(1 + a (R_r - 1))  +  (N - n_t) ln(1 - a)
//! ```
//!
//! where `R_r` is read `r`'s likelihood under `t` over its likelihood under the rest of the fit
//! (shares from the current effective abundances, `t` removed), `m_r` its multiplicity, `N` the
//! total and `n_t` the reads touching `t`. With a slab `g(a)` (the distribution of an expressed
//! transcript's share), the Bayes factor for presence is `BF_t = E_g[exp(l_t(a))]`, evaluated
//! on a fixed log-spaced grid of shares. A flat direction gets `BF ~ posterior width / prior
//! width < 1` (an Occam penalty scaled to how sharply the data constrain it), a well-supported
//! transcript `BF >> 1`, and one with an exclusive read is present outright. Then
//! `q_t = sigmoid(ln BF_t + logit pi)`. Both the base rate `pi` and the slab `g` (a discrete
//! distribution on the grid) are learned by empirical Bayes from all transcripts: expressed
//! abundances are heavy-tailed, so a slab fixed at the mean share would charge every genuinely
//! low-abundance isoform an Occam penalty. Learning `pi` is possible because `BF` can fall below
//! one.

use rayon::prelude::*;

/// Which presence evidence to use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PresenceKind {
    /// Leave-one-out maximized gain, fixed prior `rho`.
    SpikeSlab,
    /// Marginal-likelihood Bayes factor, learned base rate.
    Bayes,
}

/// Options of the presence model.
#[derive(Clone, Copy, Debug)]
pub struct PresenceOpts {
    pub kind: PresenceKind,
    pub rho: f64,
    pub warmup: u32,
    pub period: u32,
}

/// Reads as ranges of alignments: read `r` owns alignments `off[r]..off[r + 1]`, whose targets
/// are `target[k]`.
#[derive(Clone, Copy)]
pub struct ReadAlns<'a> {
    pub off: &'a [usize],
    pub target: &'a [u32],
}

/// Grid of shares `a` on which each transcript's evidence `l_t(a)` is evaluated: log-spaced
/// from `1e-10` to `0.5`.
const GRID: usize = 64;

fn grid() -> [f64; GRID] {
    let (lo, hi) = (-10.0_f64, 0.5_f64.log10());
    std::array::from_fn(|j| 10f64.powf(lo + (hi - lo) * j as f64 / (GRID - 1) as f64))
}

/// Transcript-major index of the alignments, for the per-transcript evidence pass: transcript
/// `t`'s alignments are `aln[t_off[t]..t_off[t + 1]]`, in read order; `read_of[k]` is the read
/// owning alignment `k`.
struct Index {
    t_off: Vec<usize>,
    aln: Vec<u32>,
    read_of: Vec<u32>,
}

impl Index {
    fn new(reads: ReadAlns, ntx: usize) -> Self {
        let mut t_off = vec![0usize; ntx + 1];
        for &t in reads.target {
            t_off[t as usize + 1] += 1;
        }
        for i in 0..ntx {
            t_off[i + 1] += t_off[i];
        }
        let mut fill = t_off.clone();
        let mut aln = vec![0u32; reads.target.len()];
        let mut read_of = vec![0u32; reads.target.len()];
        for (r, w) in reads.off.windows(2).enumerate() {
            for (k, &t) in reads.target.iter().enumerate().take(w[1]).skip(w[0]) {
                let t = t as usize;
                aln[fill[t]] = k as u32;
                fill[t] += 1;
                read_of[k] = r as u32;
            }
        }
        Index {
            t_off,
            aln,
            read_of,
        }
    }
}

/// Per-read clamp on the log-evidence term (a single unique read saturates here), overall clamp
/// on `delta` before the sigmoid, and the floor that keeps the ratio finite.
const READ_LOG_CAP: f64 = 50.0;
const DELTA_CAP: f64 = 40.0;
const EPS: f64 = 1e-100;

pub struct Presence {
    pub q: Vec<f64>,
    kind: PresenceKind,
    logit_rho: f64,
    warmup: u32,
    period: u32,
    /// Learned base rate and slab (Bayes-factor model; the slab is a distribution over
    /// [`grid`]).
    pi: f64,
    slab: Vec<f64>,
    index: Option<Index>,
}

impl Presence {
    pub fn new(ntx: usize, opts: PresenceOpts) -> Self {
        let p = opts.rho.clamp(1e-4, 1.0 - 1e-4);
        Presence {
            q: vec![1.0; ntx],
            kind: opts.kind,
            logit_rho: (p / (1.0 - p)).ln(),
            warmup: opts.warmup,
            period: opts.period.max(1),
            pi: 0.5,
            slab: vec![1.0 / GRID as f64; GRID],
            index: None,
        }
    }

    /// Whether `q` is refreshed before iteration `iter` (0-based).
    pub fn due(&self, iter: u32) -> bool {
        iter > self.warmup && (iter - self.warmup).is_multiple_of(self.period)
    }

    /// `out[t] = q[t] * counts[t]`.
    pub fn effective(&self, counts: &[f64], out: &mut [f64]) {
        out.par_iter_mut()
            .zip(counts.par_iter().zip(&self.q))
            .for_each(|(o, (c, q))| *o = q * c);
    }

    /// Refresh `q` from the reads at abundances `counts`: `weight(k)` is alignment `k`'s
    /// likelihood weight, `mult(r)` read `r`'s multiplicity.
    pub fn refresh<W, M>(&mut self, counts: &[f64], reads: ReadAlns, weight: W, mult: M)
    where
        W: Fn(usize) -> f64 + Sync,
        M: Fn(usize) -> f64 + Sync,
    {
        match self.kind {
            PresenceKind::SpikeSlab => self.update(counts, reads.off.len() - 1, |r, buf| {
                for k in reads.off[r]..reads.off[r + 1] {
                    buf.push((reads.target[k], weight(k)));
                }
                mult(r)
            }),
            PresenceKind::Bayes => {
                let w: Vec<f64> = (0..reads.target.len())
                    .into_par_iter()
                    .map(&weight)
                    .collect();
                self.update_bayes(counts, reads, &w, mult)
            }
        }
    }

    /// The learned base rate, for the Bayes-factor model.
    pub fn base_rate(&self) -> Option<f64> {
        (self.kind == PresenceKind::Bayes).then_some(self.pi)
    }

    /// Bayes-factor refresh (see the module documentation).
    fn update_bayes<M>(&mut self, counts: &[f64], reads: ReadAlns, w: &[f64], mult: M)
    where
        M: Fn(usize) -> f64 + Sync,
    {
        let ntx = self.q.len();
        let idx = self.index.get_or_insert_with(|| Index::new(reads, ntx));
        let eff: Vec<f64> = counts.iter().zip(&self.q).map(|(c, q)| c * q).collect();
        let e_tot: f64 = eff.iter().sum();
        let nreads = reads.off.len() - 1;
        let (d, m): (Vec<f64>, Vec<f64>) = (0..nreads)
            .into_par_iter()
            .map(|r| {
                let d: f64 = (reads.off[r]..reads.off[r + 1])
                    .map(|k| eff[reads.target[k] as usize] * w[k])
                    .sum();
                (d, mult(r))
            })
            .unzip();
        let n_tot: f64 = m.iter().sum();
        let g = grid();
        // Per transcript: l_t on the grid, and whether some read has no other explanation.
        let ev: Vec<([f64; GRID], bool)> = (0..ntx)
            .into_par_iter()
            .map(|t| {
                let mut ell = [0.0_f64; GRID];
                let (mut n_t, mut exclusive) = (0.0_f64, false);
                let rest_tot = (e_tot - eff[t]).max(0.0);
                let alns = &idx.aln[idx.t_off[t]..idx.t_off[t + 1]];
                let mut i = 0;
                while i < alns.len() {
                    // t's alignments of one read are adjacent; sum their weights
                    let r = idx.read_of[alns[i] as usize] as usize;
                    let mut w_tr = 0.0;
                    while i < alns.len() && idx.read_of[alns[i] as usize] as usize == r {
                        w_tr += w[alns[i] as usize];
                        i += 1;
                    }
                    let m_r = m[r];
                    if m_r == 0.0 || w_tr <= 0.0 {
                        continue;
                    }
                    let rest = d[r] - eff[t] * w_tr;
                    if rest <= 1e-300 || rest_tot <= 0.0 {
                        exclusive = true;
                        continue;
                    }
                    let rm1 = w_tr * rest_tot / rest - 1.0;
                    for (l, a) in ell.iter_mut().zip(&g) {
                        *l += m_r * (a * rm1).ln_1p();
                    }
                    n_t += m_r;
                }
                for (l, a) in ell.iter_mut().zip(&g) {
                    *l += (n_tot - n_t) * (-a).ln_1p();
                }
                (ell, exclusive)
            })
            .collect();
        // Empirical Bayes for the base rate and the slab, with the evidence fixed. Each
        // transcript's evidence is kept as exp(l - max l) with its max, so a slab update needs
        // only products. A transcript with an exclusive read is present, and its share is taken
        // from its current abundance.
        let tot = counts.iter().sum::<f64>().max(1e-300);
        let prep: Vec<(Vec<f64>, f64, Option<usize>)> = ev
            .par_iter()
            .zip(counts.par_iter())
            .map(|((ell, excl), &c)| {
                if *excl {
                    let a = (c / tot).max(g[0]);
                    let j = g.iter().position(|&x| x >= a).unwrap_or(GRID - 1);
                    (Vec::new(), 0.0, Some(j))
                } else {
                    let mx = ell.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                    (ell.iter().map(|l| (l - mx).exp()).collect(), mx, None)
                }
            })
            .collect();
        let logit = |p: f64| (p / (1.0 - p)).ln();
        let (mut pi, mut slab) = (self.pi, self.slab.clone());
        let mut q = vec![0.0_f64; ntx];
        for _ in 0..200 {
            let lp = logit(pi);
            let stats: Vec<(f64, Vec<f64>)> = prep
                .par_iter()
                .map(|(e, mx, excl)| {
                    if let Some(j) = excl {
                        let mut post = vec![0.0; GRID];
                        post[*j] = 1.0;
                        return (1.0, post);
                    }
                    let z: f64 = e.iter().zip(&slab).map(|(x, p)| x * p).sum();
                    let lbf = (z.max(1e-300).ln() + mx + lp).clamp(-700.0, 700.0);
                    let qt = 1.0 / (1.0 + (-lbf).exp());
                    let post = if z > 0.0 {
                        e.iter().zip(&slab).map(|(x, p)| qt * x * p / z).collect()
                    } else {
                        vec![0.0; GRID]
                    };
                    (qt, post)
                })
                .collect();
            let mut acc = vec![0.0_f64; GRID];
            for (t, (qt, post)) in stats.iter().enumerate() {
                q[t] = *qt;
                acc.iter_mut().zip(post).for_each(|(a, p)| *a += p);
            }
            let next = (q.iter().sum::<f64>() / ntx as f64).clamp(1e-4, 1.0 - 1e-4);
            // A small floor keeps every share possible.
            let z: f64 = acc.iter().sum::<f64>() + 1e-6 * GRID as f64;
            let next_slab: Vec<f64> = acc.iter().map(|a| (a + 1e-6) / z).collect();
            let moved = (next - pi).abs().max(
                next_slab
                    .iter()
                    .zip(&slab)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
            );
            pi = next;
            slab = next_slab;
            if moved < 1e-7 {
                break;
            }
        }
        self.pi = pi;
        self.slab = slab;
        self.q = q;
    }

    /// Refresh `q` from the leave-one-out evidence of `nreads` reads at abundances `counts`.
    /// `read(r, buf)` fills `buf` with read `r`'s (target, weight) pairs and returns its
    /// multiplicity (0 skips the read).
    fn update<F>(&mut self, counts: &[f64], nreads: usize, read: F)
    where
        F: Fn(usize, &mut Vec<(u32, f64)>) -> f64 + Sync,
    {
        let ntx = self.q.len();
        let q = &self.q;
        let (delta, touched) = (0..nreads)
            .into_par_iter()
            .fold(
                || (vec![0.0_f64; ntx], vec![false; ntx], Vec::new()),
                |(mut delta, mut touched, mut buf), r| {
                    buf.clear();
                    let m = read(r, &mut buf);
                    if m == 0.0 || buf.is_empty() {
                        return (delta, touched, buf);
                    }
                    let denom: f64 = buf
                        .iter()
                        .map(|&(t, w)| q[t as usize] * counts[t as usize] * w)
                        .sum();
                    for (i, &(t, w)) in buf.iter().enumerate() {
                        if buf[..i].iter().any(|&(u, _)| u == t) {
                            continue; // this target was handled for this read
                        }
                        let w_sum: f64 = w + buf[i + 1..]
                            .iter()
                            .filter(|&&(u, _)| u == t)
                            .map(|&(_, wu)| wu)
                            .sum::<f64>();
                        let t = t as usize;
                        touched[t] = true;
                        let rest = (denom - q[t] * counts[t] * w_sum).max(0.0);
                        let slab = counts[t] / q[t].max(1e-3);
                        let term = ((rest + slab * w_sum + EPS) / (rest + EPS))
                            .ln()
                            .min(READ_LOG_CAP);
                        delta[t] += m * term;
                    }
                    (delta, touched, buf)
                },
            )
            .map(|(d, t, _)| (d, t))
            .reduce(
                || (vec![0.0_f64; ntx], vec![false; ntx]),
                |(mut d1, mut t1), (d2, t2)| {
                    d1.iter_mut().zip(&d2).for_each(|(a, b)| *a += b);
                    t1.iter_mut().zip(&t2).for_each(|(a, b)| *a |= b);
                    (d1, t1)
                },
            );
        let lr = self.logit_rho;
        self.q
            .par_iter_mut()
            .zip(delta.par_iter().zip(&touched))
            .for_each(|(q, (&d, &touched))| {
                if touched {
                    let d = d.max(0.0).clamp(-DELTA_CAP, DELTA_CAP);
                    *q = 1.0 / (1.0 + (-(d + lr)).exp());
                }
            });
    }

    /// Transcripts with presence posterior below one half.
    pub fn n_absent(&self) -> usize {
        self.q.iter().filter(|&&q| q < 0.5).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A transcript whose reads are all shared with a dominant one, and that holds only a small
    /// share of them, is called absent; one with exclusive reads stays present.
    #[test]
    fn shared_only_transcript_is_suppressed() {
        // t0: 100 exclusive reads; t1: shares 50 reads with t0; t2: 10 exclusive reads.
        let mut reads: Vec<Vec<(u32, f64)>> = Vec::new();
        reads.extend((0..100).map(|_| vec![(0, 1.0)]));
        reads.extend((0..50).map(|_| vec![(0, 1.0), (1, 1.0)]));
        reads.extend((0..10).map(|_| vec![(2, 1.0)]));
        let counts = [148.0, 2.0, 10.0];
        let mut p = Presence::new(
            3,
            PresenceOpts {
                kind: PresenceKind::SpikeSlab,
                rho: 0.05,
                warmup: 0,
                period: 1,
            },
        );
        p.update(&counts, reads.len(), |r, buf| {
            buf.extend_from_slice(&reads[r]);
            1.0
        });
        assert!(p.q[0] > 0.99 && p.q[2] > 0.99, "{:?}", p.q);
        assert!(p.q[1] < 0.5, "{:?}", p.q);
    }

    /// Bayes factor: a transcript with exclusive reads is present; a nested one whose reads are
    /// all shared, at no more than its share of them, falls below one half; and one carrying a
    /// large excess of shared reads over what its sibling explains is present.
    #[test]
    fn bayes_factor_separates_supported_from_flat() {
        // t0: 400 exclusive reads + 100 shared with t1 (t1's reads fit t0 equally well).
        // t2: 300 reads shared with t3 that fit t2 4x better, and 30 exclusive to t3.
        let mut off = vec![0usize];
        let mut target = Vec::new();
        let mut wts = Vec::new();
        let mut push = |alns: &[(u32, f64)]| {
            for &(t, w) in alns {
                target.push(t);
                wts.push(w);
            }
            off.push(target.len());
        };
        for _ in 0..400 {
            push(&[(0, 1.0)]);
        }
        for _ in 0..100 {
            push(&[(0, 1.0), (1, 1.0)]);
        }
        for _ in 0..300 {
            push(&[(2, 1.0), (3, 0.25)]);
        }
        for _ in 0..30 {
            push(&[(3, 1.0)]);
        }
        let counts = [480.0, 20.0, 300.0, 30.0];
        let mut p = Presence::new(
            4,
            PresenceOpts {
                kind: PresenceKind::Bayes,
                rho: 0.05,
                warmup: 0,
                period: 1,
            },
        );
        let reads = ReadAlns {
            off: &off,
            target: &target,
        };
        p.refresh(&counts, reads, |k| wts[k], |_| 1.0);
        assert!(p.q[0] == 1.0 && p.q[2] == 1.0 && p.q[3] == 1.0, "{:?}", p.q);
        assert!(p.q[1] < 0.5, "{:?} pi {}", p.q, p.pi);
    }
}
