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

use rayon::prelude::*;

/// Options of the presence model.
#[derive(Clone, Copy, Debug)]
pub struct PresenceOpts {
    pub rho: f64,
    pub warmup: u32,
    pub period: u32,
}

/// Per-read clamp on the log-evidence term (a single unique read saturates here), overall clamp
/// on `delta` before the sigmoid, and the floor that keeps the ratio finite.
const READ_LOG_CAP: f64 = 50.0;
const DELTA_CAP: f64 = 40.0;
const EPS: f64 = 1e-100;

pub struct Presence {
    pub q: Vec<f64>,
    logit_rho: f64,
    warmup: u32,
    period: u32,
}

impl Presence {
    pub fn new(ntx: usize, opts: PresenceOpts) -> Self {
        let p = opts.rho.clamp(1e-4, 1.0 - 1e-4);
        Presence {
            q: vec![1.0; ntx],
            logit_rho: (p / (1.0 - p)).ln(),
            warmup: opts.warmup,
            period: opts.period.max(1),
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

    /// Refresh `q` from the leave-one-out evidence of `nreads` reads at abundances `counts`.
    /// `read(r, buf)` fills `buf` with read `r`'s (target, weight) pairs and returns its
    /// multiplicity (0 skips the read).
    pub fn update<F>(&mut self, counts: &[f64], nreads: usize, read: F)
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
}
