//! The EM on a component-decomposed layout (after salmon's component engine).
//!
//! The class-transcript graph splits into connected components that the EM never couples: a
//! class's update reads and writes only its own component's abundances. Each component's
//! classes are stored together with component-local transcript labels, so a component's data
//! stays in cache while it runs a tile of iterations. Single-target classes are folded into a
//! per-transcript constant (their reads go to their one target whatever the abundances), and
//! transcripts no class references are never visited (every round gives them 0).
//!
//! A class is an equivalence class (with a read count) or a single read (count 1); its entries
//! are (target, weight) pairs.
//!
//! The stopping rule stays global and exact: components advance in tiles of up to [`TILE`]
//! rounds, recording each round's convergence statistic, and if the global criterion first
//! holds partway through a tile, every component replays the tile from its start up to that
//! round. The rounds are therefore those of the plain loop, up to summation order.

use rayon::prelude::*;

use crate::util::constants;

/// Rounds a component runs before the global criterion is checked.
const TILE: u32 = 64;
/// Components with at least this many incidences step in parallel; smaller ones run serially,
/// in parallel with each other.
const LARGE: usize = 1 << 16;
/// A large component's classes are split into at most this many blocks, each accumulating into
/// its own buffer; the buffers are summed in block order, so the result does not depend on the
/// number of threads.
const BLOCKS: usize = 64;

/// Whether the plain EMs run on the component layout: yes unless `OARFISH_EM_COMPONENTS=0`.
pub fn enabled() -> bool {
    std::env::var("OARFISH_EM_COMPONENTS").map_or(true, |v| v != "0")
}

/// Read counts per class.
#[derive(Clone, Copy)]
pub enum Counts<'a> {
    /// Every class is one read.
    Unit,
    Of(&'a [f64]),
    /// Per-class multiplicities (bootstrap resamples).
    OfU32(&'a [u32]),
}

impl Counts<'_> {
    #[inline(always)]
    fn get(&self, c: u32) -> f64 {
        match self {
            Counts::Unit => 1.0,
            Counts::Of(n) => n[c as usize],
            Counts::OfU32(n) => n[c as usize] as f64,
        }
    }
}

/// One component: transcripts `tx` (global ids; the local id is the index), and classes
/// `off[k]..off[k + 1]` over `lab` (local ids) and `w` (weights), with class ids `cls[k]`.
/// Large components' classes are grouped into blocks `blk[b]..blk[b + 1]`.
struct Comp<W> {
    tx: Vec<u32>,
    off: Vec<u32>,
    lab: Vec<u32>,
    w: Vec<W>,
    cls: Vec<u32>,
    blk: Vec<u32>,
}

impl<W> Comp<W> {
    fn incidences(&self) -> usize {
        self.lab.len()
    }

    fn is_large(&self) -> bool {
        self.incidences() >= LARGE
    }
}

/// The component layout of a set of classes over `ntx` transcripts.
pub struct Plan<W> {
    comps: Vec<Comp<W>>,
    /// Single-target classes: (global transcript, class id).
    singles: Vec<(u32, u32)>,
    ntx: usize,
}

fn find(p: &mut [u32], mut x: u32) -> u32 {
    while p[x as usize] != x {
        p[x as usize] = p[p[x as usize] as usize];
        x = p[x as usize];
    }
    x
}

impl<W: Copy + Default + Into<f64> + Send + Sync> Plan<W> {
    /// Classes `c` with entries `offsets[c]..offsets[c + 1]`, entry `k` having target
    /// `target(k)` and weight `w(k)`. With `by_heaviest`, each component's classes are ordered
    /// by their heaviest target, so the classes of one transcript sit together and a round reads
    /// and writes nearby abundances; without it they keep their order (for classes already
    /// arranged so, as the store's reads are when `bulk` reorders them).
    pub fn new(
        target: impl Fn(usize) -> u32 + Sync,
        w: impl Fn(usize) -> W + Sync,
        offsets: &[usize],
        ntx: usize,
        by_heaviest: bool,
    ) -> Self {
        let ncls = offsets.len() - 1;
        let class = |c: usize| offsets[c]..offsets[c + 1];
        let multi = |c: usize| offsets[c + 1] - offsets[c] >= 2;
        let mut p: Vec<u32> = (0..ntx as u32).collect();
        let mut in_multi = vec![false; ntx];
        for c in (0..ncls).filter(|&c| multi(c)) {
            let mut a = find(&mut p, target(offsets[c]));
            for k in class(c) {
                let t = target(k);
                in_multi[t as usize] = true;
                let b = find(&mut p, t);
                if a != b {
                    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                    p[hi as usize] = lo;
                    a = lo;
                }
            }
        }
        // Components of the transcripts that multi-target classes reference, numbered by their
        // smallest transcript; local ids ascend with global ids.
        let mut comp_of_root = vec![u32::MAX; ntx];
        let mut comp_of = vec![u32::MAX; ntx];
        let mut local = vec![u32::MAX; ntx];
        let mut ntx_of: Vec<u32> = Vec::new();
        for t in 0..ntx {
            if !in_multi[t] {
                continue;
            }
            let r = find(&mut p, t as u32) as usize;
            if comp_of_root[r] == u32::MAX {
                comp_of_root[r] = ntx_of.len() as u32;
                ntx_of.push(0);
            }
            let k = comp_of_root[r];
            comp_of[t] = k;
            local[t] = ntx_of[k as usize];
            ntx_of[k as usize] += 1;
        }
        drop((p, comp_of_root, in_multi));
        // Size every component exactly before filling it.
        let ncomp = ntx_of.len();
        let (mut ncls_of, mut ninc_of) = (vec![0_usize; ncomp], vec![0_usize; ncomp]);
        let mut nsingle = 0;
        for c in 0..ncls {
            if multi(c) {
                let k = comp_of[target(offsets[c]) as usize] as usize;
                ncls_of[k] += 1;
                ninc_of[k] += offsets[c + 1] - offsets[c];
            } else if offsets[c + 1] > offsets[c] {
                nsingle += 1;
            }
        }
        let mut comps: Vec<Comp<W>> = (0..ncomp)
            .map(|k| Comp {
                tx: Vec::with_capacity(ntx_of[k] as usize),
                off: {
                    let mut v = Vec::with_capacity(ncls_of[k] + 1);
                    v.push(0);
                    v
                },
                lab: Vec::with_capacity(ninc_of[k]),
                w: Vec::with_capacity(ninc_of[k]),
                cls: Vec::with_capacity(ncls_of[k]),
                blk: Vec::new(),
            })
            .collect();
        for t in 0..ntx {
            if comp_of[t] != u32::MAX {
                comps[comp_of[t] as usize].tx.push(t as u32);
            }
        }
        // Each component's classes, in class order or by their heaviest target (then by class
        // id), packed as (heaviest target's local id) << 32 | class id.
        let mut start = vec![0_usize; ncomp + 1];
        for k in 0..ncomp {
            start[k + 1] = start[k] + ncls_of[k];
        }
        let mut order = vec![0_u64; start[ncomp]];
        let mut singles = Vec::with_capacity(nsingle);
        {
            let mut at = start[..ncomp].to_vec();
            for c in 0..ncls {
                if multi(c) {
                    let k = comp_of[target(offsets[c]) as usize] as usize;
                    order[at[k]] = c as u64;
                    at[k] += 1;
                } else if offsets[c + 1] > offsets[c] {
                    singles.push((target(offsets[c]), c as u32));
                }
            }
        }
        if by_heaviest {
            order.par_iter_mut().with_min_len(4096).for_each(|x| {
                let c = *x as usize;
                let heaviest = class(c)
                    .max_by(|&a, &b| w(a).into().total_cmp(&w(b).into()))
                    .expect("a multi-target class");
                *x |= (local[target(heaviest) as usize] as u64) << 32;
            });
            let mut rest = &mut order[..];
            let mut parts = Vec::with_capacity(ncomp);
            for &n in &ncls_of {
                let (a, b) = rest.split_at_mut(n);
                parts.push(a);
                rest = b;
            }
            parts.par_iter_mut().for_each(|p| p.sort_unstable());
        }
        // Fill each component, a large one in parallel over blocks of whole classes.
        comps.par_iter_mut().enumerate().for_each(|(k, cm)| {
            cm.cls
                .extend(order[start[k]..start[k + 1]].iter().map(|&x| x as u32));
            assert!(
                ninc_of[k] < u32::MAX as usize,
                "a component with 2^32 incidences"
            );
            let mut at = 0_u32;
            for &c in &cm.cls {
                at += (offsets[c as usize + 1] - offsets[c as usize]) as u32;
                cm.off.push(at);
            }
            cm.lab.resize(at as usize, 0);
            cm.w.resize(at as usize, W::default());
            const FILL: usize = 1 << 16;
            let mut parts = Vec::new();
            let (mut lab, mut wr) = (&mut cm.lab[..], &mut cm.w[..]);
            let mut c0 = 0;
            while c0 < cm.cls.len() {
                let mut c1 = c0;
                while c1 < cm.cls.len() && ((cm.off[c1] - cm.off[c0]) as usize) < FILL {
                    c1 += 1;
                }
                let len = (cm.off[c1] - cm.off[c0]) as usize;
                let (l, lr) = lab.split_at_mut(len);
                let (x, xr) = wr.split_at_mut(len);
                parts.push((c0..c1, l, x));
                (lab, wr) = (lr, xr);
                c0 = c1;
            }
            let cls = &cm.cls;
            parts.into_par_iter().for_each(|(r, l, x)| {
                let mut j = 0;
                for &c in &cls[r] {
                    for e in class(c as usize) {
                        l[j] = local[target(e) as usize];
                        x[j] = w(e);
                        j += 1;
                    }
                }
            });
        });
        drop(order);
        for cm in comps.iter_mut().filter(|cm| cm.is_large()) {
            // Blocks of about equal incidence counts, on class boundaries.
            let per = cm.incidences().div_ceil(BLOCKS);
            cm.blk.push(0);
            for k in 1..cm.cls.len() {
                if cm.off[k] as usize >= per * cm.blk.len() {
                    cm.blk.push(k as u32);
                }
            }
            cm.blk.push(cm.cls.len() as u32);
        }
        Plan {
            comps,
            singles,
            ntx,
        }
    }

    pub fn n_components(&self) -> usize {
        self.comps.len()
    }

    /// Whether some component is stepped in parallel.
    #[cfg(test)]
    pub fn has_large(&self) -> bool {
        self.comps.iter().any(|cm| cm.is_large())
    }
}

impl<W: Copy + Into<f64> + Send + Sync> Comp<W> {
    /// Adds class `k`'s round to `acc`: each of its targets gets `w * n / denom`.
    #[inline(always)]
    fn push(&self, k: usize, n: f64, prev: &[f64], acc: &mut [f64]) {
        if n == 0.0 {
            return;
        }
        let r = self.off[k] as usize..self.off[k + 1] as usize;
        let (lab, w) = (&self.lab[r.clone()], &self.w[r]);
        // SAFETY: every label is below the component's transcript count (`Plan::new`), and
        // `prev` and `acc` have that length.
        let d: f64 = lab
            .iter()
            .zip(w)
            .map(|(&l, &w)| unsafe { *prev.get_unchecked(l as usize) } * w.into())
            .sum();
        if d > constants::EM_DENOM_THRESH {
            let s = n / d;
            for (&l, &w) in lab.iter().zip(w) {
                unsafe { *acc.get_unchecked_mut(l as usize) += w.into() * s };
            }
        }
    }

    /// One EM round: `out[l] = base[l] + prev[l] * sum over l's classes of w * n / denom`.
    /// Returns the round's convergence statistic (largest relative increase over transcripts
    /// above [`constants::MIN_READ_THRESH`]). `acc` is scratch: one buffer of the component's
    /// size per block (one for a small component).
    fn step(&self, n: Counts, base: &[f64], prev: &[f64], out: &mut [f64], acc: &mut [f64]) -> f64 {
        self.step_by(n, base, prev, out, acc, rel_increase)
    }

    /// [`Comp::step`], with the convergence statistic the largest `stat(prev, out)`.
    fn step_by(
        &self,
        n: Counts,
        base: &[f64],
        prev: &[f64],
        out: &mut [f64],
        acc: &mut [f64],
        stat: impl Fn(f64, f64) -> f64 + Sync + Copy,
    ) -> f64 {
        let m = self.tx.len();
        if !self.is_large() {
            acc.fill(0.0);
            for k in 0..self.cls.len() {
                self.push(k, n.get(self.cls[k]), prev, acc);
            }
            let mut rd = 0.0_f64;
            for l in 0..m {
                out[l] = base[l] + prev[l] * acc[l];
                rd = rd.max(stat(prev[l], out[l]));
            }
            return rd;
        }
        acc.par_chunks_mut(m).enumerate().for_each(|(b, a)| {
            a.fill(0.0);
            for k in self.blk[b] as usize..self.blk[b + 1] as usize {
                self.push(k, n.get(self.cls[k]), prev, a);
            }
        });
        let acc = &*acc;
        out.par_iter_mut()
            .with_min_len(256)
            .enumerate()
            .map(|(l, o)| {
                let s: f64 = acc.chunks(m).map(|a| a[l]).sum();
                *o = base[l] + prev[l] * s;
                stat(prev[l], *o)
            })
            .reduce(|| 0.0, f64::max)
    }
}

/// The relative change in either direction, over abundances above
/// [`constants::MIN_READ_THRESH`] (the truncation EM's statistic).
#[inline]
fn rel_change(prev: f64, curr: f64) -> f64 {
    if prev > constants::MIN_READ_THRESH {
        (curr - prev).abs() / prev
    } else {
        0.0
    }
}

#[inline]
fn rel_increase(prev: f64, curr: f64) -> f64 {
    if prev > constants::MIN_READ_THRESH {
        (curr - prev) / prev
    } else {
        0.0
    }
}

/// A component's state: its single-target constant, abundances, scratch, and the convergence
/// statistic of each round of the current tile.
struct State {
    base: Vec<f64>,
    prev: Vec<f64>,
    curr: Vec<f64>,
    acc: Vec<f64>,
    rd: Vec<f64>,
}

/// The plain EM (no presence model) with class read counts `n`, from abundances `init`: the
/// rounds, global stopping rule and final pruning round of [`crate::eq_classes`] and
/// [`crate::em_read`]. The loop stops after the first round past `min_iter` whose
/// statistic is below `thresh`, or after `max_iter` rounds.
pub fn run<W: Copy + Into<f64> + Send + Sync>(
    plan: &Plan<W>,
    n: Counts,
    init: &[f64],
    max_iter: u32,
    thresh: f64,
    min_iter: u32,
) -> Vec<f64> {
    let mut single = vec![0.0_f64; plan.ntx];
    for &(t, c) in &plan.singles {
        single[t as usize] += n.get(c);
    }
    let mut st: Vec<State> = plan
        .comps
        .iter()
        .map(|cm| State {
            base: cm.tx.iter().map(|&t| single[t as usize]).collect(),
            prev: cm.tx.iter().map(|&t| init[t as usize]).collect(),
            curr: vec![0.0; cm.tx.len()],
            acc: vec![0.0; cm.tx.len() * cm.blk.len().saturating_sub(1).max(1)],
            rd: Vec::with_capacity(TILE as usize),
        })
        .collect();
    // Transcripts that only single-target classes reference: every round gives them their
    // reads, so their statistic is (reads - init) / init in the first round and 0 after.
    let mut in_comp = vec![false; plan.ntx];
    for cm in &plan.comps {
        for &t in &cm.tx {
            in_comp[t as usize] = true;
        }
    }
    let lone: Vec<u32> = (0..plan.ntx as u32)
        .filter(|&t| single[t as usize] > 0.0 && !in_comp[t as usize])
        .collect();
    drop(in_comp);
    let lone_rd0 = lone
        .iter()
        .map(|&t| rel_increase(init[t as usize], single[t as usize]))
        .fold(0.0_f64, f64::max);

    // Largest components first, so the long ones start before the many short ones.
    let mut order: Vec<usize> = (0..plan.comps.len()).collect();
    order.sort_by_key(|&k| std::cmp::Reverse(plan.comps[k].incidences()));
    let advance = |st: &mut [State], steps: u32| {
        let mut refs: Vec<Option<&mut State>> = st.iter_mut().map(Some).collect();
        let mut jobs: Vec<(usize, &mut State)> = order
            .iter()
            .map(|&k| (k, refs[k].take().expect("each component once")))
            .collect();
        jobs.par_iter_mut().with_max_len(1).for_each(|(k, s)| {
            let cm = &plan.comps[*k];
            s.rd.clear();
            for _ in 0..steps {
                let r = cm.step(n, &s.base, &s.prev, &mut s.curr, &mut s.acc);
                s.rd.push(r);
                std::mem::swap(&mut s.prev, &mut s.curr);
            }
        });
    };

    let mut niter = 0_u32;
    while niter < max_iter {
        let steps = TILE.min(max_iter - niter);
        let snap: Vec<Vec<f64>> = st.iter().map(|s| s.prev.clone()).collect();
        advance(&mut st, steps);
        let stop = (0..steps as usize).find(|&i| {
            let mut rd = st.iter().map(|s| s.rd[i]).fold(0.0_f64, f64::max);
            if niter == 0 && i == 0 {
                rd = rd.max(lone_rd0);
            }
            rd < thresh && niter + i as u32 > min_iter
        });
        match stop {
            Some(i) => {
                // The global criterion first holds after round `niter + i`: replay to there.
                for (s, p) in st.iter_mut().zip(snap) {
                    s.prev = p;
                }
                advance(&mut st, i as u32 + 1);
                break;
            }
            None => niter += steps,
        }
    }

    // Zero very small abundances, then one more round.
    for s in st.iter_mut() {
        for x in s.prev.iter_mut() {
            if *x < constants::MIN_READ_THRESH {
                *x = 0.0;
            }
        }
    }
    advance(&mut st, 1);
    let mut out = vec![0.0_f64; plan.ntx];
    for (cm, s) in plan.comps.iter().zip(&st) {
        for (&t, &x) in cm.tx.iter().zip(&s.prev) {
            out[t as usize] = x;
        }
    }
    for &t in &lone {
        out[t as usize] = single[t as usize];
    }
    out
}

impl<W: Copy + Into<f64> + Send + Sync> Comp<W> {
    /// The component's log-likelihood at abundances `theta`, up to a constant: over its
    /// classes (read counts `n`) and its single-target reads `base`, with `theta` normalized to
    /// proportions, as the truncation EM's `loglik_static` computes it over all reads.
    fn loglik(&self, n: Counts, base: &[f64], theta: &[f64]) -> f64 {
        let class_ll = |k: usize| -> (f64, f64) {
            let c = n.get(self.cls[k]);
            let r = self.off[k] as usize..self.off[k + 1] as usize;
            let d: f64 = self.lab[r.clone()]
                .iter()
                .zip(&self.w[r])
                .map(|(&l, &w)| theta[l as usize] * w.into())
                .sum();
            if c > 0.0 && d > constants::EM_DENOM_THRESH {
                (c * d.ln(), c)
            } else {
                (0.0, 0.0)
            }
        };
        let (mut ll, mut reads) = if self.is_large() {
            let parts: Vec<(f64, f64)> = (0..self.blk.len() - 1)
                .into_par_iter()
                .map(|b| {
                    (self.blk[b] as usize..self.blk[b + 1] as usize)
                        .map(class_ll)
                        .fold((0.0, 0.0), |a, x| (a.0 + x.0, a.1 + x.1))
                })
                .collect();
            parts
                .into_iter()
                .fold((0.0, 0.0), |a, x| (a.0 + x.0, a.1 + x.1))
        } else {
            (0..self.cls.len())
                .map(class_ll)
                .fold((0.0, 0.0), |a, x| (a.0 + x.0, a.1 + x.1))
        };
        for (&s, &t) in base.iter().zip(theta) {
            if s > 0.0 && t > constants::EM_DENOM_THRESH {
                ll += s * t.ln();
                reads += s;
            }
        }
        let total: f64 = theta.iter().sum();
        ll - reads * total.ln()
    }
}

/// A component's SQUAREM state (scheme S3, as in the truncation EM): the iterates of a cycle,
/// the bound on the extrapolation step length, and the cycles that extrapolated and kept the
/// result, rejected it, or did not extrapolate.
#[derive(Clone)]
struct Sq {
    step_max: f64,
    acc: u32,
    rej: u32,
    plain: u32,
}

const STEP_GROW: f64 = 4.0;

/// SQUAREM totals of [`run_squarem`] over the components, and the iteration count.
pub struct SquaremReport {
    pub niter: u32,
    pub acc: u64,
    pub rej: u64,
    pub plain: u64,
}

struct SqState {
    base: Vec<f64>,
    prev: Vec<f64>,
    curr: Vec<f64>,
    acc: Vec<f64>,
    t0: Vec<f64>,
    t1: Vec<f64>,
    t2: Vec<f64>,
    sq: Sq,
    rd: Vec<f64>,
}

impl<W: Copy + Into<f64> + Send + Sync> Comp<W> {
    /// One SQUAREM cycle from `s.prev`, as the truncation EM's `squarem_cycle` runs over all
    /// transcripts, on this component alone: two EM steps, a step along the extrapolated path,
    /// one EM step from there, kept if the component's likelihood is at least that of plain EM.
    /// Leaves the result in `s.prev` and returns the last EM step's statistic.
    fn squarem_cycle(&self, n: Counts, s: &mut SqState) -> f64 {
        s.t0.copy_from_slice(&s.prev);
        self.step_by(n, &s.base, &s.t0, &mut s.t1, &mut s.acc, rel_change);
        let rd2 = self.step_by(n, &s.base, &s.t1, &mut s.t2, &mut s.acc, rel_change);
        let (mut rr, mut vv) = (0.0_f64, 0.0_f64);
        for i in 0..s.t0.len() {
            let (a, b, c) = (s.t0[i], s.t1[i], s.t2[i]);
            let (r, v) = (b - a, c - 2.0 * b + a);
            rr += r * r;
            vv += v * v;
        }
        let mut st = if vv > 0.0 { (rr / vv).sqrt() } else { 1.0 };
        if !st.is_finite() || st < 1.0 {
            st = 1.0;
        }
        let hit = st > 1.0 && st >= s.sq.step_max;
        st = st.min(s.sq.step_max);
        if st == 1.0 {
            if hit {
                s.sq.step_max *= STEP_GROW;
            }
            s.sq.plain += 1;
            s.prev.copy_from_slice(&s.t2);
            return rd2;
        }
        for i in 0..s.t0.len() {
            let (a, b, c) = (s.t0[i], s.t1[i], s.t2[i]);
            let (r, v) = (b - a, c - 2.0 * b + a);
            s.prev[i] = (a + 2.0 * st * r + st * st * v).max(0.0);
        }
        let rd3 = self.step_by(n, &s.base, &s.prev, &mut s.curr, &mut s.acc, rel_change);
        let (l3, l2) = (
            self.loglik(n, &s.base, &s.curr),
            self.loglik(n, &s.base, &s.t2),
        );
        if l3.is_finite() && l3 >= l2 {
            s.sq.acc += 1;
            if hit {
                s.sq.step_max *= STEP_GROW;
            }
            std::mem::swap(&mut s.prev, &mut s.curr);
            rd3
        } else {
            s.sq.rej += 1;
            s.sq.step_max = (s.sq.step_max / STEP_GROW).max(1.0);
            s.prev.copy_from_slice(&s.t2);
            rd2
        }
    }
}

/// The truncation EM's frozen-weight loop on the component layout, with SQUAREM run on each
/// component separately: from iteration `niter0` and abundances `init`, plain EM steps while
/// the iteration count is below `plain_first`, then SQUAREM cycles (each counted as three
/// iterations), stopping at the first iteration past 1 where the largest relative change in
/// either direction, over all transcripts, is below `thresh`, or at `max_iter`; then very small
/// abundances are zeroed and one more EM step is taken.
///
/// Components are independent fixed-point problems, so a SQUAREM per component is SQUAREM on
/// each; the step lengths and acceptance tests are each component's own, where the loop over
/// all transcripts uses one step length and the total likelihood.
#[allow(clippy::too_many_arguments)]
pub fn run_squarem<W: Copy + Into<f64> + Send + Sync>(
    plan: &Plan<W>,
    n: Counts,
    init: &[f64],
    niter0: u32,
    max_iter: u32,
    thresh: f64,
    plain_first: u32,
) -> (Vec<f64>, SquaremReport) {
    let mut single = vec![0.0_f64; plan.ntx];
    for &(t, c) in &plan.singles {
        single[t as usize] += n.get(c);
    }
    let mut st: Vec<SqState> = plan
        .comps
        .iter()
        .map(|cm| {
            let m = cm.tx.len();
            SqState {
                base: cm.tx.iter().map(|&t| single[t as usize]).collect(),
                prev: cm.tx.iter().map(|&t| init[t as usize]).collect(),
                curr: vec![0.0; m],
                acc: vec![0.0; m * cm.blk.len().saturating_sub(1).max(1)],
                t0: vec![0.0; m],
                t1: vec![0.0; m],
                t2: vec![0.0; m],
                sq: Sq {
                    step_max: 1.0,
                    acc: 0,
                    rej: 0,
                    plain: 0,
                },
                rd: Vec::with_capacity(TILE as usize),
            }
        })
        .collect();
    let mut in_comp = vec![false; plan.ntx];
    for cm in &plan.comps {
        for &t in &cm.tx {
            in_comp[t as usize] = true;
        }
    }
    // Transcripts only single-target classes reference reach their reads in one step and stay.
    let lone: Vec<u32> = (0..plan.ntx as u32)
        .filter(|&t| single[t as usize] > 0.0 && !in_comp[t as usize])
        .collect();
    drop(in_comp);
    let lone_rd0 = lone
        .iter()
        .map(|&t| rel_change(init[t as usize], single[t as usize]))
        .fold(0.0_f64, f64::max);

    // The iterations of the loop: (count before it, whether it is a SQUAREM cycle).
    let kind = |niter: u32| niter >= plain_first;
    let mut order: Vec<usize> = (0..plan.comps.len()).collect();
    order.sort_by_key(|&k| std::cmp::Reverse(plan.comps[k].incidences()));
    let advance = |st: &mut [SqState], from: u32, loops: u32| {
        let mut refs: Vec<Option<&mut SqState>> = st.iter_mut().map(Some).collect();
        let mut jobs: Vec<(usize, &mut SqState)> = order
            .iter()
            .map(|&k| (k, refs[k].take().expect("each component once")))
            .collect();
        jobs.par_iter_mut().with_max_len(1).for_each(|(k, s)| {
            let cm = &plan.comps[*k];
            s.rd.clear();
            let mut niter = from;
            for _ in 0..loops {
                let r = if kind(niter) {
                    niter += 3;
                    cm.squarem_cycle(n, s)
                } else {
                    niter += 1;
                    let r = cm.step_by(n, &s.base, &s.prev, &mut s.curr, &mut s.acc, rel_change);
                    std::mem::swap(&mut s.prev, &mut s.curr);
                    r
                };
                s.rd.push(r);
            }
        });
    };

    let mut niter = niter0;
    let mut first = true;
    'outer: while niter < max_iter {
        // Up to TILE loop iterations, all starting below max_iter.
        let mut starts = Vec::with_capacity(TILE as usize);
        let mut x = niter;
        while starts.len() < TILE as usize && x < max_iter {
            starts.push(x);
            x += if kind(x) { 3 } else { 1 };
        }
        let snap: Vec<(Vec<f64>, Sq)> = st.iter().map(|s| (s.prev.clone(), s.sq.clone())).collect();
        advance(&mut st, niter, starts.len() as u32);
        for (i, &x0) in starts.iter().enumerate() {
            let mut rd = st.iter().map(|s| s.rd[i]).fold(0.0_f64, f64::max);
            if first && i == 0 {
                rd = rd.max(lone_rd0);
            }
            // The loop's count at its test: after a cycle's two extra iterations.
            let at = if kind(x0) { x0 + 2 } else { x0 };
            if rd < thresh && at > 1 {
                for (s, (p, q)) in st.iter_mut().zip(snap) {
                    s.prev = p;
                    s.sq = q;
                }
                advance(&mut st, niter, i as u32 + 1);
                niter = at;
                break 'outer;
            }
        }
        first = false;
        niter = x;
    }

    let report = SquaremReport {
        niter,
        acc: st.iter().map(|s| s.sq.acc as u64).sum(),
        rej: st.iter().map(|s| s.sq.rej as u64).sum(),
        plain: st.iter().map(|s| s.sq.plain as u64).sum(),
    };
    // Zero very small abundances, then one more EM step.
    let mut out = vec![0.0_f64; plan.ntx];
    plan.comps
        .par_iter()
        .zip(st.par_iter_mut())
        .for_each(|(cm, s)| {
            for x in s.prev.iter_mut() {
                if *x < constants::MIN_READ_THRESH {
                    *x = 0.0;
                }
            }
            cm.step_by(n, &s.base, &s.prev, &mut s.curr, &mut s.acc, rel_change);
        });
    for (cm, s) in plan.comps.iter().zip(&st) {
        for (&t, &x) in cm.tx.iter().zip(&s.curr) {
            out[t as usize] = x;
        }
    }
    for &t in &lone {
        out[t as usize] = single[t as usize];
    }
    (out, report)
}
