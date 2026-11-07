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
}

impl Counts<'_> {
    #[inline(always)]
    fn get(&self, c: u32) -> f64 {
        match self {
            Counts::Unit => 1.0,
            Counts::Of(n) => n[c as usize],
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

impl<W: Copy + Into<f64> + Send + Sync> Plan<W> {
    /// Classes `c` with entries `offsets[c]..offsets[c + 1]`, entry `k` having target
    /// `target(k)` and weight `w(k)`.
    pub fn new(
        target: impl Fn(usize) -> u32 + Sync,
        w: impl Fn(usize) -> W + Sync,
        offsets: &[usize],
        ntx: usize,
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
        // Each component's classes, ordered by their heaviest target (then by class id): the
        // classes of one transcript sit together, so a round reads and writes nearby
        // abundances. Packed as (heaviest target's local id) << 32 | class id.
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
        order.par_iter_mut().with_min_len(4096).for_each(|x| {
            let c = *x as usize;
            let heaviest = class(c)
                .max_by(|&a, &b| w(a).into().total_cmp(&w(b).into()))
                .expect("a multi-target class");
            *x |= (local[target(heaviest) as usize] as u64) << 32;
        });
        {
            let mut rest = &mut order[..];
            let mut parts = Vec::with_capacity(ncomp);
            for &n in &ncls_of {
                let (a, b) = rest.split_at_mut(n);
                parts.push(a);
                rest = b;
            }
            parts.par_iter_mut().for_each(|p| p.sort_unstable());
        }
        comps.par_iter_mut().enumerate().for_each(|(k, cm)| {
            for &x in &order[start[k]..start[k + 1]] {
                let c = (x & u32::MAX as u64) as usize;
                for e in class(c) {
                    cm.lab.push(local[target(e) as usize]);
                    cm.w.push(w(e));
                }
                cm.off.push(cm.lab.len() as u32);
                cm.cls.push(c as u32);
            }
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
        let m = self.tx.len();
        if !self.is_large() {
            acc.fill(0.0);
            for k in 0..self.cls.len() {
                self.push(k, n.get(self.cls[k]), prev, acc);
            }
            let mut rd = 0.0_f64;
            for l in 0..m {
                out[l] = base[l] + prev[l] * acc[l];
                rd = rd.max(rel_increase(prev[l], out[l]));
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
                rel_increase(prev[l], *o)
            })
            .reduce(|| 0.0, f64::max)
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
/// [`crate::em_read`].
pub fn run<W: Copy + Into<f64> + Send + Sync>(
    plan: &Plan<W>,
    n: Counts,
    init: &[f64],
    max_iter: u32,
    thresh: f64,
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
            rd < thresh && niter + i as u32 > 1
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
