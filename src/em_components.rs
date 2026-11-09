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

use std::sync::atomic::{AtomicBool, AtomicU32, Ordering::Relaxed};

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
// Block ids are stored as `u8`.
const _: () = assert!(BLOCKS <= 256);
/// Entry counts at or above this sort as one size in [`Plan::new_by_size`].
const SIZE_CAP: usize = 255;

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
/// Large components' classes are grouped into blocks `blk[b]..blk[b + 1]`; block `b` touches
/// the labels `touch[touch_off[b]..touch_off[b + 1]]` (ascending), and label `l` is touched by
/// the blocks `by[by_off[l]..by_off[l + 1]]` (ascending).
struct Comp<W> {
    tx: Vec<u32>,
    off: Vec<u32>,
    lab: Vec<u32>,
    w: Vec<W>,
    cls: Vec<u32>,
    blk: Vec<u32>,
    touch_off: Vec<usize>,
    touch: Vec<u32>,
    by_off: Vec<u32>,
    by: Vec<u8>,
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

/// The root of `x` in the shared forest `p`, halving the path on the way. A parent only ever
/// moves to a smaller id that is an ancestor, so the forest stays acyclic under any interleaving.
fn find_shared(p: &[AtomicU32], mut x: u32) -> u32 {
    loop {
        let px = p[x as usize].load(Relaxed);
        if px == x {
            return x;
        }
        let gp = p[px as usize].load(Relaxed);
        if gp != px {
            let _ = p[x as usize].compare_exchange_weak(px, gp, Relaxed, Relaxed);
        }
        x = gp;
    }
}

/// Joins the trees of `a` and `b` in the shared forest `p`, always linking the larger root under
/// the smaller, so every root is its tree's smallest id whatever the order of the joins. Returns
/// a member of the joined tree.
fn union_shared(p: &[AtomicU32], a: u32, b: u32) -> u32 {
    let (mut a, mut b) = (a, b);
    loop {
        a = find_shared(p, a);
        b = find_shared(p, b);
        if a == b {
            return a;
        }
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        // Fails only if `hi` stopped being a root meanwhile; then retry from the new roots.
        if p[hi as usize]
            .compare_exchange(hi, lo, Relaxed, Relaxed)
            .is_ok()
        {
            return lo;
        }
        (a, b) = (lo, hi);
    }
}

/// A buffer that parallel tasks write at disjoint positions.
struct Disjoint<T>(*mut T);
// SAFETY: used only for writes to positions no two tasks share.
unsafe impl<T: Send> Sync for Disjoint<T> {}
impl<T> Disjoint<T> {
    /// # Safety
    /// `i` is in bounds and no other task writes it.
    #[inline(always)]
    unsafe fn write(&self, i: usize, x: T) {
        unsafe { self.0.add(i).write(x) }
    }
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
        Self::build(target, w, offsets, ntx, by_heaviest, false)
    }

    /// As [`Plan::new`], but each component's classes are first ordered by their entry count
    /// (capped at [`SIZE_CAP`]), stably, so the order `new` would give holds among classes of one
    /// size. Classes of equal size sit together, and the variable-length loops of a round over
    /// them exit at predictable trip counts (for classes of reads, whose alignment counts vary).
    pub fn new_by_size(
        target: impl Fn(usize) -> u32 + Sync,
        w: impl Fn(usize) -> W + Sync,
        offsets: &[usize],
        ntx: usize,
        by_heaviest: bool,
    ) -> Self {
        Self::build(target, w, offsets, ntx, by_heaviest, true)
    }

    fn build(
        target: impl Fn(usize) -> u32 + Sync,
        w: impl Fn(usize) -> W + Sync,
        offsets: &[usize],
        ntx: usize,
        by_heaviest: bool,
        by_size: bool,
    ) -> Self {
        let ncls = offsets.len() - 1;
        let class = |c: usize| offsets[c]..offsets[c + 1];
        let multi = |c: usize| offsets[c + 1] - offsets[c] >= 2;
        // The class-transcript graph's components, joined concurrently; each root is its
        // component's smallest transcript, so the result is that of a serial pass.
        let (mut p, in_multi) = {
            let p: Vec<AtomicU32> = (0..ntx as u32).map(AtomicU32::new).collect();
            let in_multi: Vec<AtomicBool> = (0..ntx).map(|_| AtomicBool::new(false)).collect();
            (0..ncls)
                .into_par_iter()
                .with_min_len(4096)
                .filter(|&c| multi(c))
                .for_each(|c| {
                    let mut a = target(offsets[c]);
                    for k in class(c) {
                        let t = target(k);
                        if !in_multi[t as usize].load(Relaxed) {
                            in_multi[t as usize].store(true, Relaxed);
                        }
                        a = union_shared(&p, a, t);
                    }
                });
            let p: Vec<u32> = p.into_iter().map(AtomicU32::into_inner).collect();
            let in_multi: Vec<bool> = in_multi.into_iter().map(AtomicBool::into_inner).collect();
            (p, in_multi)
        };
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
        // Count each component's classes (and incidences) per chunk of classes, and each chunk's
        // single-target classes; the counts' running sums, in component then chunk order, place
        // each chunk's classes after those of the chunks before it, so every component keeps
        // class order. Every buffer that outlives this is allocated on the calling thread, as
        // the serial passes did (an allocator with per-thread heaps reuses that thread's freed
        // memory, not the workers').
        let ncomp = ntx_of.len();
        let comp = |c: usize| comp_of[target(offsets[c]) as usize] as usize;
        let width = ncomp.max(1);
        let nchunk = ncls
            .div_ceil(1 << 16)
            .min(4 * rayon::current_num_threads())
            .min((1 << 21) / width)
            .max(1);
        let per = ncls.div_ceil(nchunk).max(1);
        let chunk = |j: usize| (j * per).min(ncls)..((j + 1) * per).min(ncls);
        let mut at = vec![0_usize; nchunk * width];
        let mut inc = vec![0_usize; nchunk * width];
        let mut sat = vec![0_usize; nchunk];
        at.par_chunks_mut(width)
            .zip(inc.par_chunks_mut(width))
            .zip(sat.par_iter_mut())
            .enumerate()
            .for_each(|(j, ((n, i), s))| {
                for c in chunk(j) {
                    let len = offsets[c + 1] - offsets[c];
                    if len >= 2 {
                        let k = comp(c);
                        n[k] += 1;
                        i[k] += len;
                    } else if len == 1 {
                        *s += 1;
                    }
                }
            });
        let mut start = vec![0_usize; ncomp + 1];
        let (mut ncls_of, mut ninc_of) = (vec![0_usize; ncomp], vec![0_usize; ncomp]);
        for k in 0..ncomp {
            let mut run = start[k];
            for j in 0..nchunk {
                let n = at[j * width + k];
                at[j * width + k] = run;
                run += n;
                ninc_of[k] += inc[j * width + k];
            }
            ncls_of[k] = run - start[k];
            start[k + 1] = run;
        }
        drop(inc);
        let mut nsingle = 0;
        for s in sat.iter_mut() {
            (*s, nsingle) = (nsingle, nsingle + *s);
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
                touch_off: Vec::new(),
                touch: Vec::new(),
                by_off: Vec::new(),
                by: Vec::new(),
            })
            .collect();
        for t in 0..ntx {
            if comp_of[t] != u32::MAX {
                comps[comp_of[t] as usize].tx.push(t as u32);
            }
        }
        // Each component's classes, in class order or by their heaviest target (then by class
        // id), packed as (heaviest target's local id) << 32 | class id.
        let mut order = vec![0_u64; start[ncomp]];
        let mut singles = vec![(0_u32, 0_u32); nsingle];
        {
            let out = Disjoint(order.as_mut_ptr());
            let one = Disjoint(singles.as_mut_ptr());
            at.par_chunks_mut(width)
                .zip(sat.par_iter_mut())
                .enumerate()
                .for_each(|(j, (at, s))| {
                    for c in chunk(j) {
                        let len = offsets[c + 1] - offsets[c];
                        if len >= 2 {
                            let k = comp(c);
                            // SAFETY: chunk `j`'s classes of component `k` own positions
                            // `at[k]..` up to the next chunk's start, all below `start[ncomp]`.
                            unsafe { out.write(at[k], c as u64) };
                            at[k] += 1;
                        } else if len == 1 {
                            // SAFETY: as above, chunk `j`'s single-target classes own
                            // positions `*s..` up to the next chunk's, all below `nsingle`.
                            unsafe { one.write(*s, (target(offsets[c]), c as u32)) };
                            *s += 1;
                        }
                    }
                });
        }
        drop((at, sat));
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
            let part = &order[start[k]..start[k + 1]];
            if by_size {
                let len = |x: u64| {
                    let c = x as u32 as usize;
                    offsets[c + 1] - offsets[c]
                };
                // A stable counting sort by capped entry count, from `p` straight into `out`.
                let sort = |p: &[u64], out: &mut [u32]| {
                    let mut at = [0_usize; SIZE_CAP + 2];
                    for &x in p {
                        at[len(x).min(SIZE_CAP) + 1] += 1;
                    }
                    for s in 1..at.len() {
                        at[s] += at[s - 1];
                    }
                    for &x in p {
                        let s = len(x).min(SIZE_CAP);
                        out[at[s]] = x as u32;
                        at[s] += 1;
                    }
                };
                cm.cls.resize(part.len(), 0);
                if ninc_of[k] >= LARGE {
                    // Sort within runs of about the incidences of one block, so each block keeps
                    // the classes (and transcripts) the order above gives it.
                    let per = ninc_of[k].div_ceil(BLOCKS);
                    let mut segs = Vec::with_capacity(BLOCKS + 1);
                    let (mut rest, mut out) = (part, &mut cm.cls[..]);
                    let (mut run, mut c0, mut nseg) = (0_usize, 0_usize, 1_usize);
                    for (c, &x) in part.iter().enumerate() {
                        if c > 0 && run >= per * nseg {
                            let (a, b) = rest.split_at(c - c0);
                            let (o, ob) = out.split_at_mut(c - c0);
                            segs.push((a, o));
                            (rest, out, c0, nseg) = (b, ob, c, nseg + 1);
                        }
                        run += len(x);
                    }
                    segs.push((rest, out));
                    segs.into_par_iter().for_each(|(p, o)| sort(p, o));
                } else {
                    sort(part, &mut cm.cls[..]);
                }
            } else {
                cm.cls.extend(part.iter().map(|&x| x as u32));
            }
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
            // The labels each block touches, and for each label the blocks that touch it.
            let nblk = cm.blk.len() - 1;
            let m = cm.tx.len();
            let lists: Vec<Vec<u32>> = (0..nblk)
                .into_par_iter()
                .map(|b| {
                    let r = cm.off[cm.blk[b] as usize] as usize
                        ..cm.off[cm.blk[b + 1] as usize] as usize;
                    let mut seen = vec![false; m];
                    for &l in &cm.lab[r] {
                        seen[l as usize] = true;
                    }
                    (0..m as u32).filter(|&l| seen[l as usize]).collect()
                })
                .collect();
            cm.touch_off = Vec::with_capacity(nblk + 1);
            cm.touch_off.push(0);
            let mut by_off = vec![0_u32; m + 1];
            for v in &lists {
                cm.touch_off.push(cm.touch_off.last().unwrap() + v.len());
                for &l in v {
                    by_off[l as usize + 1] += 1;
                }
            }
            for l in 0..m {
                by_off[l + 1] += by_off[l];
            }
            let mut at = by_off[..m].to_vec();
            cm.by = vec![0_u8; by_off[m] as usize];
            for (b, v) in lists.iter().enumerate() {
                for &l in v {
                    cm.by[at[l as usize] as usize] = b as u8;
                    at[l as usize] += 1;
                }
            }
            cm.by_off = by_off;
            cm.touch = lists.concat();
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
    #[cfg(test)]
    fn step(&self, n: Counts, base: &[f64], prev: &[f64], out: &mut [f64], acc: &mut [f64]) -> f64 {
        self.step_by(self.count_of(n), base, prev, out, acc, rel_increase)
    }

    /// The read count of each of this component's classes, in component order.
    #[inline(always)]
    fn count_of<'a>(&'a self, n: Counts<'a>) -> impl Fn(usize) -> f64 + Sync + Copy + 'a {
        move |k| n.get(self.cls[k])
    }

    /// One EM round as `Comp::step` takes it, with the convergence statistic the largest
    /// `stat(prev, out)`, and `n(k)` the read count of the component's class `k` (local index).
    fn step_by(
        &self,
        n: impl Fn(usize) -> f64 + Sync + Copy,
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
                self.push(k, n(k), prev, acc);
            }
            let mut rd = 0.0_f64;
            for l in 0..m {
                out[l] = base[l] + prev[l] * acc[l];
                rd = rd.max(stat(prev[l], out[l]));
            }
            return rd;
        }
        // Each block zeroes and fills only the labels it touches, and each label sums only the
        // blocks that touch it, in block order: the skipped addends are the +0.0 of untouched
        // entries, so the sums are those over every block.
        acc.par_chunks_mut(m).enumerate().for_each(|(b, a)| {
            // SAFETY: the labels are below `m`, the length of `a`.
            for &l in &self.touch[self.touch_off[b]..self.touch_off[b + 1]] {
                unsafe { *a.get_unchecked_mut(l as usize) = 0.0 };
            }
            for k in self.blk[b] as usize..self.blk[b + 1] as usize {
                self.push(k, n(k), prev, a);
            }
        });
        let acc = &*acc;
        out.par_iter_mut()
            .with_min_len(256)
            .enumerate()
            .map(|(l, o)| {
                let bs = &self.by[self.by_off[l] as usize..self.by_off[l + 1] as usize];
                // SAFETY: each block id is below the block count, and `acc` holds that many
                // buffers of length `m`.
                let s: f64 = bs
                    .iter()
                    .map(|&b| unsafe { *acc.get_unchecked(b as usize * m + l) })
                    .sum();
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

/// A component's state: its single-target constant, its classes' read counts gathered into
/// component order (empty when every class is one read), abundances, scratch, and the
/// convergence statistic of each round of the current tile.
struct State {
    base: Vec<f64>,
    cnt: Vec<f64>,
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
            cnt: match n {
                Counts::Unit => Vec::new(),
                _ => cm.cls.iter().map(|&c| n.get(c)).collect(),
            },
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
            let (unit, cnt) = (matches!(n, Counts::Unit), &s.cnt[..]);
            let nk = move |k: usize| if unit { 1.0 } else { cnt[k] };
            for _ in 0..steps {
                let r = cm.step_by(nk, &s.base, &s.prev, &mut s.curr, &mut s.acc, rel_increase);
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
    #[cfg(test)]
    fn loglik(&self, n: Counts, base: &[f64], theta: &[f64]) -> f64 {
        self.loglik_pair(n, base, theta, theta).0
    }

    /// [`Comp::loglik`] at two abundance vectors `a` and `b` in one pass over the classes:
    /// each of the two is summed in the order a separate call would sum it, so both results
    /// are those of separate calls.
    fn loglik_pair(&self, n: Counts, base: &[f64], a: &[f64], b: &[f64]) -> (f64, f64) {
        type Acc = (f64, f64, f64, f64);
        let add = |x: Acc, y: Acc| (x.0 + y.0, x.1 + y.1, x.2 + y.2, x.3 + y.3);
        let class_ll = |k: usize| -> Acc {
            let c = n.get(self.cls[k]);
            let r = self.off[k] as usize..self.off[k + 1] as usize;
            let (mut da, mut db) = (0.0_f64, 0.0_f64);
            for (&l, &w) in self.lab[r.clone()].iter().zip(&self.w[r]) {
                let w: f64 = w.into();
                da += a[l as usize] * w;
                db += b[l as usize] * w;
            }
            let (la, ra) = if c > 0.0 && da > constants::EM_DENOM_THRESH {
                (c * da.ln(), c)
            } else {
                (0.0, 0.0)
            };
            let (lb, rb) = if c > 0.0 && db > constants::EM_DENOM_THRESH {
                (c * db.ln(), c)
            } else {
                (0.0, 0.0)
            };
            (la, ra, lb, rb)
        };
        let zero = (0.0, 0.0, 0.0, 0.0);
        let (mut lla, mut ra, mut llb, mut rb) = if self.is_large() {
            let parts: Vec<Acc> = (0..self.blk.len() - 1)
                .into_par_iter()
                .map(|b| {
                    (self.blk[b] as usize..self.blk[b + 1] as usize)
                        .map(class_ll)
                        .fold(zero, add)
                })
                .collect();
            parts.into_iter().fold(zero, add)
        } else {
            (0..self.cls.len()).map(class_ll).fold(zero, add)
        };
        for ((&s, &ta), &tb) in base.iter().zip(a).zip(b) {
            if s > 0.0 {
                if ta > constants::EM_DENOM_THRESH {
                    lla += s * ta.ln();
                    ra += s;
                }
                if tb > constants::EM_DENOM_THRESH {
                    llb += s * tb.ln();
                    rb += s;
                }
            }
        }
        let (mut ta, mut tb) = (0.0_f64, 0.0_f64);
        for (&x, &y) in a.iter().zip(b) {
            ta += x;
            tb += y;
        }
        (lla - ra * ta.ln(), llb - rb * tb.ln())
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
        let nk = self.count_of(n);
        self.step_by(nk, &s.base, &s.t0, &mut s.t1, &mut s.acc, rel_change);
        let rd2 = self.step_by(nk, &s.base, &s.t1, &mut s.t2, &mut s.acc, rel_change);
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
        let rd3 = self.step_by(nk, &s.base, &s.prev, &mut s.curr, &mut s.acc, rel_change);
        let (l3, l2) = self.loglik_pair(n, &s.base, &s.curr, &s.t2);
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
                    let r = cm.step_by(
                        cm.count_of(n),
                        &s.base,
                        &s.prev,
                        &mut s.curr,
                        &mut s.acc,
                        rel_change,
                    );
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
            cm.step_by(
                cm.count_of(n),
                &s.base,
                &s.prev,
                &mut s.curr,
                &mut s.acc,
                rel_change,
            );
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::em_tests::{
        Classes, Read, RefEm, Rng, assert_bits_eq, assert_close, check_invariants, mixed_reads, on,
        ref_em,
    };

    /// Classes in CSR form: class `c` is entries `off[c]..off[c + 1]` of `t` and `w`.
    struct Csr {
        off: Vec<usize>,
        t: Vec<u32>,
        w: Vec<f64>,
    }

    fn csr(reads: &[Read]) -> Csr {
        let mut c = Csr {
            off: vec![0],
            t: vec![],
            w: vec![],
        };
        for r in reads {
            for &(t, p) in r {
                c.t.push(t);
                c.w.push(p as f64);
            }
            c.off.push(c.t.len());
        }
        c
    }

    fn plan_of(c: &Csr, ntx: usize, by_heaviest: bool) -> Plan<f64> {
        Plan::new(|k| c.t[k], |k| c.w[k], &c.off, ntx, by_heaviest)
    }

    /// Each component's transcripts and its classes (class id, then (global target, weight
    /// bits) entries), and the single-target classes.
    type View = (
        Vec<(Vec<u32>, Vec<(u32, Vec<(u32, u64)>)>)>,
        Vec<(u32, u32)>,
    );

    fn view(p: &Plan<f64>) -> View {
        let comps = p
            .comps
            .iter()
            .map(|cm| {
                let cls = (0..cm.cls.len())
                    .map(|k| {
                        let r = cm.off[k] as usize..cm.off[k + 1] as usize;
                        let e = r
                            .map(|i| (cm.tx[cm.lab[i] as usize], cm.w[i].to_bits()))
                            .collect();
                        (cm.cls[k], e)
                    })
                    .collect();
                (cm.tx.clone(), cls)
            })
            .collect();
        (comps, p.singles.clone())
    }

    /// Reads over 96 transcripts in overlapping windows of 8 (one in ten reaching into another
    /// window, so windows chain into larger components), with repeated targets, single-target
    /// reads on transcripts 100..104 only, and nothing on the rest of 96..120; weights drawn
    /// from a continuum, so a class's heaviest target is unique.
    fn chained_reads(seed: u64, n: usize) -> Vec<Read> {
        let mut rng = Rng::new(seed);
        let mut reads: Vec<Read> = (0..n)
            .map(|_| {
                let g = rng.below(12) * 8;
                let k = 1 + rng.below(3);
                let mut r: Read = (0..k)
                    .map(|_| {
                        (
                            (g + rng.below(8)) as u32,
                            0.05 + rng.below(1 << 20) as f32 / (1 << 20) as f32,
                        )
                    })
                    .collect();
                if rng.below(10) == 0 {
                    r.push((rng.below(96) as u32, 0.3));
                }
                r
            })
            .collect();
        reads.extend((100..104).flat_map(|t| std::iter::repeat_n(vec![(t, 0.7)], 3)));
        rng.shuffle(&mut reads);
        reads
    }

    /// Transcript labels after merging every class's targets (a plain fixed-point flood fill).
    fn flood(c: &Csr, ntx: usize) -> Vec<usize> {
        let mut lab: Vec<usize> = (0..ntx).collect();
        loop {
            let mut changed = false;
            for w in c.off.windows(2).filter(|w| w[1] - w[0] >= 2) {
                let m = c.t[w[0]..w[1]]
                    .iter()
                    .map(|&t| lab[t as usize])
                    .min()
                    .unwrap();
                for &t in &c.t[w[0]..w[1]] {
                    if lab[t as usize] != m {
                        lab[t as usize] = m;
                        changed = true;
                    }
                }
            }
            if !changed {
                return lab;
            }
        }
    }

    /// The layout partitions the multi-target classes by connected components of the
    /// class-transcript graph: every such class sits whole, in order, in the one component that
    /// holds its targets; single-target classes are kept apart; empty classes and transcripts no
    /// multi-target class touches are in no component. The layout is the same on any number of
    /// threads.
    #[test]
    fn plan_partitions_classes_by_connected_component() {
        let mut reads = chained_reads(11, 600);
        let ntx = 120;
        let mut c = csr(&reads);
        // Empty classes at a few positions.
        for at in [0, 17, 300, 600] {
            c.off.insert(at, c.off[at]);
            reads.insert(at, vec![]);
        }
        let lab = flood(&c, ntx);
        for by_heaviest in [false, true] {
            let plan = plan_of(&c, ntx, by_heaviest);
            let (comps, singles) = view(&plan);
            assert_eq!(comps.len(), plan.n_components());
            let mut comp_of = vec![usize::MAX; ntx];
            let mut placed = vec![0; reads.len()];
            for (k, (tx, cls)) in comps.iter().enumerate() {
                assert!(tx.windows(2).all(|w| w[0] < w[1]), "local ids ascend");
                for &t in tx {
                    assert_eq!(comp_of[t as usize], usize::MAX, "t{t} in two components");
                    comp_of[t as usize] = k;
                }
                for (cid, e) in cls {
                    placed[*cid as usize] += 1;
                    let want: Vec<(u32, u64)> = reads[*cid as usize]
                        .iter()
                        .map(|&(t, p)| (t, (p as f64).to_bits()))
                        .collect();
                    assert_eq!(e, &want, "class {cid}");
                    assert!(e.iter().all(|&(t, _)| tx.contains(&t)));
                }
            }
            for (cid, r) in reads.iter().enumerate() {
                assert_eq!(placed[cid], usize::from(r.len() >= 2), "class {cid}");
            }
            let want_singles: Vec<(u32, u32)> = reads
                .iter()
                .enumerate()
                .filter(|(_, r)| r.len() == 1)
                .map(|(cid, r)| (r[0].0, cid as u32))
                .collect();
            assert_eq!(singles, want_singles);
            let in_multi: Vec<bool> = (0..ntx)
                .map(|t| {
                    reads
                        .iter()
                        .any(|r| r.len() >= 2 && r.iter().any(|x| x.0 == t as u32))
                })
                .collect();
            for t in 0..ntx {
                assert_eq!(comp_of[t] != usize::MAX, in_multi[t], "t{t}");
                for u in 0..ntx {
                    if in_multi[t] && in_multi[u] {
                        assert_eq!(comp_of[t] == comp_of[u], lab[t] == lab[u], "t{t} t{u}");
                    }
                }
            }
            for threads in [1, 4, 16] {
                assert_eq!(
                    on(threads, || view(&plan_of(&c, ntx, by_heaviest))),
                    (comps.clone(), singles.clone())
                );
            }
        }
    }

    /// With `by_heaviest`, a component's classes are ordered by their heaviest target, then by
    /// class id; without it they keep class order.
    #[test]
    fn plan_orders_classes_by_heaviest_target() {
        let reads = chained_reads(12, 800);
        let c = csr(&reads);
        for by_heaviest in [false, true] {
            let (comps, _) = view(&plan_of(&c, 120, by_heaviest));
            for (tx, cls) in &comps {
                let key: Vec<(usize, u32)> = cls
                    .iter()
                    .map(|(cid, e)| {
                        let (t, _) = e
                            .iter()
                            .max_by(|a, b| f64::from_bits(a.1).total_cmp(&f64::from_bits(b.1)))
                            .unwrap();
                        let local = tx.binary_search(t).unwrap();
                        (if by_heaviest { local } else { 0 }, *cid)
                    })
                    .collect();
                assert!(
                    key.windows(2).all(|w| w[0] < w[1]),
                    "{by_heaviest}: {key:?}"
                );
            }
        }
    }

    /// `run` takes exactly the reference loop's rounds, wherever the stop falls relative to the
    /// tiles: with an infinite threshold it stops at the first round past `min_iter`, with a zero
    /// threshold at `max_iter`.
    #[test]
    fn run_takes_the_reference_rounds_across_tiles() {
        let reads = chained_reads(13, 500);
        let ntx = 120;
        let c = csr(&reads);
        let plan = plan_of(&c, ntx, true);
        let cl = Classes::of_reads(&reads, ntx);
        let init = cl.uniform();
        let mut cases: Vec<RefEm> = [0, 1, 2, 62, 63, 64, 65, 127, 128, 200]
            .map(|m| RefEm::new(1000, f64::INFINITY, m))
            .to_vec();
        cases.extend([0, 1, 2, 63, 64, 65, 129].map(|m| RefEm::new(m, 0.0, 1)));
        for o in cases {
            let (want, _) = ref_em(&cl, &init, o);
            let got = run(&plan, Counts::Unit, &init, o.max_iter, o.thresh, o.min_iter);
            let what = format!(
                "max_iter {} thresh {} min_iter {}",
                o.max_iter, o.thresh, o.min_iter
            );
            assert_close(&got, &want, 1e-11, &what);
        }
    }

    /// A class with no reads changes nothing: the run with some counts at zero is the run
    /// without those classes, and transcripts only they touch get 0.
    #[test]
    fn zero_count_classes_are_ignored() {
        let reads = chained_reads(14, 400);
        let ntx = 120;
        let mut rng = Rng::new(15);
        let counts: Vec<f64> = (0..reads.len())
            .map(|_| [0.0, 1.0, 2.5, 7.0][rng.below(4)])
            .collect();
        let kept: Vec<Read> = reads
            .iter()
            .zip(&counts)
            .filter(|(_, n)| **n > 0.0)
            .map(|(r, _)| r.clone())
            .collect();
        let kept_n: Vec<f64> = counts.iter().copied().filter(|&n| n > 0.0).collect();
        let cl = Classes::with_counts(&kept, &kept_n, ntx);
        let init = cl.uniform();
        let (want, _) = ref_em(&cl, &init, RefEm::new(2000, 1e-10, 1));
        let plan = plan_of(&csr(&reads), ntx, false);
        let got = run(&plan, Counts::Of(&counts), &init, 2000, 1e-10, 1);
        assert_close(&got, &want, 1e-7, "zero counts");
        let total: f64 = kept_n.iter().sum();
        assert!((got.iter().sum::<f64>() - total).abs() < 1e-9 * total);
    }

    /// With `plain_first` at or past `max_iter`, `run_squarem` is the plain loop with the
    /// truncation EM's statistic (relative change in either direction), from iteration
    /// `niter0`, and reports that loop's iteration count and no SQUAREM cycles.
    #[test]
    fn run_squarem_without_cycles_is_the_plain_loop() {
        let reads = mixed_reads(16);
        let ntx = crate::em_tests::MIXED_NTX;
        let plan = plan_of(&csr(&reads), ntx, true);
        let cl = Classes::of_reads(&reads, ntx);
        let init = cl.uniform();
        let mut cases = vec![];
        for niter0 in [0, 1, 2, 5, 70] {
            cases.push((niter0, 1000, f64::INFINITY, 1e-11));
            cases.push((niter0, 3000, 1e-9, 1e-7));
        }
        for max_iter in [0, 1, 2, 64, 65, 130] {
            cases.push((0, max_iter, 0.0, 1e-11));
        }
        for (niter0, max_iter, thresh, rtol) in cases {
            let o = RefEm {
                max_iter,
                thresh,
                min_iter: 1,
                niter0,
                abs: true,
            };
            let (want, wn) = ref_em(&cl, &init, o);
            let (got, rep) = run_squarem(
                &plan,
                Counts::Unit,
                &init,
                niter0,
                max_iter,
                thresh,
                max_iter,
            );
            let what = format!("niter0 {niter0} max_iter {max_iter} thresh {thresh}");
            assert_close(&got, &want, rtol, &what);
            if rtol < 1e-9 {
                assert_eq!(rep.niter, wn, "{what}");
            }
            assert_eq!((rep.acc, rep.rej, rep.plain), (0, 0, 0), "{what}");
        }
    }

    /// Nested, nearly indistinguishable transcripts (plain EM's slow case), next to the mixed
    /// fixture's components.
    fn slow_reads() -> (Vec<Read>, usize) {
        let mut reads = mixed_reads(17);
        let b = crate::em_tests::MIXED_NTX as u32;
        for i in 0..60 {
            reads.push(vec![(b, 1.0), (b + 1, 1.0), (b + 2, 0.9)]);
            if i % 3 == 0 {
                reads.push(vec![(b, 1.0)]);
            }
            if i % 5 == 0 {
                reads.push(vec![(b + 3, 1.0), (b + 2, 0.5)]);
            }
        }
        (reads, b as usize + 4)
    }

    /// SQUAREM per component reaches the plain EM's fixed point, extrapolating where plain EM
    /// is slow, from the start or after plain steps, and keeps the invariants.
    #[test]
    fn run_squarem_reaches_the_fixed_point() {
        let (reads, ntx) = slow_reads();
        let plan = plan_of(&csr(&reads), ntx, true);
        let cl = Classes::of_reads(&reads, ntx);
        let init = cl.uniform();
        let mut o = RefEm::new(200_000, 1e-14, 1);
        o.abs = true;
        let (fixed, _) = ref_em(&cl, &init, o);
        for plain_first in [0, 10] {
            let (got, rep) = run_squarem(&plan, Counts::Unit, &init, 0, 20_000, 1e-10, plain_first);
            assert_close(&got, &fixed, 1e-5, &format!("plain_first {plain_first}"));
            check_invariants(&got, &reads, None, ntx, "run_squarem");
            assert!(rep.acc > 0, "never extrapolated");
            assert!(rep.niter < 20_000, "did not converge");
        }
        // Multiplicities are repeated reads.
        let mult = Rng::new(18).resample(reads.len());
        let dup: Vec<Read> = reads
            .iter()
            .zip(&mult)
            .flat_map(|(r, &m)| std::iter::repeat_n(r.clone(), m as usize))
            .collect();
        let dc = Classes::of_reads(&dup, ntx);
        let (fixed, _) = ref_em(&dc, &init, o);
        let (got, _) = run_squarem(&plan, Counts::OfU32(&mult), &init, 0, 20_000, 1e-10, 0);
        assert_close(&got, &fixed, 1e-5, "multiplicities");
        check_invariants(&got, &reads, Some(&mult), ntx, "run_squarem OfU32");
    }

    /// The single-target reads of component `k`'s transcripts, by local id.
    fn base_of(plan: &Plan<f64>, k: usize, n: Counts) -> Vec<f64> {
        let cm = &plan.comps[k];
        let mut base = vec![0.0; cm.tx.len()];
        for &(t, c) in &plan.singles {
            if let Ok(l) = cm.tx.binary_search(&t) {
                base[l] += n.get(c);
            }
        }
        base
    }

    /// Checks component `k`'s log-likelihood against the dense formula over its classes and
    /// single-target reads (up to the constant the latter's weights add), its invariance to scaling, and that an EM round does not lower it.
    fn check_loglik(plan: &Plan<f64>, k: usize, reads: &[Read], counts: &[f64], theta: &[f64]) {
        let cm = &plan.comps[k];
        let n = Counts::Of(counts);
        let base = base_of(plan, k, n);
        let local: Vec<f64> = cm.tx.iter().map(|&t| theta[t as usize]).collect();
        let mut masked = vec![0.0; theta.len()];
        for &t in &cm.tx {
            masked[t as usize] = theta[t as usize];
        }
        let touching: Vec<usize> = (0..reads.len())
            .filter(|&r| reads[r].iter().any(|x| cm.tx.binary_search(&x.0).is_ok()))
            .collect();
        let sub = Classes::with_counts(
            &touching
                .iter()
                .map(|&r| reads[r].clone())
                .collect::<Vec<_>>(),
            &touching.iter().map(|&r| counts[r]).collect::<Vec<_>>(),
            theta.len(),
        );
        let ll = cm.loglik(n, &base, &local);
        // The single-target reads enter as `n ln theta`: their weights are a constant.
        let consts: f64 = touching
            .iter()
            .filter(|&&r| reads[r].len() == 1 && counts[r] > 0.0)
            .map(|&r| counts[r] * (reads[r][0].1 as f64).ln())
            .sum();
        let want = sub.loglik(&masked) - consts;
        assert!(
            (ll - want).abs() <= 1e-10 * want.abs().max(1.0),
            "{ll} vs {want}"
        );
        let scaled: Vec<f64> = local.iter().map(|x| 2.5 * x).collect();
        let ls = cm.loglik(n, &base, &scaled);
        assert!((ls - ll).abs() <= 1e-10 * ll.abs().max(1.0), "{ls} vs {ll}");
        let mut next = vec![0.0; local.len()];
        let mut acc = vec![0.0; local.len() * cm.blk.len().saturating_sub(1).max(1)];
        cm.step(n, &base, &local, &mut next, &mut acc);
        let ln = cm.loglik(n, &base, &next);
        assert!(
            ln >= ll - 1e-9 * ll.abs(),
            "EM round lowered the likelihood: {ll} -> {ln}"
        );
    }

    #[test]
    fn loglik_matches_the_dense_formula() {
        let (reads, ntx) = slow_reads();
        let mut rng = Rng::new(19);
        let counts: Vec<f64> = (0..reads.len()).map(|_| rng.below(4) as f64).collect();
        let plan = plan_of(&csr(&reads), ntx, true);
        let theta: Vec<f64> = (0..ntx).map(|_| 0.5 + rng.below(100) as f64).collect();
        for k in 0..plan.n_components() {
            check_loglik(&plan, k, &reads, &counts, &theta);
        }
    }

    #[test]
    fn convergence_statistics() {
        let tiny = constants::MIN_READ_THRESH / 2.0;
        assert_eq!(rel_change(tiny, 5.0), 0.0);
        assert_eq!(rel_increase(tiny, 5.0), 0.0);
        assert_eq!(rel_change(2.0, 1.0), 0.5);
        assert_eq!(rel_increase(2.0, 1.0), -0.5);
        assert_eq!(rel_change(1.0, 3.0), 2.0);
        assert_eq!(rel_increase(1.0, 3.0), 2.0);
        assert_eq!(Counts::Unit.get(7), 1.0);
        assert_eq!(Counts::Of(&[0.5, 2.0]).get(1), 2.0);
        assert_eq!(Counts::OfU32(&[3, 4]).get(0), 3.0);
    }

    /// A component of at least [`LARGE`] incidences is split into at most [`BLOCKS`] blocks of
    /// whole classes and stepped in parallel; the layout, `run` and `run_squarem` give the same
    /// bits on 1, 4 and 16 threads, `run` takes the reference rounds, and the blocked
    /// log-likelihood matches the dense formula.
    #[test]
    fn large_component_is_blocked_and_independent_of_threads() {
        let mut rng = Rng::new(20);
        // 40,000 reads over 200 transcripts: one component of about 90,000 incidences, next to
        // the mixed fixture's small ones (moved to transcripts 200..264).
        let mut reads = crate::em_tests::grouped_reads(&mut rng, 40_000, 1, 200, false);
        reads.extend(
            mixed_reads(21)
                .into_iter()
                .map(|r| r.into_iter().map(|(t, p)| (t + 200, p)).collect::<Read>()),
        );
        rng.shuffle(&mut reads);
        let ntx = 200 + crate::em_tests::MIXED_NTX;
        let c = csr(&reads);
        let plan = on(1, || plan_of(&c, ntx, true));
        assert!(plan.has_large());
        for cm in plan.comps.iter().filter(|cm| cm.is_large()) {
            let b = &cm.blk;
            assert!(
                b.len() >= 3 && b.len() - 1 <= BLOCKS,
                "{} blocks",
                b.len() - 1
            );
            assert_eq!((b[0], *b.last().unwrap()), (0, cm.cls.len() as u32));
            assert!(b.windows(2).all(|w| w[0] < w[1]));
        }
        let v = view(&plan);
        for threads in [4, 16] {
            assert!(
                on(threads, || view(&plan_of(&c, ntx, true))) == v,
                "layout on {threads} threads"
            );
        }
        let cl = Classes::of_reads(&reads, ntx);
        let init = cl.uniform();
        let runs: Vec<Vec<f64>> = [1, 4, 16]
            .map(|t| on(t, || run(&plan, Counts::Unit, &init, 40, 0.0, 1)))
            .to_vec();
        let (want, _) = ref_em(&cl, &init, RefEm::new(40, 0.0, 1));
        assert_close(&runs[0], &want, 1e-10, "large run");
        let sq: Vec<Vec<f64>> = [1, 4, 16]
            .map(|t| {
                on(t, || {
                    run_squarem(&plan, Counts::Unit, &init, 0, 30, 0.0, 0).0
                })
            })
            .to_vec();
        for i in 1..3 {
            assert_bits_eq(&runs[i], &runs[0], "run");
            assert_bits_eq(&sq[i], &sq[0], "run_squarem");
        }
        check_invariants(&sq[0], &reads, None, ntx, "large run_squarem");
        let k = (0..plan.n_components())
            .find(|&k| plan.comps[k].is_large())
            .unwrap();
        let theta: Vec<f64> = (0..ntx).map(|_| 0.5 + rng.below(100) as f64).collect();
        check_loglik(&plan, k, &reads, &vec![1.0; reads.len()], &theta);
    }
}
