//! Equivalence-class EM for when no per-read model term is active.
//!
//! Without the coverage model and the KDE length model, a read enters the EM only through the
//! multiset of its (target, alignment-score probability) pairs. Reads with the same multiset
//! contribute identical terms, so they are collapsed into one class weighted by its read count.
//! Probabilities are keyed by their exact bits, so the class EM computes the same estimates as
//! the per-read EM (up to floating-point summation order) while touching each distinct class once
//! per iteration. Each iteration uses the CSR layout of `crate::em_read`: one pass over classes
//! computes each class's scale, and one pass over transcripts sums the scales through a
//! transpose, so no update is shared between threads.
//!
//! Without the presence model (which couples all transcripts through its base rate), the same
//! rounds run on the component layout of [`crate::em_components`] instead.

use hashbrown::HashTable;
use num_format::{Locale, ToFormattedString};
use rand::rng as trng;
use rayon::prelude::*;
use rustc_hash::FxBuildHasher;
use std::hash::BuildHasher;
use tracing::{info, span, trace};

use crate::bootstrap;
use crate::em_components::{self, Counts, Plan};
use crate::presence::{Presence, PresenceOpts};
use crate::util::constants;
use crate::util::oarfish_types::{
    InMemoryAlignmentStore, InMemoryAlignmentStoreIter, TranscriptInfo, release_freed_memory,
};
use crate::util::run_warnings::EmConvergence;

/// Reads collapsed by their (target, probability) multisets.
pub struct EqClasses {
    /// Targets of class `c`: `targets[offsets[c]..offsets[c + 1]]`, sorted.
    targets: Vec<u32>,
    /// Alignment-score probabilities, parallel to `targets`.
    probs: Vec<f32>,
    offsets: Vec<usize>,
    /// Reads per class.
    counts: Vec<f64>,
    /// Class of each read in the store (for bootstrap resampling).
    read_class: Vec<u32>,
    /// Transpose: transcript `t`'s (class, probability) entries are
    /// `t_class[t_off[t]..t_off[t + 1]]` and the same range of `t_prob`.
    t_off: Vec<usize>,
    t_class: Vec<u32>,
    t_prob: Vec<f32>,
}

/// Reads per chunk below which [`EqClasses::from_store`] adds no further chunk.
const MIN_CHUNK_READS: usize = 16_384;

/// One chunk's classes, numbered in order of first occurrence within the chunk.
struct Local {
    /// Sorted (target, probability bits) key of local class `c`: `keys[offsets[c]..offsets[c + 1]]`.
    keys: Vec<(u32, u32)>,
    offsets: Vec<usize>,
    /// Reads per local class.
    counts: Vec<u32>,
    /// Local class of each of the chunk's reads.
    read_class: Vec<u32>,
}

impl Local {
    fn build<R, I>(reads: R, nreads: usize) -> Self
    where
        R: Iterator<Item = I>,
        I: Iterator<Item = (u32, f32)>,
    {
        // Local class ids, hashed and compared through their keys in `local.keys`, so a class's
        // key is stored once, flat, rather than also as an owned key per map entry.
        let mut ids: HashTable<u32> = HashTable::new();
        let mut local = Local {
            keys: Vec::new(),
            offsets: vec![0],
            counts: Vec::new(),
            read_class: Vec::with_capacity(nreads),
        };
        let mut key: Vec<(u32, u32)> = Vec::new();
        for read in reads {
            key.clear();
            key.extend(read.map(|(t, p)| (t, p.to_bits())));
            key.sort_unstable();
            let h = FxBuildHasher.hash_one(key.as_slice());
            let c = match ids.find(h, |&c| local.key(c as usize) == key.as_slice()) {
                Some(&c) => c,
                None => {
                    let c = local.counts.len() as u32;
                    local.keys.extend_from_slice(&key);
                    local.offsets.push(local.keys.len());
                    local.counts.push(0);
                    ids.insert_unique(h, c, |&c| FxBuildHasher.hash_one(local.key(c as usize)));
                    c
                }
            };
            local.counts[c as usize] += 1;
            local.read_class.push(c);
        }
        // The chunks are keyed concurrently, each on its own worker; return this worker's
        // freed memory (the map, the keys' regrowth) before the merge allocates its own.
        drop(ids);
        release_freed_memory();
        local
    }

    fn len(&self) -> usize {
        self.counts.len()
    }

    fn key(&self, c: usize) -> &[(u32, u32)] {
        &self.keys[self.offsets[c]..self.offsets[c + 1]]
    }
}

impl EqClasses {
    /// Classes of the store's reads. The reads are keyed in chunks in parallel on the current
    /// rayon pool, then merged in chunk order, so the classes, their order (first occurrence
    /// in the store), counts and `read_class` do not depend on the number of threads.
    pub fn from_store(store: &InMemoryAlignmentStore, ntx: usize) -> Self {
        // Memory the parse freed is still resident until the allocator purges it.
        release_freed_memory();
        let n = store.len();
        let nchunks = rayon::current_num_threads()
            .min(n.div_ceil(MIN_CHUNK_READS))
            .max(1);
        let locals: Vec<Local> = (0..nchunks)
            .into_par_iter()
            .map(|i| {
                let (lo, hi) = (i * n / nchunks, (i + 1) * n / nchunks);
                let reads = InMemoryAlignmentStoreIter { store, idx: lo }
                    .take(hi - lo)
                    .map(|(alns, probs, _)| alns.iter().zip(probs).map(|(a, &p)| (a.ref_id, p)));
                Local::build(reads, hi - lo)
            })
            .collect();
        Self::merge(locals, ntx)
    }

    /// Classes from each read's (target, probability) pairs, over `ntx` transcripts.
    #[cfg(test)]
    fn from_reads<R, I>(reads: R, nreads: usize, ntx: usize) -> Self
    where
        R: Iterator<Item = I>,
        I: Iterator<Item = (u32, f32)>,
    {
        Self::merge(vec![Local::build(reads, nreads)], ntx)
    }

    /// The classes of consecutive chunks of reads: a class's global id is its first
    /// occurrence over the chunks in order, so the result is that of one chunk of all reads.
    fn merge(mut locals: Vec<Local>, ntx: usize) -> Self {
        let mut eqc = if locals.len() == 1 {
            let l = locals.pop().unwrap();
            EqClasses {
                targets: l.keys.iter().map(|&(t, _)| t).collect(),
                probs: l.keys.iter().map(|&(_, p)| f32::from_bits(p)).collect(),
                offsets: l.offsets,
                counts: l.counts.iter().map(|&k| k as f64).collect(),
                read_class: l.read_class,
                t_off: Vec::new(),
                t_class: Vec::new(),
                t_prob: Vec::new(),
            }
        } else {
            Self::merge_sharded(&locals)
        };
        drop(locals);
        release_freed_memory();
        eqc.transpose(ntx);
        eqc
    }

    /// [`EqClasses::merge`] of several chunks, in parallel: the local classes are split into
    /// shards by key hash, each shard numbers its distinct keys in (chunk, local class) order,
    /// and a class's global id is the rank of its first (chunk, local class) among all firsts.
    fn merge_sharded(locals: &[Local]) -> Self {
        let nshards = 4 * rayon::current_num_threads();
        let shard_of = |key: &[(u32, u32)]| {
            let h = FxBuildHasher.hash_one(key);
            ((h.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 32) as usize) % nshards
        };
        // The local classes of each chunk in each shard, in order.
        let buckets: Vec<Vec<Vec<u32>>> = locals
            .par_iter()
            .map(|l| {
                let mut b = vec![Vec::new(); nshards];
                for c in 0..l.len() {
                    b[shard_of(l.key(c))].push(c as u32);
                }
                b
            })
            .collect();
        // Each shard's distinct keys: the (chunk, local class) of the first, their reads, and
        // the shard class of each bucket entry, per chunk.
        struct Shard {
            first: Vec<(u32, u32)>,
            counts: Vec<u64>,
            ids: Vec<Vec<u32>>,
        }
        let shards: Vec<Shard> = (0..nshards)
            .into_par_iter()
            .map(|s| {
                // Shard class ids, keyed through their first (chunk, local class); sized for
                // every bucket entry, so the table never grows.
                let entries: usize = buckets.iter().map(|b| b[s].len()).sum();
                let mut map: HashTable<u32> = HashTable::with_capacity(entries);
                let key_of = |first: &[(u32, u32)], id: u32| {
                    let (k, c) = first[id as usize];
                    locals[k as usize].key(c as usize)
                };
                let mut sh = Shard {
                    first: Vec::new(),
                    counts: Vec::new(),
                    ids: Vec::with_capacity(locals.len()),
                };
                for (k, l) in locals.iter().enumerate() {
                    let bucket = &buckets[k][s];
                    let mut ids = Vec::with_capacity(bucket.len());
                    for &c in bucket {
                        let key = l.key(c as usize);
                        let h = FxBuildHasher.hash_one(key);
                        let id = match map.find(h, |&id| key_of(&sh.first, id) == key) {
                            Some(&id) => id,
                            None => {
                                let id = sh.first.len() as u32;
                                sh.first.push((k as u32, c));
                                sh.counts.push(0);
                                map.insert_unique(h, id, |&id| {
                                    FxBuildHasher.hash_one(key_of(&sh.first, id))
                                });
                                id
                            }
                        };
                        sh.counts[id as usize] += u64::from(l.counts[c as usize]);
                        ids.push(id);
                    }
                    sh.ids.push(ids);
                }
                sh
            })
            .collect();
        // Per chunk: each local class's (shard, shard class), and the local classes that
        // are first occurrences, in order.
        type ShardRefs = Vec<(u32, u32)>;
        let (refs, firsts): (Vec<ShardRefs>, Vec<Vec<u32>>) = locals
            .par_iter()
            .enumerate()
            .map(|(k, l)| {
                let mut refs = vec![(0_u32, 0_u32); l.len()];
                let mut is_first = vec![false; l.len()];
                for (s, sh) in shards.iter().enumerate() {
                    for (&c, &id) in buckets[k][s].iter().zip(&sh.ids[k]) {
                        refs[c as usize] = (s as u32, id);
                        is_first[c as usize] = sh.first[id as usize] == (k as u32, c);
                    }
                }
                let firsts = (0..l.len() as u32)
                    .filter(|&c| is_first[c as usize])
                    .collect();
                (refs, firsts)
            })
            .unzip();
        drop(buckets);
        // Chunk `k`'s first occurrences are classes `base[k]..base[k + 1]`.
        let mut base = vec![0_usize; locals.len() + 1];
        for k in 0..locals.len() {
            base[k + 1] = base[k] + firsts[k].len();
        }
        let nclasses = base[locals.len()];
        let mut offsets = Vec::with_capacity(nclasses + 1);
        offsets.push(0_usize);
        for (l, f) in locals.iter().zip(&firsts) {
            for &c in f {
                let c = c as usize;
                offsets.push(offsets.last().unwrap() + (l.offsets[c + 1] - l.offsets[c]));
            }
        }
        // Global id of each shard class, from the chunk that holds its first occurrence.
        let global: Vec<Vec<u32>> = shards
            .par_iter()
            .map(|sh| {
                let mut rank: Vec<u32> = Vec::with_capacity(sh.first.len());
                for &(k, c) in &sh.first {
                    let f = &firsts[k as usize];
                    let i = f.binary_search(&c).expect("a first occurrence is listed");
                    rank.push((base[k as usize] + i) as u32);
                }
                rank
            })
            .collect();
        let nentries = *offsets.last().unwrap();
        let mut targets = vec![0_u32; nentries];
        let mut probs = vec![0_f32; nentries];
        let mut counts = vec![0_f64; nclasses];
        {
            let mut parts = Vec::with_capacity(locals.len());
            let (mut t_rest, mut p_rest, mut c_rest) = (
                targets.as_mut_slice(),
                probs.as_mut_slice(),
                counts.as_mut_slice(),
            );
            for k in 0..locals.len() {
                let len = offsets[base[k + 1]] - offsets[base[k]];
                let (t, tr) = t_rest.split_at_mut(len);
                let (p, pr) = p_rest.split_at_mut(len);
                let (c, cr) = c_rest.split_at_mut(base[k + 1] - base[k]);
                parts.push((t, p, c));
                (t_rest, p_rest, c_rest) = (tr, pr, cr);
            }
            parts
                .into_par_iter()
                .enumerate()
                .for_each(|(k, (ts, ps, cs))| {
                    let l = &locals[k];
                    let mut e = 0;
                    for (&c, cnt) in firsts[k].iter().zip(cs.iter_mut()) {
                        let key = l.key(c as usize);
                        for &(t, p) in key {
                            ts[e] = t;
                            ps[e] = f32::from_bits(p);
                            e += 1;
                        }
                        let (s, id) = refs[k][c as usize];
                        *cnt = shards[s as usize].counts[id as usize] as f64;
                    }
                });
        }
        let n: usize = locals.iter().map(|l| l.read_class.len()).sum();
        let mut read_class = vec![0_u32; n];
        let mut parts: Vec<&mut [u32]> = Vec::with_capacity(locals.len());
        let mut rest = read_class.as_mut_slice();
        for l in locals {
            let (a, b) = rest.split_at_mut(l.read_class.len());
            parts.push(a);
            rest = b;
        }
        parts
            .into_par_iter()
            .zip(locals.par_iter().zip(refs.par_iter()))
            .for_each(|(out, (l, refs))| {
                for (o, &c) in out.iter_mut().zip(&l.read_class) {
                    let (s, id) = refs[c as usize];
                    *o = global[s as usize][id as usize];
                }
            });
        EqClasses {
            targets,
            probs,
            offsets,
            counts,
            read_class,
            t_off: Vec::new(),
            t_class: Vec::new(),
            t_prob: Vec::new(),
        }
    }

    /// Fill `t_off`, `t_class` and `t_prob` by a counting sort of the entries on target.
    fn transpose(&mut self, ntx: usize) {
        let mut off = vec![0_usize; ntx + 1];
        for &t in &self.targets {
            off[t as usize + 1] += 1;
        }
        for t in 0..ntx {
            off[t + 1] += off[t];
        }
        let mut next = off[..ntx].to_vec();
        let mut class = vec![0_u32; self.targets.len()];
        let mut prob = vec![0_f32; self.targets.len()];
        for c in 0..self.len() {
            let (ts, ps) = self.class(c);
            for (&t, &p) in ts.iter().zip(ps) {
                let k = &mut next[t as usize];
                class[*k] = c as u32;
                prob[*k] = p;
                *k += 1;
            }
        }
        (self.t_off, self.t_class, self.t_prob) = (off, class, prob);
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// The component layout of the classes (see [`crate::em_components`]).
    fn plan(&self, ntx: usize) -> Plan<f32> {
        Plan::new(
            |k| self.targets[k],
            |k| self.probs[k],
            &self.offsets,
            ntx,
            false,
        )
    }

    fn class(&self, c: usize) -> (&[u32], &[f32]) {
        let r = self.offsets[c]..self.offsets[c + 1];
        (&self.targets[r.clone()], &self.probs[r])
    }
}

/// One EM round: distribute each class's reads over its targets in proportion to
/// `prev[t] * prob`, writing `curr`. `scale[c]` receives class `c`'s reads over its normalizer.
fn m_step(eqc: &EqClasses, counts: &[f64], prev: &[f64], scale: &mut [f64], curr: &mut [f64]) {
    scale.par_iter_mut().enumerate().for_each(|(c, s)| {
        let n = counts[c];
        *s = 0.0;
        if n == 0.0 {
            return;
        }
        let (ts, ps) = eqc.class(c);
        let denom: f64 = ts
            .iter()
            .zip(ps)
            .map(|(&t, &p)| prev[t as usize] * p as f64)
            .sum();
        if denom > constants::EM_DENOM_THRESH {
            *s = n / denom;
        }
    });
    curr.par_iter_mut()
        .with_min_len(1024)
        .enumerate()
        .for_each(|(t, x)| {
            let r = eqc.t_off[t]..eqc.t_off[t + 1];
            let sum: f64 = eqc.t_class[r.clone()]
                .iter()
                .zip(&eqc.t_prob[r])
                .map(|(&c, &p)| p as f64 * scale[c as usize])
                .sum();
            *x = prev[t] * sum;
        });
}

/// The initial abundances of [`run`]: `init_abundances`, or the reads spread evenly.
fn initial(counts: &[f64], ntx: usize, init_abundances: Option<&Vec<f64>>) -> Vec<f64> {
    match init_abundances {
        Some(init) => init.clone(),
        None => vec![counts.iter().sum::<f64>() / ntx as f64; ntx],
    }
}

/// The EM of [`crate::em::em_par`] over classes with read `counts` (the class counts, or a
/// bootstrap resample of them): same initialization, convergence rule and final pruning round.
/// Also returns how the loop stopped.
#[allow(clippy::too_many_arguments)]
fn run(
    eqc: &EqClasses,
    counts: &[f64],
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    presence: Option<PresenceOpts>,
    do_log: bool,
) -> (Vec<f64>, EmConvergence) {
    let mut prev = initial(counts, tinfo.len(), init_abundances);
    let mut curr = vec![0.0_f64; tinfo.len()];
    let mut scale = vec![0.0_f64; eqc.len()];
    // With the presence model the EM sees `q * prev` (see `crate::presence`).
    let mut pres = presence.map(|o| Presence::new(tinfo.len(), o));
    let mut eff = if pres.is_some() {
        vec![0.0_f64; tinfo.len()]
    } else {
        Vec::new()
    };

    let mut niter = 0_u32;
    let (mut converged, mut last_rel_diff) = (false, f64::INFINITY);
    while niter < max_iter {
        match &mut pres {
            Some(p) => {
                if p.due(niter) {
                    let reads = crate::presence::ReadAlns {
                        off: &eqc.offsets,
                        target: &|k| eqc.targets[k],
                    };
                    p.refresh(&prev, reads, |k| eqc.probs[k] as f64, |c| counts[c]);
                }
                p.effective(&prev, &mut eff);
                m_step(eqc, counts, &eff, &mut scale, &mut curr);
            }
            None => m_step(eqc, counts, &prev, &mut scale, &mut curr),
        }
        let mut rel_diff = 0.0_f64;
        for (&pc, &c) in prev.iter().zip(&curr) {
            if pc > constants::MIN_READ_THRESH {
                rel_diff = rel_diff.max((c - pc) / pc);
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        last_rel_diff = rel_diff;
        if rel_diff < convergence_thresh && niter > 1 {
            converged = true;
            break;
        }
        niter += 1;
        if do_log && niter.is_multiple_of(10) {
            if niter.is_multiple_of(100) {
                info!(
                    "iteration {}; rel diff {}",
                    niter.to_formatted_string(&Locale::en),
                    rel_diff
                );
            } else {
                trace!(
                    "iteration {}; rel diff {}",
                    niter.to_formatted_string(&Locale::en),
                    rel_diff
                );
            }
        }
    }

    // Zero very small abundances, then one more round, as the per-read EM does.
    for x in &mut prev {
        if *x < constants::MIN_READ_THRESH {
            *x = 0.0;
        }
    }
    match &pres {
        Some(p) => {
            if do_log {
                info!(
                    "presence model: {} transcripts called absent (q < 0.5){}",
                    p.n_absent().to_formatted_string(&Locale::en),
                    p.base_rate()
                        .map_or(String::new(), |pi| format!("; base rate {pi:.3}"))
                );
            }
            p.effective(&prev, &mut eff);
            m_step(eqc, counts, &eff, &mut scale, &mut curr);
        }
        None => m_step(eqc, counts, &prev, &mut scale, &mut curr),
    }
    let conv = EmConvergence {
        niter,
        max_iter,
        converged,
        rel_diff: last_rel_diff,
        thresh: convergence_thresh,
    };
    (curr, conv)
}

/// Abundance estimates (expected reads per target), and how the EM loop stopped.
pub fn em(
    eqc: &EqClasses,
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    presence: Option<PresenceOpts>,
    nthreads: usize,
) -> (Vec<f64>, EmConvergence) {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        if presence.is_none() && em_components::enabled() {
            let t0 = std::time::Instant::now();
            let plan = eqc.plan(tinfo.len());
            info!(
                "{} connected components, laid out in {:.2}s",
                plan.n_components().to_formatted_string(&Locale::en),
                t0.elapsed().as_secs_f64()
            );
            let init = initial(&eqc.counts, tinfo.len(), init_abundances);
            return em_components::run_counted(
                &plan,
                Counts::Of(&eqc.counts),
                &init,
                max_iter,
                convergence_thresh,
                1,
            );
        }
        run(
            eqc,
            &eqc.counts,
            tinfo,
            max_iter,
            convergence_thresh,
            init_abundances,
            presence,
            true,
        )
    })
}

/// Bootstrap replicates: each resamples the reads with replacement (as the per-read bootstrap
/// does) and runs the EM on the resulting class counts.
#[allow(clippy::too_many_arguments)]
pub fn bootstrap(
    eqc: &EqClasses,
    tinfo: &[TranscriptInfo],
    max_iter: u32,
    convergence_thresh: f64,
    init_abundances: Option<&Vec<f64>>,
    presence: Option<PresenceOpts>,
    num_boot: u32,
    nthreads: usize,
) -> Vec<Vec<f64>> {
    let span = span!(tracing::Level::INFO, "bootstrap");
    let _guard = span.enter();
    info!("will collection {num_boot} bootstraps");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    pool.install(|| {
        // One component layout serves every replicate; only the class counts differ.
        let plan = (presence.is_none() && em_components::enabled()).then(|| eqc.plan(tinfo.len()));
        (0..num_boot)
            .into_par_iter()
            .map(|i| {
                info!("evaluating bootstrap replicate {}", i);
                let mut rng = trng();
                let mut counts = vec![0.0_f64; eqc.len()];
                for r in bootstrap::get_sample_inds(eqc.read_class.len(), &mut rng) {
                    counts[eqc.read_class[r] as usize] += 1.0;
                }
                if let Some(plan) = &plan {
                    let init = initial(&counts, tinfo.len(), init_abundances);
                    return em_components::run(
                        plan,
                        Counts::Of(&counts),
                        &init,
                        max_iter,
                        convergence_thresh,
                        1,
                    );
                }
                run(
                    eqc,
                    &counts,
                    tinfo,
                    max_iter,
                    convergence_thresh,
                    init_abundances,
                    presence,
                    false,
                )
                .0
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reads over `ntx` transcripts in groups of `group` (each read's targets drawn from one
    /// group, so the groups are separate components), with some single-target reads, from a
    /// fixed pseudo-random sequence. The last transcript is never a target.
    fn grouped_reads(nreads: usize, ntx: usize, group: usize) -> Vec<Vec<(u32, f32)>> {
        let mut x: u64 = 0x9e3779b97f4a7c15;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let ngroups = (ntx - 1) / group;
        (0..nreads)
            .map(|_| {
                let g = (next() as usize) % ngroups;
                let k = 1 + (next() as usize) % 3;
                let mut ts: Vec<u32> = (0..k)
                    .map(|_| (g * group + (next() as usize) % group) as u32)
                    .collect();
                ts.sort_unstable();
                ts.dedup();
                ts.into_iter()
                    .map(|t| (t, [1.0_f32, 0.5, 0.25][(next() as usize) % 3]))
                    .collect()
            })
            .collect()
    }

    fn tinfo(ntx: usize) -> Vec<TranscriptInfo> {
        (0..ntx)
            .map(|_| TranscriptInfo::with_len(std::num::NonZeroUsize::new(1000).unwrap()))
            .collect()
    }

    /// The component layout runs the class loop's rounds: the same estimates, up to summation
    /// order, whether the global criterion stops it inside a tile or the iteration cap does.
    #[test]
    fn component_em_matches_class_em() {
        let (ntx, reads) = (301, grouped_reads(5000, 301, 6));
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), ntx);
        let ti = tinfo(ntx);
        let plan = eqc.plan(ntx);
        assert!(plan.n_components() > 10);
        for (max_iter, thresh) in [(1000, 1e-3), (1000, 1e-9), (150, 1e-12)] {
            let want = run(&eqc, &eqc.counts, &ti, max_iter, thresh, None, None, false).0;
            let init = initial(&eqc.counts, ntx, None);
            let got =
                em_components::run(&plan, Counts::Of(&eqc.counts), &init, max_iter, thresh, 1);
            for (t, (g, w)) in got.iter().zip(&want).enumerate() {
                assert!(
                    (g - w).abs() <= 1e-9 * w.max(1.0),
                    "{max_iter} {thresh} t{t}: {g} vs {w}"
                );
            }
            assert_eq!(got[ntx - 1], 0.0);
        }
    }

    /// The component EM, a large component's parallel rounds included, gives the same bits on
    /// any number of threads.
    #[test]
    fn component_em_is_independent_of_threads() {
        // Two groups of 150 transcripts: components large enough to step in parallel.
        let (ntx, reads) = (301, grouped_reads(200_000, 301, 150));
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), ntx);
        let plan = eqc.plan(ntx);
        assert!(plan.has_large());
        let init = initial(&eqc.counts, ntx, None);
        let on = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| em_components::run(&plan, Counts::Of(&eqc.counts), &init, 300, 1e-6, 1))
        };
        let one = on(1);
        for threads in [3, 8] {
            assert!(
                one.iter()
                    .zip(on(threads))
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
        }
    }

    /// The class EM equals a plain per-read EM on the same reads.
    #[test]
    fn class_em_matches_per_read_em() {
        // Reads as (target, prob) lists; several repeat, in different orders.
        let reads: Vec<Vec<(u32, f32)>> = vec![
            vec![(0, 1.0), (1, 0.5)],
            vec![(1, 0.5), (0, 1.0)],
            vec![(1, 1.0), (2, 1.0)],
            vec![(0, 1.0)],
            vec![(2, 0.25), (0, 1.0), (1, 1.0)],
            vec![(1, 1.0), (2, 1.0)],
            vec![(2, 1.0)],
        ];
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), 3);
        assert_eq!(eqc.len(), 5);
        let tinfo: Vec<TranscriptInfo> = (0..3)
            .map(|_| TranscriptInfo::with_len(std::num::NonZeroUsize::new(1000).unwrap()))
            .collect();
        let got = run(&eqc, &eqc.counts, &tinfo, 200, 1e-12, None, None, false).0;

        // Per-read EM with the same initialization, rounds and final pruning round.
        let mut prev = vec![reads.len() as f64 / 3.0; 3];
        let step = |prev: &[f64]| {
            let mut curr = vec![0.0; 3];
            for r in &reads {
                let d: f64 = r.iter().map(|&(t, p)| prev[t as usize] * p as f64).sum();
                for &(t, p) in r {
                    curr[t as usize] += prev[t as usize] * p as f64 / d;
                }
            }
            curr
        };
        for niter in 0..200u32 {
            let curr = step(&prev);
            let rd = prev
                .iter()
                .zip(&curr)
                .filter(|(p, _)| **p > constants::MIN_READ_THRESH)
                .map(|(p, c)| (c - p) / p)
                .fold(0.0_f64, f64::max);
            prev = curr;
            if rd < 1e-12 && niter > 1 {
                break;
            }
        }
        for x in &mut prev {
            if *x < constants::MIN_READ_THRESH {
                *x = 0.0;
            }
        }
        let want = step(&prev);
        for (g, w) in got.iter().zip(&want) {
            assert!((g - w).abs() < 1e-9, "{got:?} vs {want:?}");
        }
        assert!((got.iter().sum::<f64>() - reads.len() as f64).abs() < 1e-9);
    }

    use crate::em_tests::{
        Classes, MIXED_NTX, Read, RefEm, Rng, assert_close, check_invariants, mixed_reads, on,
        ref_em,
    };
    use crate::presence::PresenceKind;

    fn classes_of(reads: &[Read], ntx: usize) -> EqClasses {
        EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), ntx)
    }

    /// The classes as reference-EM input.
    fn ref_classes(eqc: &EqClasses, counts: &[f64], ntx: usize) -> Classes {
        let reads: Vec<Read> = (0..eqc.len())
            .map(|c| {
                let (ts, ps) = eqc.class(c);
                ts.iter().copied().zip(ps.iter().copied()).collect()
            })
            .collect();
        Classes::with_counts(&reads, counts, ntx)
    }

    /// A read's (target, probability bits) multiset, sorted.
    fn key(r: &[(u32, f32)]) -> Vec<(u32, u32)> {
        let mut k: Vec<(u32, u32)> = r.iter().map(|&(t, p)| (t, p.to_bits())).collect();
        k.sort_unstable();
        k
    }

    /// Reads collapse into one class exactly when their (target, probability-bits) multisets
    /// agree: alignment order does not matter, a probability one ulp away makes a new class,
    /// and a target repeated within a read is kept. Each class's targets are sorted, its count
    /// is its reads, and the transpose lists each transcript's (class, probability) entries.
    #[test]
    fn classes_are_the_exact_multisets_of_reads() {
        let mut rng = Rng::new(31);
        let mut reads = mixed_reads(30);
        for r in reads.iter_mut() {
            rng.shuffle(r);
        }
        let ulp = f32::from_bits(0.5_f32.to_bits() + 1);
        reads.extend([
            vec![(3, 0.5), (4, 1.0)],
            vec![(4, 1.0), (3, 0.5)],
            vec![(3, ulp), (4, 1.0)],
            vec![(3, 0.5), (3, 0.5), (4, 1.0)],
        ]);
        let ntx = MIXED_NTX;
        let eqc = classes_of(&reads, ntx);
        let n = reads.len();
        assert_eq!(eqc.read_class.len(), n);
        assert_eq!(eqc.counts.iter().sum::<f64>(), n as f64);
        let keys: Vec<Vec<(u32, u32)>> = (0..eqc.len())
            .map(|c| {
                let (ts, ps) = eqc.class(c);
                assert!(ts.windows(2).all(|w| w[0] <= w[1]), "class {c} unsorted");
                key(&ts
                    .iter()
                    .copied()
                    .zip(ps.iter().copied())
                    .collect::<Vec<_>>())
            })
            .collect();
        let mut distinct = keys.clone();
        distinct.sort();
        distinct.dedup();
        assert_eq!(distinct.len(), eqc.len(), "two classes share a multiset");
        let mut count = vec![0.0; eqc.len()];
        for (r, read) in reads.iter().enumerate() {
            let c = eqc.read_class[r] as usize;
            assert_eq!(keys[c], key(read), "read {r}");
            count[c] += 1.0;
        }
        assert_eq!(count, eqc.counts);
        let (a, b) = (eqc.read_class[n - 4], eqc.read_class[n - 3]);
        assert_eq!(a, b);
        assert_ne!(eqc.read_class[n - 2], a);
        assert_ne!(eqc.read_class[n - 1], a);
        // The transpose.
        assert_eq!(eqc.t_off.len(), ntx + 1);
        for t in 0..ntx {
            let want: Vec<(u32, u32)> = (0..eqc.len())
                .flat_map(|c| {
                    let (ts, ps) = eqc.class(c);
                    ts.iter()
                        .zip(ps)
                        .filter(|x| *x.0 == t as u32)
                        .map(|x| (c as u32, x.1.to_bits()))
                        .collect::<Vec<_>>()
                })
                .collect();
            let r = eqc.t_off[t]..eqc.t_off[t + 1];
            let got: Vec<(u32, u32)> = eqc.t_class[r.clone()]
                .iter()
                .zip(&eqc.t_prob[r])
                .map(|(&c, p)| (c, p.to_bits()))
                .collect();
            assert_eq!(got, want, "t{t}");
        }
        // The store's classes are those of its reads.
        let h = noodles_sam::Header::default();
        let mut txps = tinfo(ntx);
        let store = crate::em_tests::store_of(&h, &reads, &mut txps);
        let es = EqClasses::from_store(&store, ntx);
        assert_eq!(
            (&es.targets, &es.offsets, &es.counts, &es.read_class),
            (&eqc.targets, &eqc.offsets, &eqc.counts, &eqc.read_class)
        );
        assert_eq!(
            es.probs.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
            eqc.probs.iter().map(|p| p.to_bits()).collect::<Vec<_>>()
        );
    }

    /// The class loop (the path `OARFISH_EM_COMPONENTS=0` selects) and the component path of
    /// `em` take the reference EM's rounds over the classes: fixed round counts (an infinite
    /// or zero threshold), a converged run, and a run from given initial abundances.
    #[test]
    fn class_em_takes_the_reference_rounds() {
        let reads = mixed_reads(32);
        let ntx = MIXED_NTX;
        let eqc = classes_of(&reads, ntx);
        let ti = tinfo(ntx);
        let cl = ref_classes(&eqc, &eqc.counts, ntx);
        let mut rng = Rng::new(33);
        let given: Vec<f64> = (0..ntx).map(|_| 0.1 + rng.below(50) as f64).collect();
        for (max_iter, thresh, init, rtol) in [
            (1000, f64::INFINITY, None, 1e-11),
            (0, 0.0, None, 1e-11),
            (1, 0.0, None, 1e-11),
            (70, 0.0, None, 1e-11),
            (3000, 1e-10, None, 1e-7),
            (3000, 1e-10, Some(&given), 1e-7),
        ] {
            let start = init.cloned().unwrap_or_else(|| cl.uniform());
            let (want, _) = ref_em(&cl, &start, RefEm::new(max_iter, thresh, 1));
            let what = format!("{max_iter} {thresh} {}", init.is_some());
            let got = run(&eqc, &eqc.counts, &ti, max_iter, thresh, init, None, false).0;
            assert_close(&got, &want, rtol, &format!("class loop {what}"));
            let got = em(&eqc, &ti, max_iter, thresh, init, None, 2).0;
            assert_close(&got, &want, rtol, &format!("em {what}"));
        }
        // Per read: the classes' reference is the reads' reference.
        let (a, _) = ref_em(&cl, &cl.uniform(), RefEm::new(3000, 1e-10, 1));
        let rc = Classes::of_reads(&reads, ntx);
        let (b, _) = ref_em(&rc, &rc.uniform(), RefEm::new(3000, 1e-10, 1));
        assert_close(&a, &b, 1e-7, "classes vs reads");
    }

    fn presence(kind: PresenceKind) -> PresenceOpts {
        PresenceOpts {
            kind,
            rho: 0.05,
            warmup: 5,
            period: 5,
        }
    }

    /// With a presence model the EM still assigns every read: estimates are non-negative and
    /// sum to the read count, transcripts without reads get 0 and those with only
    /// single-target reads get them; the result does not depend on the number of threads
    /// beyond summation order.
    #[test]
    fn presence_em_keeps_the_invariants() {
        let reads = mixed_reads(34);
        let ntx = MIXED_NTX;
        let eqc = classes_of(&reads, ntx);
        let ti = tinfo(ntx);
        for kind in [PresenceKind::SpikeSlab, PresenceKind::Bayes] {
            let p = Some(presence(kind));
            let one = em(&eqc, &ti, 500, 1e-8, None, p, 1).0;
            check_invariants(&one, &reads, None, ntx, &format!("{kind:?}"));
            let four = em(&eqc, &ti, 500, 1e-8, None, p, 4).0;
            assert_close(&four, &one, 1e-7, &format!("{kind:?} threads"));
            let direct = on(4, || {
                run(&eqc, &eqc.counts, &ti, 500, 1e-8, None, p, false).0
            });
            assert_close(&direct, &one, 1e-7, &format!("{kind:?} class loop"));
            // Every read single-target: each transcript gets its reads, whatever q.
            let singles: Vec<Read> = (0..40).map(|i| vec![((i % 7) as u32, 0.5)]).collect();
            let es = classes_of(&singles, 9);
            let got = em(&es, &tinfo(9), 500, 1e-8, None, p, 2).0;
            let want: Vec<f64> = (0..9)
                .map(|t| {
                    if t < 7 {
                        [6.0, 6.0, 6.0, 6.0, 6.0, 5.0, 5.0][t]
                    } else {
                        0.0
                    }
                })
                .collect();
            assert_close(&got, &want, 1e-9, &format!("{kind:?} singles"));
        }
    }

    /// The class counts of a resample are the duplicated reads: the class loop and the
    /// component layout on them take the reference rounds on the duplicated reads; and every
    /// bootstrap replicate, with or without a presence model, keeps the invariants.
    #[test]
    fn bootstrap_counts_are_duplicated_reads() {
        let reads = mixed_reads(35);
        let ntx = MIXED_NTX;
        let eqc = classes_of(&reads, ntx);
        let ti = tinfo(ntx);
        let mult = Rng::new(36).resample(reads.len());
        let mut counts = vec![0.0; eqc.len()];
        for (r, &m) in mult.iter().enumerate() {
            counts[eqc.read_class[r] as usize] += m as f64;
        }
        let dup: Vec<Read> = reads
            .iter()
            .zip(&mult)
            .flat_map(|(r, &m)| std::iter::repeat_n(r.clone(), m as usize))
            .collect();
        let dc = Classes::of_reads(&dup, ntx);
        let (want, _) = ref_em(&dc, &dc.uniform(), RefEm::new(3000, 1e-10, 1));
        let got = run(&eqc, &counts, &ti, 3000, 1e-10, None, None, false).0;
        assert_close(&got, &want, 1e-7, "class loop");
        let init = initial(&counts, ntx, None);
        let got = em_components::run(&eqc.plan(ntx), Counts::Of(&counts), &init, 3000, 1e-10, 1);
        assert_close(&got, &want, 1e-7, "components");
        for p in [None, Some(presence(PresenceKind::SpikeSlab))] {
            let reps = bootstrap(&eqc, &ti, 300, 1e-6, None, p, 2, 2);
            assert_eq!(reps.len(), 2);
            for b in reps {
                assert!(b.iter().all(|&x| x.is_finite() && x >= 0.0));
                let sum: f64 = b.iter().sum();
                assert!((sum - reads.len() as f64).abs() < 1e-9 * reads.len() as f64);
                assert_eq!((b[62], b[63]), (0.0, 0.0));
            }
        }
    }

    /// A transcript that starts at 0 (no short-read abundance) keeps the reads only it explains
    /// once the start is floored, and the class loop and the component layout then agree.
    #[test]
    fn zero_start_keeps_unique_reads() {
        // 5 reads unique to t0; t1 and t2 share 4 reads, and t2 has 3 of its own.
        let mut reads: Vec<Vec<(u32, f32)>> = vec![vec![(0, 1.0)]; 5];
        reads.extend(vec![vec![(1, 1.0), (2, 1.0)]; 4]);
        reads.extend(vec![vec![(2, 1.0)]; 3]);
        let eqc = EqClasses::from_reads(reads.iter().map(|r| r.iter().copied()), reads.len(), 3);
        let mut init = vec![0.0, 3.0, 3.0];
        crate::util::read_function::floor_zero_abundances(&mut init);
        let want = run(
            &eqc,
            &eqc.counts,
            &tinfo(3),
            1000,
            1e-9,
            Some(&init),
            None,
            false,
        )
        .0;
        let got = em_components::run(&eqc.plan(3), Counts::Of(&eqc.counts), &init, 1000, 1e-9, 1);
        for est in [&want, &got] {
            assert!(
                (est.iter().sum::<f64>() - reads.len() as f64).abs() < 1e-9,
                "{est:?}"
            );
            assert!((est[0] - 5.0).abs() < 1e-9, "{est:?}");
        }
        for (g, w) in got.iter().zip(&want) {
            assert!((g - w).abs() <= 1e-6 * w.max(1.0), "{got:?} vs {want:?}");
        }
    }
}
