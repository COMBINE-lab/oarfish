//! Fixtures shared by the EM test suites, a dense reference EM, and tests that run every EM
//! entry point on the same inputs.
//!
//! The reference EM is a plain loop over classes and their targets with the stopping rule the
//! EMs document: after each round, the largest relative increase over transcripts above
//! [`constants::MIN_READ_THRESH`] (or, for the truncation EM's frozen loop, the largest
//! relative change in either direction); stop after the first round whose statistic is below the
//! threshold once the round's 0-based index exceeds `min_iter`, or after `max_iter` rounds; then
//! zero abundances below [`constants::MIN_READ_THRESH`] and run one more round.
//!
//! Expected values come from the reference, from closed forms, or from invariants (mass
//! conservation, non-negativity, transcripts no read touches get exactly 0, transcripts that
//! only single-target reads touch get their read count). Tolerances allow for any summation
//! order; only thread-independence of the component layouts is checked bit for bit.

use std::num::NonZeroUsize;

use bio_types::strand::Strand;
use noodles_sam::Header;

use crate::em_components::{self, Counts, Plan};
use crate::em_read::{self, Impl};
use crate::eq_classes::{self, EqClasses};
use crate::util::constants;
use crate::util::oarfish_types::{
    AlignmentFilters, AlnInfo, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};

/// A fixed-seed xorshift generator, so fixtures are the same on every run and platform.
pub(crate) struct Rng(u64);

impl Rng {
    pub(crate) fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1)
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    pub(crate) fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    pub(crate) fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i + 1);
            v.swap(i, j);
        }
    }

    /// A random permutation of `0..n`.
    pub(crate) fn perm(&mut self, n: usize) -> Vec<u32> {
        let mut p: Vec<u32> = (0..n as u32).collect();
        self.shuffle(&mut p);
        p
    }

    /// Multiplicities of `n` draws with replacement from `0..n` (they sum to `n`).
    pub(crate) fn resample(&mut self, n: usize) -> Vec<u32> {
        let mut m = vec![0_u32; n];
        for _ in 0..n {
            m[self.below(n)] += 1;
        }
        m
    }
}

/// A read: its (target, alignment-score probability) pairs.
pub(crate) type Read = Vec<(u32, f32)>;

/// Probabilities fixture reads draw from (few, so that reads repeat as classes).
pub(crate) const PROBS: [f32; 5] = [1.0, 0.8, 0.5, 0.25, 0.1];

/// `nreads` reads over `ngroups` groups of `group` consecutive transcripts: each read has 1-4
/// targets from one group, probabilities from [`PROBS`]. With `dups`, a read may align to one
/// target more than once (as a read with two alignments to one transcript does).
pub(crate) fn grouped_reads(
    rng: &mut Rng,
    nreads: usize,
    ngroups: usize,
    group: usize,
    dups: bool,
) -> Vec<Read> {
    (0..nreads)
        .map(|_| {
            let g = rng.below(ngroups);
            let k = 1 + rng.below(4);
            let mut r: Read = (0..k)
                .map(|_| {
                    (
                        (g * group + rng.below(group)) as u32,
                        PROBS[rng.below(PROBS.len())],
                    )
                })
                .collect();
            if !dups {
                r.sort_by_key(|x| x.0);
                r.dedup_by_key(|x| x.0);
            }
            r
        })
        .collect()
}

/// Each read `n` times.
pub(crate) fn repeated(reads: &[(Read, usize)]) -> Vec<Read> {
    reads
        .iter()
        .flat_map(|(r, n)| std::iter::repeat_n(r.clone(), *n))
        .collect()
}

/// Number of transcripts of [`mixed_reads`].
pub(crate) const MIXED_NTX: usize = 64;

/// A mixed fixture over [`MIXED_NTX`] transcripts: 1,500 grouped reads over transcripts 0..60
/// (10 groups of 6, with repeated targets), 5 and 9 single-target reads on transcripts 60 and
/// 61 (which no other read touches), nothing on 62 and 63; shuffled.
pub(crate) fn mixed_reads(seed: u64) -> Vec<Read> {
    let mut rng = Rng::new(seed);
    let mut reads = grouped_reads(&mut rng, 1500, 10, 6, true);
    reads.extend((0..5).map(|_| vec![(60, 0.8)]));
    reads.extend((0..9).map(|_| vec![(61, 1.0)]));
    rng.shuffle(&mut reads);
    reads
}

/// Classes (entries and read count) over `ntx` transcripts, for the reference EM.
pub(crate) struct Classes {
    pub(crate) ntx: usize,
    pub(crate) cls: Vec<(Vec<(u32, f64)>, f64)>,
}

impl Classes {
    /// Each read its own class of count 1.
    pub(crate) fn of_reads(reads: &[Read], ntx: usize) -> Self {
        Self::with_counts(reads, &vec![1.0; reads.len()], ntx)
    }

    /// Read `r` with count `counts[r]`.
    pub(crate) fn with_counts(reads: &[Read], counts: &[f64], ntx: usize) -> Self {
        Classes {
            ntx,
            cls: reads
                .iter()
                .zip(counts)
                .map(|(r, &n)| (r.iter().map(|&(t, p)| (t, p as f64)).collect(), n))
                .collect(),
        }
    }

    pub(crate) fn total(&self) -> f64 {
        self.cls.iter().map(|c| c.1).sum()
    }

    /// The even initial abundances every EM uses by default.
    pub(crate) fn uniform(&self) -> Vec<f64> {
        vec![self.total() / self.ntx as f64; self.ntx]
    }

    /// One EM round: each class's reads spread over its targets in proportion to
    /// `prev[t] * w`.
    pub(crate) fn step(&self, prev: &[f64]) -> Vec<f64> {
        let mut curr = vec![0.0; self.ntx];
        for (e, n) in &self.cls {
            if *n == 0.0 {
                continue;
            }
            let d: f64 = e.iter().map(|&(t, w)| prev[t as usize] * w).sum();
            if d > constants::EM_DENOM_THRESH {
                for &(t, w) in e {
                    curr[t as usize] += n * prev[t as usize] * w / d;
                }
            }
        }
        curr
    }

    /// The log-likelihood, up to a constant, of abundances `theta` normalized to proportions.
    pub(crate) fn loglik(&self, theta: &[f64]) -> f64 {
        let total: f64 = theta.iter().sum();
        let (mut ll, mut reads) = (0.0, 0.0);
        for (e, n) in &self.cls {
            let d: f64 = e.iter().map(|&(t, w)| theta[t as usize] * w).sum();
            if *n > 0.0 && d > constants::EM_DENOM_THRESH {
                ll += n * d.ln();
                reads += n;
            }
        }
        ll - reads * total.ln()
    }
}

/// The reference EM's options.
#[derive(Clone, Copy)]
pub(crate) struct RefEm {
    pub(crate) max_iter: u32,
    pub(crate) thresh: f64,
    pub(crate) min_iter: u32,
    /// Iteration count the loop starts from.
    pub(crate) niter0: u32,
    /// Relative change in either direction (else relative increase only).
    pub(crate) abs: bool,
}

impl RefEm {
    pub(crate) fn new(max_iter: u32, thresh: f64, min_iter: u32) -> Self {
        RefEm {
            max_iter,
            thresh,
            min_iter,
            niter0: 0,
            abs: false,
        }
    }
}

/// The reference EM from `init`; returns the estimates and the iteration count at exit.
pub(crate) fn ref_em(c: &Classes, init: &[f64], o: RefEm) -> (Vec<f64>, u32) {
    let mut prev = init.to_vec();
    let mut niter = o.niter0;
    while niter < o.max_iter {
        let curr = c.step(&prev);
        let mut rd = 0.0_f64;
        for (&p, &x) in prev.iter().zip(&curr) {
            if p > constants::MIN_READ_THRESH {
                let d = if o.abs { (x - p).abs() } else { x - p };
                rd = rd.max(d / p);
            }
        }
        prev = curr;
        if rd < o.thresh && niter > o.min_iter {
            break;
        }
        niter += 1;
    }
    for x in &mut prev {
        if *x < constants::MIN_READ_THRESH {
            *x = 0.0;
        }
    }
    (c.step(&prev), niter)
}

/// Asserts `got` and `want` agree to `rtol` relative to `max(|want|, 1)`.
#[track_caller]
pub(crate) fn assert_close(got: &[f64], want: &[f64], rtol: f64, what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: lengths");
    for (t, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(
            (g - w).abs() <= rtol * w.abs().max(1.0),
            "{what}: t{t}: {g} vs {w}"
        );
    }
}

#[track_caller]
pub(crate) fn assert_bits_eq(a: &[f64], b: &[f64], what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: lengths");
    for (t, (x, y)) in a.iter().zip(b).enumerate() {
        assert_eq!(x.to_bits(), y.to_bits(), "{what}: t{t}: {x} vs {y}");
    }
}

/// Checks what every EM guarantees of estimates for `reads` (each counted `mult[r]` times, or
/// once): they are finite and non-negative and sum to the read count; transcripts no read
/// touches get exactly 0; transcripts only single-target reads touch get those reads.
#[track_caller]
pub(crate) fn check_invariants(
    est: &[f64],
    reads: &[Read],
    mult: Option<&[u32]>,
    ntx: usize,
    what: &str,
) {
    assert_eq!(est.len(), ntx, "{what}: length");
    assert!(
        est.iter().all(|&x| x.is_finite() && x >= 0.0),
        "{what}: {est:?}"
    );
    let m = |r: usize| mult.map_or(1.0, |m| m[r] as f64);
    let n: f64 = (0..reads.len()).map(m).sum();
    let sum: f64 = est.iter().sum();
    assert!(
        (sum - n).abs() <= 1e-9 * n.max(1.0),
        "{what}: total {sum} vs {n} reads"
    );
    let (mut seen, mut shared) = (vec![false; ntx], vec![false; ntx]);
    let mut single = vec![0.0; ntx];
    for (r, read) in reads.iter().enumerate() {
        let t0 = read[0].0;
        let one = read.iter().all(|&(t, _)| t == t0);
        for &(t, _) in read {
            seen[t as usize] = true;
            shared[t as usize] |= !one;
        }
        if one {
            single[t0 as usize] += m(r);
        }
    }
    for t in 0..ntx {
        if !seen[t] {
            assert_eq!(est[t], 0.0, "{what}: t{t} has no reads");
        } else if !shared[t] {
            assert!(
                (est[t] - single[t]).abs() <= 1e-9 * single[t].max(1.0),
                "{what}: t{t}: {} vs its {} single-target reads",
                est[t],
                single[t]
            );
        }
    }
}

/// Filters that keep everything (the EMs never look at them but for `model_coverage`).
pub(crate) fn filters(model_coverage: bool) -> AlignmentFilters {
    AlignmentFilters::builder()
        .five_prime_clip(u32::MAX)
        .three_prime_clip(i64::MAX)
        .score_threshold(0.0)
        .min_aligned_fraction(0.0)
        .min_aligned_len(1)
        .which_strand(Strand::Unknown)
        .model_coverage(model_coverage)
        .logistic_growth_rate(2.0)
        .write_assignment_probs(false)
        .write_assignment_probs_type(None)
        .build()
}

/// Transcripts of the given lengths.
pub(crate) fn tinfo(lens: &[usize]) -> Vec<TranscriptInfo> {
    lens.iter()
        .map(|&l| TranscriptInfo::with_len(NonZeroUsize::new(l).unwrap()))
        .collect()
}

/// A store holding `reads`, each alignment given as (alignment, probability).
pub(crate) fn store_of_alns<'h>(
    h: &'h Header,
    reads: &[Vec<(AlnInfo, f32)>],
    txps: &mut [TranscriptInfo],
    model_coverage: bool,
) -> InMemoryAlignmentStore<'h> {
    let mut store = InMemoryAlignmentStore::new(filters(model_coverage), h);
    for r in reads {
        let alns: Vec<AlnInfo> = r.iter().map(|x| x.0).collect();
        let probs: Vec<f32> = r.iter().map(|x| x.1).collect();
        assert!(store.add_filtered_group(&alns, &probs, txps));
    }
    store
}

/// A store holding `reads` (every alignment spanning 0..100 of its target).
pub(crate) fn store_of<'h>(
    h: &'h Header,
    reads: &[Read],
    txps: &mut [TranscriptInfo],
) -> InMemoryAlignmentStore<'h> {
    let alns: Vec<Vec<(AlnInfo, f32)>> = reads
        .iter()
        .map(|r| {
            r.iter()
                .map(|&(t, p)| {
                    (
                        AlnInfo {
                            ref_id: t,
                            start: 0,
                            end: 100,
                            strand: Strand::Forward,
                        },
                        p,
                    )
                })
                .collect()
        })
        .collect();
    store_of_alns(h, &alns, txps, false)
}

pub(crate) fn em_info<'a, 'b, 'h>(
    store: &'a InMemoryAlignmentStore<'h>,
    txps: &'b [TranscriptInfo],
    max_iter: u32,
    thresh: f64,
) -> EMInfo<'a, 'b, 'h> {
    EMInfo {
        eq_map: store,
        txp_info: txps,
        max_iter,
        convergence_thresh: thresh,
        init_abundances: None,
        kde_model: None,
    }
}

/// Runs `f` on a pool of `threads` threads.
pub(crate) fn on<T: Send>(threads: usize, f: impl FnOnce() -> T + Send) -> T {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(f)
}

/// The component layout of a store's reads, as the per-read EM builds it (without coverage or
/// length terms).
fn read_plan(store: &InMemoryAlignmentStore, ntx: usize, by_heaviest: bool) -> Plan<f64> {
    Plan::new(
        |k| store.alignments[k].ref_id,
        |k| store.as_probabilities[k] as f64,
        store.boundaries(),
        ntx,
        by_heaviest,
    )
}

/// Every plain-EM entry point on one store, named, with the `min_iter` its loop uses. The
/// equivalence-class EM, the per-read EMs and the component layout should all run the
/// reference's rounds.
fn every_path(emi: &EMInfo, ntx: usize) -> Vec<(String, u32, Vec<f64>)> {
    let store = emi.eq_map;
    let (max_iter, thresh) = (emi.max_iter, emi.convergence_thresh);
    let eqc = EqClasses::from_store(store, ntx);
    let mut out = vec![
        ("em::em".to_string(), 50, crate::em::em(emi, 1)),
        ("em::em_par(4)".to_string(), 1, crate::em::em_par(emi, 4)),
    ];
    for threads in [1, 4] {
        out.push((
            format!("em_read Local({threads})"),
            1,
            em_read::em(emi, threads, Impl::Local),
        ));
    }
    for threads in [1, 4, 16] {
        out.push((
            format!("em_read Csr({threads})"),
            1,
            em_read::em(emi, threads, Impl::Csr),
        ));
        out.push((
            format!("eq_classes::em({threads})"),
            1,
            eq_classes::em(&eqc, emi.txp_info, max_iter, thresh, None, None, threads),
        ));
    }
    for by_heaviest in [false, true] {
        let plan = read_plan(store, ntx, by_heaviest);
        let init = vec![store.len() as f64 / ntx as f64; ntx];
        out.push((
            format!("em_components::run(by_heaviest {by_heaviest})"),
            1,
            em_components::run(&plan, Counts::Unit, &init, max_iter, thresh, 1),
        ));
    }
    out
}

/// Every EM path reaches the reference EM's estimates (run with the path's own `min_iter`),
/// and each satisfies the invariants.
#[test]
fn every_em_path_reaches_the_reference_estimates() {
    let reads = mixed_reads(1);
    let ntx = MIXED_NTX;
    let h = Header::default();
    let mut txps = tinfo(&vec![1000; ntx]);
    let store = store_of(&h, &reads, &mut txps);
    let emi = em_info(&store, &txps, 3000, 1e-10);
    let c = Classes::of_reads(&reads, ntx);
    let want = |min_iter| ref_em(&c, &c.uniform(), RefEm::new(3000, 1e-10, min_iter)).0;
    let (want1, want50) = (want(1), want(50));
    check_invariants(&want1, &reads, None, ntx, "reference");
    for (name, min_iter, got) in every_path(&emi, ntx) {
        let want = if min_iter == 50 { &want50 } else { &want1 };
        assert_close(&got, want, 1e-7, &name);
        check_invariants(&got, &reads, None, ntx, &name);
    }
}

/// Relabelling the transcripts and reordering the reads (and each read's alignments) relabels
/// the estimates.
#[test]
fn em_paths_are_permutation_invariant() {
    let reads = mixed_reads(2);
    let ntx = MIXED_NTX;
    let mut rng = Rng::new(3);
    let pi = rng.perm(ntx);
    let mut moved: Vec<Read> = reads
        .iter()
        .map(|r| {
            let mut r: Read = r.iter().map(|&(t, p)| (pi[t as usize], p)).collect();
            rng.shuffle(&mut r);
            r
        })
        .collect();
    rng.shuffle(&mut moved);
    let h = Header::default();
    let (mut ta, mut tb) = (tinfo(&vec![1000; ntx]), tinfo(&vec![1000; ntx]));
    let (sa, sb) = (store_of(&h, &reads, &mut ta), store_of(&h, &moved, &mut tb));
    let (ea, eb) = (
        em_info(&sa, &ta, 3000, 1e-10),
        em_info(&sb, &tb, 3000, 1e-10),
    );
    for ((name, _, a), (_, _, b)) in every_path(&ea, ntx).into_iter().zip(every_path(&eb, ntx)) {
        let relabelled: Vec<f64> = (0..ntx).map(|t| b[pi[t] as usize]).collect();
        assert_close(&relabelled, &a, 1e-7, &name);
    }
}

/// The component EM and the transposed layouts write each slot once and sum in a fixed order:
/// the per-read EM (`Csr`) and the equivalence-class EM give the same bits on any number of
/// threads.
#[test]
fn csr_and_class_em_are_bitwise_independent_of_threads() {
    let reads = mixed_reads(4);
    let ntx = MIXED_NTX;
    let h = Header::default();
    let mut txps = tinfo(&vec![1000; ntx]);
    let store = store_of(&h, &reads, &mut txps);
    let emi = em_info(&store, &txps, 500, 1e-8);
    let eqc = EqClasses::from_store(&store, ntx);
    let read1 = em_read::em(&emi, 1, Impl::Csr);
    let cls1 = eq_classes::em(&eqc, &txps, 500, 1e-8, None, None, 1);
    for threads in [4, 16] {
        assert_bits_eq(
            &em_read::em(&emi, threads, Impl::Csr),
            &read1,
            &format!("em_read Csr {threads}"),
        );
        assert_bits_eq(
            &eq_classes::em(&eqc, &txps, 500, 1e-8, None, None, threads),
            &cls1,
            &format!("eq_classes {threads}"),
        );
    }
}

/// Two transcripts share a class: with `a` reads unique to t0, `b` unique to t1 and `c` shared
/// at equal weights, the fixed point is `theta0 = a N / (a + b)`, `theta1 = b N / (a + b)`
/// (N = a + b + c). With `a` reads unique to t0 and `c` shared at weights (w0, w1), w1 > w0, the
/// fixed point with both present solves `w0 theta0 + w1 theta1 = c w1`:
/// `theta1 = (c w1 - w0 (a + c)) / (w1 - w0)`.
#[test]
fn two_transcripts_sharing_a_class_reach_the_closed_form() {
    let h = Header::default();
    let cases: Vec<(Vec<Read>, [f64; 2])> = vec![
        (
            repeated(&[
                (vec![(0, 1.0)], 3),
                (vec![(1, 1.0)], 1),
                (vec![(0, 1.0), (1, 1.0)], 4),
            ]),
            [6.0, 2.0],
        ),
        (
            repeated(&[(vec![(0, 1.0)], 2), (vec![(0, 0.5), (1, 1.0)], 10)]),
            [4.0, 8.0],
        ),
    ];
    for (reads, want) in cases {
        let mut txps = tinfo(&[1000, 1000]);
        let store = store_of(&h, &reads, &mut txps);
        let emi = em_info(&store, &txps, 10_000, 1e-14);
        for (name, _, got) in every_path(&emi, 2) {
            assert_close(&got, &want, 1e-8, &name);
        }
        let plan = read_plan(&store, 2, true);
        let init = vec![reads.len() as f64 / 2.0; 2];
        for plain_first in [0, 5, u32::MAX] {
            let (got, _) = em_components::run_squarem(
                &plan,
                Counts::Unit,
                &init,
                0,
                10_000,
                1e-14,
                plain_first,
            );
            assert_close(&got, &want, 1e-8, &format!("run_squarem {plain_first}"));
        }
    }
}

/// A degenerate input: its name, transcript count, reads and expected estimates.
type Case = (&'static str, usize, Vec<Read>, Option<Vec<f64>>);

/// Degenerate inputs: no reads; one transcript; one read (its mass goes to its best target);
/// every read in one class (equal weights keep the even split); every class single-target.
#[test]
fn degenerate_inputs() {
    let h = Header::default();
    // (name, transcripts, reads, expected estimates)
    let cases: Vec<Case> = vec![
        ("no reads", 3, vec![], Some(vec![0.0; 3])),
        (
            "one transcript",
            1,
            (0..50).map(|i| vec![(0, PROBS[i % 5])]).collect(),
            Some(vec![50.0]),
        ),
        (
            "one read",
            4,
            vec![vec![(2, 0.5), (0, 1.0), (1, 0.25)]],
            Some(vec![1.0, 0.0, 0.0, 0.0]),
        ),
        (
            "one class",
            4,
            vec![vec![(0, 0.5), (1, 0.5), (3, 0.5)]; 40],
            Some(vec![40.0 / 3.0, 40.0 / 3.0, 0.0, 40.0 / 3.0]),
        ),
        (
            "all singles",
            5,
            [(0, 3), (1, 1), (3, 7), (4, 2)]
                .into_iter()
                .flat_map(|(t, n)| std::iter::repeat_n(vec![(t, 0.5)], n))
                .collect(),
            Some(vec![3.0, 1.0, 0.0, 7.0, 2.0]),
        ),
    ];
    for (what, ntx, reads, want) in cases {
        let mut txps = tinfo(&vec![1000; ntx]);
        let store = store_of(&h, &reads, &mut txps);
        let emi = em_info(&store, &txps, 10_000, 1e-12);
        for (name, _, got) in every_path(&emi, ntx) {
            let name = format!("{what}: {name}");
            check_invariants(&got, &reads, None, ntx, &name);
            if let Some(want) = &want {
                assert_close(&got, want, 1e-8, &name);
            }
        }
    }
}

/// A multiplicity of k is k copies of the read, in the component EM (`Counts::Of`,
/// `Counts::OfU32`), in the score-only EM over a resampling index list, and in the per-read EM
/// over a store holding the copies.
#[test]
fn multiplicities_equal_duplicated_reads() {
    let reads = mixed_reads(5);
    let ntx = MIXED_NTX;
    let mult = Rng::new(6).resample(reads.len());
    let dup: Vec<Read> = reads
        .iter()
        .zip(&mult)
        .flat_map(|(r, &m)| std::iter::repeat_n(r.clone(), m as usize))
        .collect();
    let h = Header::default();
    let mut txps = tinfo(&vec![1000; ntx]);
    let store = store_of(&h, &reads, &mut txps);
    let dstore = store_of(&h, &dup, &mut txps);
    assert_eq!(dstore.len(), store.len());
    let c = Classes::of_reads(&dup, ntx);
    let init = c.uniform();
    let multf: Vec<f64> = mult.iter().map(|&m| m as f64).collect();
    let plan = read_plan(&store, ntx, true);
    for min_iter in [1, 50] {
        let (want, _) = ref_em(&c, &init, RefEm::new(3000, 1e-10, min_iter));
        check_invariants(&want, &reads, Some(&mult), ntx, "reference");
        for (name, n) in [("Of", Counts::Of(&multf)), ("OfU32", Counts::OfU32(&mult))] {
            let got = em_components::run(&plan, n, &init, 3000, 1e-10, min_iter);
            assert_close(&got, &want, 1e-7, &format!("{name} min_iter {min_iter}"));
        }
        let demi = em_info(&dstore, &txps, 3000, 1e-10);
        let got = if min_iter == 50 {
            crate::em::em(&demi, 1)
        } else {
            em_read::em(&demi, 2, Impl::Csr)
        };
        assert_close(&got, &want, 1e-7, &format!("duplicated store {min_iter}"));
    }
    // The score-only EM's bootstrap iterates over drawn read indices.
    let emi = em_info(&store, &txps, 3000, 1e-10);
    let inds: Vec<usize> = (0..reads.len())
        .flat_map(|r| std::iter::repeat_n(r, mult[r] as usize))
        .collect();
    let got = crate::em::do_em(&emi, || store.random_sampling_iter(&inds), false);
    let (want, _) = ref_em(&c, &init, RefEm::new(3000, 1e-10, 50));
    assert_close(&got, &want, 1e-7, "do_em over drawn indices");
}

/// The per-read bootstrap's EM with every read drawn once is the point estimate: the component
/// run with all-ones multiplicities matches the unit-count run, and at `min_iter` 50 the
/// score-only EM.
#[test]
fn all_ones_multiplicities_give_the_point_estimate() {
    let reads = mixed_reads(7);
    let ntx = MIXED_NTX;
    let h = Header::default();
    let mut txps = tinfo(&vec![1000; ntx]);
    let store = store_of(&h, &reads, &mut txps);
    let emi = em_info(&store, &txps, 3000, 1e-10);
    let plan = read_plan(&store, ntx, true);
    let init = vec![reads.len() as f64 / ntx as f64; ntx];
    let (ones, ones32) = (vec![1.0; reads.len()], vec![1_u32; reads.len()]);
    for min_iter in [1, 50] {
        let unit = em_components::run(&plan, Counts::Unit, &init, 3000, 1e-10, min_iter);
        for n in [Counts::Of(&ones), Counts::OfU32(&ones32)] {
            let got = em_components::run(&plan, n, &init, 3000, 1e-10, min_iter);
            assert_close(&got, &unit, 1e-12, "all-ones multiplicities");
        }
        let point = if min_iter == 50 {
            crate::em::em(&emi, 1)
        } else {
            em_read::em(&emi, 4, Impl::Csr)
        };
        assert_close(&unit, &point, 1e-7, &format!("point estimate {min_iter}"));
    }
}

/// Every bootstrap replicate is an EM on a resample of the reads: finite, non-negative, summing
/// to the read count, 0 where no read lies, and a whole number of reads on transcripts only
/// single-target reads touch.
#[test]
fn bootstrap_replicates_conserve_reads() {
    let reads = mixed_reads(8);
    let ntx = MIXED_NTX;
    let h = Header::default();
    let mut txps = tinfo(&vec![1000; ntx]);
    let store = store_of(&h, &reads, &mut txps);
    let emi = em_info(&store, &txps, 200, 1e-6);
    let eqc = EqClasses::from_store(&store, ntx);
    let reps = [
        (
            "eq_classes",
            eq_classes::bootstrap(&eqc, &txps, 200, 1e-6, None, None, 3, 4),
        ),
        ("em_read", em_read::bootstrap(&emi, 3, 4)),
        ("em", crate::em::bootstrap(&emi, 2, 2)),
    ];
    for (name, reps) in reps {
        assert!(reps.len() >= 2, "{name}");
        for b in reps {
            assert_eq!(b.len(), ntx);
            assert!(b.iter().all(|&x| x.is_finite() && x >= 0.0), "{name}");
            let sum: f64 = b.iter().sum();
            assert!(
                (sum - reads.len() as f64).abs() < 1e-9 * reads.len() as f64,
                "{name}: {sum}"
            );
            assert_eq!((b[62], b[63]), (0.0, 0.0), "{name}");
            for t in [60, 61] {
                assert!((b[t] - b[t].round()).abs() < 1e-9, "{name}: {}", b[t]);
            }
        }
    }
}

/// Reordering the store's reads (as `bulk` does for the per-read EMs and the truncation EM)
/// leaves every EM's estimates unchanged.
#[test]
fn reordering_the_store_leaves_the_estimates() {
    let reads = mixed_reads(9);
    let ntx = MIXED_NTX;
    let h = Header::default();
    let mut txps = tinfo(&vec![1000; ntx]);
    let store = store_of(&h, &reads, &mut txps);
    let mut moved = store_of(&h, &reads, &mut txps);
    let order = em_read::read_order(moved.len(), |r| {
        let (alns, probs, _) = moved.read(r);
        alns.iter()
            .zip(probs)
            .max_by(|x, y| x.1.total_cmp(y.1))
            .map_or(0, |x| x.0.ref_id)
    });
    moved.reorder_reads(&order);
    let a = every_path(&em_info(&store, &txps, 3000, 1e-10), ntx);
    let b = every_path(&em_info(&moved, &txps, 3000, 1e-10), ntx);
    for ((name, _, a), (_, _, b)) in a.into_iter().zip(b) {
        assert_close(&b, &a, 1e-7, &name);
    }
}
