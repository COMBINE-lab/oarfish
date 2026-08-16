use crate::util::constants;
use crate::util::oarfish_types::{AlnInfo, EMInfo};
use itertools::izip;
use num_format::{Locale, ToFormattedString};
use rand::rng as trng;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use rayon::prelude::IntoParallelRefMutIterator;
use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use tracing::{info, span, trace};

use crate::bootstrap;

struct PreparedEq<'a> {
    alignments: &'a [AlnInfo],
    weight_start: usize,
    multiplicity: u32,
}

/// Per-class novel-state data, held in an array parallel to the prepared
/// classes rather than inline. Carrying these 16 bytes inside `PreparedEq`
/// would cost cache footprint in the M-step hot loop on every run, including
/// the overwhelmingly common one where no novel state exists; the array is left
/// empty unless `--model-unannotated-isoforms` created states.
#[derive(Clone, Copy)]
struct NovelSlot {
    /// Index of this read's novel latent state (past the annotated
    /// transcripts), or `usize::MAX` when the read shows no annotation-gap
    /// evidence.
    id: usize,
    /// Weight of the novel candidate, on the same scale as the annotated ones.
    weight: f64,
}

const NO_NOVEL: NovelSlot = NovelSlot {
    id: usize::MAX,
    weight: 0.0,
};

fn prepare_equivalence_classes<'a>(
    em_info: &'a EMInfo<'_, '_, '_>,
    collapse_exact: bool,
) -> (Vec<PreparedEq<'a>>, Vec<f64>, Vec<NovelSlot>) {
    let model_coverage = em_info.eq_map.filter_opts.model_coverage;
    let density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.0,
        }
    };
    let mut weights = Vec::with_capacity(em_info.eq_map.alignments.len());
    let mut classes: Vec<PreparedEq<'a>> = Vec::with_capacity(em_info.eq_map.len());
    let mut by_hash: HashMap<u64, Vec<usize>> = HashMap::new();
    let n_txps = em_info.txp_info.len();
    // Only materialised when the feature actually created states.
    let modelling_novel = em_info.novel_loci > 0;
    let mut novel: Vec<NovelSlot> = if modelling_novel {
        Vec::with_capacity(em_info.eq_map.len())
    } else {
        Vec::new()
    };
    for (read_index, (alns, probs, coverage_probs)) in em_info.eq_map.iter().enumerate() {
        let weight_start = weights.len();
        weights.extend(izip!(alns, probs, coverage_probs).map(|(a, p, cp)| {
            let target_id = a.ref_id as usize;
            let txp_len = em_info.txp_info[target_id].lenf as usize;
            let cov_prob = if model_coverage { *cp } else { 1.0 };
            *p as f64 * cov_prob * density_fn(txp_len, a.alignment_span() as usize)
        }));
        // A read whose splice structure disagrees with every annotated
        // candidate gets an extra candidate: a hypothetical transcript at this
        // locus that would match. Its weight is the best annotated weight
        // scaled by the odds of each unmatched junction, so the EM can decline
        // to force the read onto an annotated transcript.
        let slot = if !modelling_novel {
            NO_NOVEL
        } else {
            match em_info.novel_locus.get(read_index) {
                Some(&locus) if locus >= 0 => {
                    let best = weights[weight_start..]
                        .iter()
                        .copied()
                        .fold(0.0_f64, f64::max);
                    let misses = em_info
                        .eq_map
                        .min_junc_misses
                        .get(read_index)
                        .copied()
                        .unwrap_or(0)
                        .max(1) as f64;
                    NovelSlot {
                        id: n_txps + locus as usize,
                        weight: best * em_info.novel_odds_per_miss.powf(misses),
                    }
                }
                _ => NO_NOVEL,
            }
        };
        if collapse_exact {
            let local_weights = &weights[weight_start..];
            let mut hasher = DefaultHasher::new();
            alns.len().hash(&mut hasher);
            for (alignment, weight) in alns.iter().zip(local_weights) {
                alignment.ref_id.hash(&mut hasher);
                weight.to_bits().hash(&mut hasher);
            }
            if modelling_novel {
                slot.id.hash(&mut hasher);
                slot.weight.to_bits().hash(&mut hasher);
            }
            let hash = hasher.finish();
            let duplicate = by_hash.get(&hash).and_then(|indices| {
                indices.iter().copied().find(|&index| {
                    let existing = &classes[index];
                    (!modelling_novel
                        || (novel[index].id == slot.id
                            && novel[index].weight.to_bits() == slot.weight.to_bits()))
                        && existing.alignments.len() == alns.len()
                        && existing
                            .alignments
                            .iter()
                            .zip(alns)
                            .all(|(left, right)| left.ref_id == right.ref_id)
                        && weights[existing.weight_start..existing.weight_start + alns.len()]
                            .iter()
                            .zip(local_weights)
                            .all(|(left, right)| left.to_bits() == right.to_bits())
                })
            });
            if let Some(index) = duplicate {
                weights.truncate(weight_start);
                classes[index].multiplicity += 1;
                continue;
            }
            by_hash.entry(hash).or_default().push(classes.len());
        }
        classes.push(PreparedEq {
            alignments: alns,
            weight_start,
            multiplicity: 1,
        });
        if modelling_novel {
            novel.push(slot);
        }
    }
    debug_assert!(!modelling_novel || novel.len() == classes.len());
    (classes, weights, novel)
}

/// `NOVEL` is a const generic rather than a runtime flag so that the
/// unannotated-isoform path costs nothing when it is off: with `NOVEL == false`
/// every `novel_id`/`novel_weight` test below folds away at compile time and the
/// loop is byte-for-byte the pre-feature one, including the unconditional
/// single-candidate fast path.
#[inline]
fn m_step_prepared<const NOVEL: bool>(
    eq_iterates: &[PreparedEq],
    weights: &[f64],
    novel: &[NovelSlot],
    prev_count: &[f64],
    curr_counts: &mut [f64],
) {
    for (index, eq) in eq_iterates.iter().enumerate() {
        let read_count = eq.multiplicity as f64;
        let alns = eq.alignments;
        // Short-circuits before the load when NOVEL is false, so the disabled
        // path neither touches `novel` nor keeps the index alive.
        let slot = if NOVEL { novel[index] } else { NO_NOVEL };
        let has_novel = NOVEL && slot.id != usize::MAX;
        if !has_novel && let [alignment] = alns {
            // With no novel alternative the single-candidate weight cancels.
            curr_counts[alignment.ref_id as usize] += read_count;
            continue;
        }
        let eq_weights = &weights[eq.weight_start..eq.weight_start + alns.len()];
        let mut denom: f64 = alns
            .iter()
            .zip(eq_weights)
            .map(|(a, weight)| prev_count[a.ref_id as usize] * weight)
            .sum();
        if has_novel {
            // The novel state competes for this read, so the annotated
            // candidates now receive less than one full read between them.
            denom += prev_count[slot.id] * slot.weight;
        }
        if denom > constants::EM_DENOM_THRESH {
            let scale = read_count / denom;
            for (a, weight) in alns.iter().zip(eq_weights) {
                let target_id = a.ref_id as usize;
                curr_counts[target_id] += prev_count[target_id] * weight * scale;
            }
            if has_novel {
                curr_counts[slot.id] += prev_count[slot.id] * slot.weight * scale;
            }
        }
    }
}

#[inline]
fn m_step_prepared_counts(
    eq_iterates: &[PreparedEq],
    weights: &[f64],
    multiplicities: &[u32],
    prev_count: &[f64],
    curr_counts: &mut [f64],
) {
    for (eq, &multiplicity) in eq_iterates.iter().zip(multiplicities) {
        if multiplicity == 0 {
            continue;
        }
        let read_count = multiplicity as f64;
        let alns = eq.alignments;
        if let [alignment] = alns {
            curr_counts[alignment.ref_id as usize] += read_count;
            continue;
        }
        let eq_weights = &weights[eq.weight_start..eq.weight_start + alns.len()];
        let denom: f64 = alns
            .iter()
            .zip(eq_weights)
            .map(|(a, weight)| prev_count[a.ref_id as usize] * weight)
            .sum();
        if denom > constants::EM_DENOM_THRESH {
            let scale = read_count / denom;
            for (a, weight) in alns.iter().zip(eq_weights) {
                let target_id = a.ref_id as usize;
                curr_counts[target_id] += prev_count[target_id] * weight * scale;
            }
        }
    }
}

/// Performs one iteration of the EM algorithm by looping over all
/// alignments and computing their estimated probability of being
/// the true alignment (using the abunance estimates from `prev_counts`).
/// Then, `curr_counts` is computed by summing over the expected assignment
/// likelihood for all reads mapping to each target.
/// As with [`m_step_prepared`], `NOVEL` is a const generic so the disabled path
/// is exactly the original loop. Shards are sized by the caller to cover the
/// novel states, which live past the annotated transcripts, so a novel id
/// indexes `local` directly and the reduction below picks it up unchanged.
#[inline]
fn m_step_par<const NOVEL: bool>(
    eq_iterates: &[PreparedEq],
    weights: &[f64],
    novel: &[NovelSlot],
    prev_count: &[f64],
    curr_counts: &mut [f64],
    shards: &mut [Vec<f64>],
) {
    let chunk_size = eq_iterates.len().div_ceil(shards.len().max(1));
    shards
        .par_iter_mut()
        .enumerate()
        .for_each(|(shard, local)| {
            local.fill(0.0);
            let start = shard * chunk_size;
            let end = ((shard + 1) * chunk_size).min(eq_iterates.len());
            for (offset, eq) in eq_iterates[start..end].iter().enumerate() {
                let alns = eq.alignments;
                let slot = if NOVEL {
                    novel[start + offset]
                } else {
                    NO_NOVEL
                };
                let has_novel = NOVEL && slot.id != usize::MAX;
                // A read with one candidate contributes exactly one count: its
                // alignment/coverage weight cancels between numerator and
                // denominator. Avoid evaluating that weight twice. With a novel
                // competitor present the weight no longer cancels.
                if !has_novel && let [alignment] = alns {
                    local[alignment.ref_id as usize] += eq.multiplicity as f64;
                    continue;
                }
                let eq_weights = &weights[eq.weight_start..eq.weight_start + alns.len()];
                let mut denom = 0.0_f64;
                for (a, weight) in alns.iter().zip(eq_weights) {
                    let target_id = a.ref_id as usize;
                    denom += prev_count[target_id] * weight;
                }
                if has_novel {
                    denom += prev_count[slot.id] * slot.weight;
                }
                if denom > constants::EM_DENOM_THRESH {
                    let scale = eq.multiplicity as f64 / denom;
                    for (a, weight) in alns.iter().zip(eq_weights) {
                        let target_id = a.ref_id as usize;
                        local[target_id] += prev_count[target_id] * weight * scale;
                    }
                    if has_novel {
                        local[slot.id] += prev_count[slot.id] * slot.weight * scale;
                    }
                }
            }
        });
    curr_counts
        .par_iter_mut()
        .enumerate()
        .for_each(|(target_id, count)| {
            *count = shards.iter().map(|local| local[target_id]).sum();
        });
}

#[derive(Debug, Clone)]
pub struct EMResult {
    pub counts: Vec<f64>,
    pub evaluations: u32,
    pub converged: bool,
    /// Per-transcript presence probability `q_t` when `--presence-model` is
    /// active; `None` otherwise. Indexed like `counts` (novel states pinned
    /// at 1.0).
    pub presence: Option<Vec<f64>>,
}

/// Spike-and-slab presence machinery (`--presence-model spike-slab`).
///
/// Each transcript carries a Bernoulli presence posterior `q_t`. The M-step
/// sees the *effective* abundance `q_t * counts[t]`, which is the standard
/// zero-inflated EM update: a transcript with low presence probability
/// competes for shared reads at a discount and its ambiguous mass drains to
/// its supported competitors, while reads unique to a transcript always keep
/// it (their evidence drives `q_t -> 1`). `q` is refreshed every
/// `period` evaluations (after `warmup`) from a leave-one-out evidence pass:
///
///   delta_t = sum over reads touching t of
///               m_r * ln( (rest_r + slab_t·w) / rest_r ),
///   q_t     = sigmoid( delta_t + logit(rho) ),
///
/// where `rest_r` is the read's denominator without t and `slab_t =
/// counts[t]/max(q_t, 1e-3)` (the abundance conditional on presence).
///
/// `delta_t` is a leave-one-out log-likelihood *gain* and is therefore always
/// >= 0 (adding a mixture component never lowers a read's likelihood), so the
/// prior `rho` must act as a fixed penalty: an empirical-Bayes update of the
/// form `rho = mean(q)` only ratchets upward and suppresses nothing. `rho`
/// is a fixed prior (`--presence-rho`, default 0.05, logit ~= -2.94). The
/// units are interpretable: for a transcript holding small posterior shares,
/// `delta_t ~= sum of its read shares ~= its assigned read count`, so the
/// default demands roughly three effective exclusive reads for q > 0.5 —
/// smoothly, and self-reinforcing through the q-discounted M-step in both
/// directions (a suppressed transcript's shares shrink, a supported one's
/// grow).
const POS_BINS: usize = 40;
const POS_STRATA: [f64; 3] = [750.0, 1250.0, 2000.0];

/// Length-stratified endpoint-gap model fitted from unique reads: per stratum,
/// binned distributions of the relative left/right alignment gaps
/// (start/len, (len-end)/len). Stored as `ln(p_bin * POS_BINS)` so each value
/// is a log-likelihood ratio against a uniform positional reference: positive
/// = the endpoint sits where this library's reads typically sit, negative =
/// atypical placement. Measured (F4.1, zone-sep): the per-transcript
/// responsibility-weighted sum of these terms separates truly-expressed
/// overlap-zone transcripts from false positives at AUC 0.76-0.85 across all
/// six Panel B samples.
struct PositionalModel {
    log_lr: [[[f64; POS_BINS]; 2]; 4],
}

impl PositionalModel {
    fn stratum(tlen: f64) -> usize {
        POS_STRATA.iter().filter(|b| tlen >= **b).count()
    }

    fn fit(em_info: &EMInfo) -> Self {
        let mut hist = [[[1.0f64; POS_BINS]; 2]; 4]; // +1 smoothing
        for (alns, _, _) in em_info.eq_map.iter() {
            if alns.len() != 1 {
                continue;
            }
            let a = &alns[0];
            let tlen = em_info.txp_info[a.ref_id as usize].lenf;
            if tlen <= 0.0 {
                continue;
            }
            let st = Self::stratum(tlen);
            let lg = ((a.start.saturating_sub(1)) as f64 / tlen).clamp(0.0, 1.0);
            let rg = ((tlen - a.end as f64).max(0.0) / tlen).clamp(0.0, 1.0);
            hist[st][0][((lg * POS_BINS as f64) as usize).min(POS_BINS - 1)] += 1.0;
            hist[st][1][((rg * POS_BINS as f64) as usize).min(POS_BINS - 1)] += 1.0;
        }
        let mut log_lr = [[[0.0f64; POS_BINS]; 2]; 4];
        for st in 0..4 {
            for side in 0..2 {
                let tot: f64 = hist[st][side].iter().sum();
                for bin in 0..POS_BINS {
                    log_lr[st][side][bin] = (hist[st][side][bin] / tot * POS_BINS as f64).ln();
                }
            }
        }
        Self { log_lr }
    }

    #[inline]
    fn read_log_lr(&self, aln: &AlnInfo, tlen: f64) -> f64 {
        if tlen <= 0.0 {
            return 0.0;
        }
        let st = Self::stratum(tlen);
        let lg = ((aln.start.saturating_sub(1)) as f64 / tlen).clamp(0.0, 1.0);
        let rg = ((tlen - aln.end as f64).max(0.0) / tlen).clamp(0.0, 1.0);
        self.log_lr[st][0][((lg * POS_BINS as f64) as usize).min(POS_BINS - 1)]
            + self.log_lr[st][1][((rg * POS_BINS as f64) as usize).min(POS_BINS - 1)]
    }
}

struct PresenceState {
    q: Vec<f64>,
    eff: Vec<f64>,
    delta: Vec<f64>,
    /// Responsibility-weighted positional log-LR accumulator (endpoint-gap
    /// plausibility evidence), scaled by `endpoint_alpha` into `delta`.
    pos: Vec<f64>,
    positional: Option<PositionalModel>,
    endpoint_alpha: f64,
    touched: Vec<bool>,
    rho: f64,
    calls: u32,
    warmup: u32,
    period: u32,
    n_txps: usize,
}

impl PresenceState {
    fn new(
        n_txps: usize,
        n_states: usize,
        warmup: u32,
        period: u32,
        rho: f64,
        positional: Option<PositionalModel>,
        endpoint_alpha: f64,
    ) -> Self {
        Self {
            q: vec![1.0; n_states],
            eff: vec![0.0; n_states],
            delta: vec![0.0; n_txps],
            pos: vec![0.0; n_txps],
            positional,
            endpoint_alpha,
            touched: vec![false; n_txps],
            rho,
            calls: 0,
            warmup,
            period: period.max(1),
            n_txps,
        }
    }

    /// Per-read clamp on the log-evidence term (a single unique read saturates
    /// here) and overall clamp on delta before the sigmoid.
    const READ_LOG_CAP: f64 = 50.0;
    const DELTA_CAP: f64 = 40.0;
    const EPS: f64 = 1e-100;

    fn update_q(
        &mut self,
        classes: &[PreparedEq],
        weights: &[f64],
        novel: &[NovelSlot],
        counts: &[f64],
        txp_lens: &[f64],
    ) {
        self.delta.fill(0.0);
        self.pos.fill(0.0);
        let modelling_novel = !novel.is_empty();
        for (index, class) in classes.iter().enumerate() {
            let m = class.multiplicity as f64;
            let ws = &weights[class.weight_start..class.weight_start + class.alignments.len()];
            let mut denom = 0.0;
            for (a, w) in class.alignments.iter().zip(ws) {
                let t = a.ref_id as usize;
                denom += self.q[t] * counts[t] * w;
            }
            if modelling_novel {
                let slot = novel[index];
                if slot.id != usize::MAX {
                    denom += counts[slot.id] * slot.weight;
                }
            }
            // Positional plausibility: responsibility-weighted endpoint-gap
            // log-LR against a uniform reference (see PositionalModel).
            if denom > 0.0
                && let Some(pm) = &self.positional
            {
                for (a, w) in class.alignments.iter().zip(ws) {
                    let t = a.ref_id as usize;
                    let gam = self.q[t] * counts[t] * w / denom * m;
                    if gam > 1e-9 {
                        self.pos[t] += gam * pm.read_log_lr(a, txp_lens[t]);
                    }
                }
            }
            // Per-tid leave-one-out: contribution of t at slab strength vs
            // its current effective contribution. Candidate lists are short,
            // so the quadratic-in-duplicates handling below (skip repeated
            // tids) costs nothing in practice.
            for (i, (a, w)) in class.alignments.iter().zip(ws).enumerate() {
                let t = a.ref_id as usize;
                if class.alignments[..i].iter().any(|b| b.ref_id == a.ref_id) {
                    continue; // this tid already handled for this class
                }
                self.touched[t] = true;
                // summed weight of t's alignments in this class
                let mut w_sum = *w;
                for (b, wb) in class.alignments.iter().zip(ws).skip(i + 1) {
                    if b.ref_id == a.ref_id {
                        w_sum += wb;
                    }
                }
                let eff_contrib = self.q[t] * counts[t] * w_sum;
                let rest = (denom - eff_contrib).max(0.0);
                let slab = counts[t] / self.q[t].max(1e-3);
                let term = ((rest + slab * w_sum + Self::EPS) / (rest + Self::EPS))
                    .ln()
                    .min(Self::READ_LOG_CAP);
                self.delta[t] += m * term;
            }
        }
        let logit_rho = (self.rho / (1.0 - self.rho)).ln();
        for t in 0..self.n_txps {
            if !self.touched[t] {
                continue;
            }
            // Leave-one-out gain is >= 0; the positional term is signed
            // (typical placements add evidence, atypical placements subtract).
            let d = (self.delta[t].max(0.0) + self.endpoint_alpha * self.pos[t])
                .clamp(-Self::DELTA_CAP, Self::DELTA_CAP);
            self.q[t] = 1.0 / (1.0 + (-(d + logit_rho)).exp());
        }
    }
}

/// Run the EM driver, optionally wrapping the fixed point with the
/// spike-and-slab presence model. `classes`/`weights`/`novel` are the same
/// prepared structures the fixed point iterates; they are needed again for
/// the leave-one-out evidence pass.
fn drive_with_presence(
    em_info: &EMInfo,
    classes: &[PreparedEq],
    weights: &[f64],
    novel: &[NovelSlot],
    mut inner: impl FnMut(&[f64], &mut [f64]),
    do_log: bool,
    allow_presence: bool,
) -> EMResult {
    use crate::prog_opts::PresenceModel;
    if !allow_presence || em_info.presence_model == PresenceModel::None {
        return run_driver(em_info, &mut inner, do_log);
    }
    let n_txps = em_info.txp_info.len();
    let n_states = n_txps + em_info.novel_loci;
    let endpoint_alpha = em_info.presence_endpoint_alpha;
    let positional = (endpoint_alpha != 0.0).then(|| PositionalModel::fit(em_info));
    let txp_lens: Vec<f64> = em_info.txp_info.iter().map(|t| t.lenf).collect();
    let mut st = PresenceState::new(
        n_txps,
        n_states,
        em_info.presence_warmup,
        em_info.presence_period,
        em_info.presence_rho,
        positional,
        endpoint_alpha,
    );
    let mut wrapped = |src: &[f64], dst: &mut [f64]| {
        st.calls += 1;
        if st.calls > st.warmup && (st.calls - st.warmup).is_multiple_of(st.period) {
            st.update_q(classes, weights, novel, src, &txp_lens);
        }
        for i in 0..src.len() {
            st.eff[i] = st.q[i] * src[i];
        }
        inner(&st.eff, dst);
    };
    let mut result = run_driver(em_info, &mut wrapped, do_log);
    result.presence = Some(std::mem::take(&mut st.q));
    result
}

fn convergence_distance(previous: &[f64], current: &[f64]) -> f64 {
    // Symmetric denominator: dividing by the *current* value alone makes a
    // geometrically decaying parameter (c = r·p) report a constant relative
    // change (1-r)/r no matter how small it gets, so a single dying transcript
    // can hold the whole EM above threshold until it crosses the active-count
    // gate. max(p, c) bounds the ratio by 1 and decays with the parameter.
    previous
        .iter()
        .zip(current)
        .filter(|(p, c)| (**p).max(**c) >= constants::MIN_ACTIVE_COUNT)
        .map(|(p, c)| (c - p).abs() / p.max(*c).max(constants::EM_DENOM_THRESH))
        .fold(0.0, f64::max)
}

/// Global relative L1 change between successive iterates, `Σ|Δc| / Σc`.
/// Logged as a mass-weighted companion to the max-relative criterion.
fn relative_l1(previous: &[f64], current: &[f64]) -> f64 {
    let (num, den) = previous
        .iter()
        .zip(current)
        .fold((0.0_f64, 0.0_f64), |(num, den), (p, c)| {
            (num + (c - p).abs(), den + c)
        });
    num / den.max(constants::EM_DENOM_THRESH)
}

fn run_driver(
    em_info: &EMInfo,
    mut fixed_point: impl FnMut(&[f64], &mut [f64]),
    do_log: bool,
) -> EMResult {
    const MIN_EVAL: u32 = 50;
    // Novel latent states live past the annotated transcripts in the same
    // vector, so downstream consumers indexing by reference id are unaffected.
    let n = em_info.txp_info.len() + em_info.novel_loci;
    let avg = if n == 0 {
        0.0
    } else {
        em_info.eq_map.num_aligned_reads() as f64 / n as f64
    };
    let has_initial_abundances = em_info
        .init_abundances
        .as_ref()
        .is_some_and(|v| v.len() == n);
    let mut counts = em_info
        .init_abundances
        .as_ref()
        .filter(|v| v.len() == n)
        .cloned()
        .unwrap_or_else(|| vec![avg; n]);
    // Exact zeros are absorbing states for multiplicative EM updates.  A
    // warm-start estimate may contain zeros even though the new likelihood
    // supports those transcripts, so retain a tiny positive route back into
    // the active set.  This also makes accelerator choice independent of the
    // warm-start's zero pattern.
    if has_initial_abundances {
        counts
            .iter_mut()
            .for_each(|x| *x = x.max(10.0 * constants::MIN_READ_THRESH));
    }
    let mut next = vec![0.0; n];
    let mut evaluations = 0;
    // Reserve one evaluation for the historical post-threshold redistribution
    // step so --max-em-iter is a strict total M-step budget.
    let iteration_budget = em_info.max_iter.saturating_sub(1);

    // Combined stopping rule. The per-transcript L-infinity criterion demands
    // that every active transcript's relative change fall below the threshold;
    // at ~200k transcripts there is essentially always some low-count
    // transcript in slow geometric decay (~0.1-1%/iteration), so that bar is
    // unreachable within any practical budget even under acceleration, while
    // the aggregate estimate is long since stable. The mass-weighted
    // relative-L1 (sum |delta| / sum counts) captures the latter: when it
    // falls below `convergence_l1_thresh`, total count movement per iteration
    // is a negligible fraction of the library and we declare convergence.
    let l1_thresh = em_info.convergence_l1_thresh;
    let distance = |p: &[f64], c: &[f64]| -> f64 {
        if l1_thresh > 0.0 && relative_l1(p, c) < l1_thresh {
            return 0.0;
        }
        convergence_distance(p, c)
    };

    let converged = match em_info.accel {
        crate::prog_opts::EmAccel::None => {
            let mut done = false;
            while evaluations < iteration_budget {
                fixed_point(&counts, &mut next);
                evaluations += 1;
                let d = distance(&counts, &next);
                std::mem::swap(&mut counts, &mut next);
                if do_log && evaluations.is_multiple_of(10) {
                    if evaluations.is_multiple_of(100) {
                        info!(
                            "iteration {}; rel diff {}; rel l1 {}",
                            evaluations.to_formatted_string(&Locale::en),
                            d,
                            relative_l1(&next, &counts)
                        );
                    } else {
                        trace!(
                            "iteration {}; rel diff {}",
                            evaluations.to_formatted_string(&Locale::en),
                            d
                        );
                    }
                }
                if evaluations >= MIN_EVAL.min(iteration_budget)
                    && d.is_finite()
                    && d < em_info.convergence_thresh
                {
                    done = true;
                    break;
                }
            }
            done
        }
        crate::prog_opts::EmAccel::Squarem => crate::em_accel::squarem(
            &mut fixed_point,
            &mut counts,
            &mut next,
            iteration_budget,
            MIN_EVAL.min(iteration_budget),
            &mut evaluations,
            em_info.convergence_thresh,
            &distance,
        ),
        crate::prog_opts::EmAccel::Daarem => crate::em_accel::daarem(
            &mut fixed_point,
            &mut counts,
            &mut next,
            iteration_budget,
            MIN_EVAL.min(iteration_budget),
            &mut evaluations,
            em_info.convergence_thresh,
            &distance,
        ),
    };

    // Zero out estimates below the count floor, then let one final fixed-point
    // evaluation redistribute their mass.  The historical floor (1e-5) only
    // removes numerical dust; a higher floor (e.g. 0.5 reads) additionally
    // suppresses low-mass false-positive transcripts at the cost of the ~35%
    // of such transcripts that are real.
    let count_floor = em_info.count_floor.max(constants::MIN_READ_THRESH);
    counts.iter_mut().for_each(|x| {
        if *x < count_floor {
            *x = 0.0;
        }
    });
    if evaluations < em_info.max_iter {
        fixed_point(&counts, &mut next);
        evaluations += 1;
    } else {
        next.copy_from_slice(&counts);
    }
    if do_log {
        info!(evaluations, converged, accelerator = ?em_info.accel, "EM completed");
    }
    EMResult {
        counts: next,
        evaluations,
        converged,
        presence: None,
    }
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
#[allow(dead_code)]
pub fn em(em_info: &EMInfo, _nthreads: usize) -> EMResult {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    let (eq_iterates, weights, novel) = prepare_equivalence_classes(em_info, true);
    // Dispatch once, outside the fixed point, so the hot loop is monomorphic.
    if em_info.novel_loci > 0 {
        let n_txps = em_info.txp_info.len();
        let fixed_point = |src: &[f64], dst: &mut [f64]| {
            dst.fill(0.0);
            m_step_prepared::<true>(&eq_iterates, &weights, &novel, src, dst);
            // Projection-failed reads are deterministic single-candidate
            // reads of their locus's novel state (--novel-from-failures).
            for (locus, mass) in em_info.novel_base_mass.iter().enumerate() {
                dst[n_txps + locus] += mass;
            }
        };
        drive_with_presence(em_info, &eq_iterates, &weights, &novel, fixed_point, true, true)
    } else {
        let fixed_point = |src: &[f64], dst: &mut [f64]| {
            dst.fill(0.0);
            m_step_prepared::<false>(&eq_iterates, &weights, &novel, src, dst);
        };
        drive_with_presence(em_info, &eq_iterates, &weights, &novel, fixed_point, true, true)
    }
}

fn do_bootstrap_prepared(
    em_info: &EMInfo,
    eq_iterates: &[PreparedEq],
    weights: &[f64],
) -> Vec<f64> {
    let mut rng = trng();
    let n = em_info.eq_map.len();
    let inds = bootstrap::get_sample_inds(n, &mut rng);
    let mut multiplicities = vec![0_u32; n];
    inds.iter().for_each(|&index| multiplicities[index] += 1);
    let mut fixed_point = |src: &[f64], dst: &mut [f64]| {
        dst.fill(0.0);
        m_step_prepared_counts(eq_iterates, weights, &multiplicities, src, dst);
    };
    run_driver(em_info, &mut fixed_point, false).counts
}

pub fn bootstrap(em_info: &EMInfo, num_boot: u32, nthreads: usize) -> Vec<Vec<f64>> {
    let span = span!(tracing::Level::INFO, "bootstrap");
    let _guard = span.enter();

    info!("will collection {num_boot} bootstraps");

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    let (eq_iterates, weights, _novel) = prepare_equivalence_classes(em_info, false);

    pool.install(|| {
        (0..num_boot)
            .into_par_iter()
            .map(|i| {
                let span = span!(tracing::Level::INFO, "bootstrap");
                let _guard = span.enter();
                info!("evaluating bootstrap replicate {}", i);
                do_bootstrap_prepared(em_info, &eq_iterates, &weights)
            })
            .collect()
    })
}

/// Perform the EM algorithm to estimate the abundances of the
/// target sequences.  The return value is a `Vec` of f64 values,
/// each of which is the estimated number of fragments arising from
/// each target.
pub fn em_par(em_info: &EMInfo, nthreads: usize) -> EMResult {
    let span = span!(tracing::Level::INFO, "em");
    let _guard = span.enter();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    // Pack invariant alignment, coverage and density terms once. They are
    // otherwise recomputed twice per alignment on every EM iteration.
    let (eq_iterates, weights, novel) = prepare_equivalence_classes(em_info, false);
    // Reuse private dense accumulators across all fixed-point evaluations.
    // Capping the shard count avoids excessive memory use on high-core hosts.
    // Shards must span the novel states as well, since those ids index past the
    // annotated transcripts into the same accumulator.
    let shard_count = nthreads.clamp(1, 64);
    let n_states = em_info.txp_info.len() + em_info.novel_loci;
    let mut shards = vec![vec![0.0; n_states]; shard_count];
    // Dispatch once, outside the fixed point, so the hot loop is monomorphic.
    if em_info.novel_loci > 0 {
        let n_txps = em_info.txp_info.len();
        let fixed_point = |src: &[f64], dst: &mut [f64]| {
            m_step_par::<true>(&eq_iterates, &weights, &novel, src, dst, &mut shards);
            // Projection-failed reads are deterministic single-candidate
            // reads of their locus's novel state (--novel-from-failures).
            for (locus, mass) in em_info.novel_base_mass.iter().enumerate() {
                dst[n_txps + locus] += mass;
            }
        };
        pool.install(|| {
            drive_with_presence(em_info, &eq_iterates, &weights, &novel, fixed_point, true, true)
        })
    } else {
        let fixed_point = |src: &[f64], dst: &mut [f64]| {
            m_step_par::<false>(&eq_iterates, &weights, &novel, src, dst, &mut shards);
        };
        pool.install(|| {
            drive_with_presence(em_info, &eq_iterates, &weights, &novel, fixed_point, true, true)
        })
    }
}

#[cfg(test)]
mod size_tests {
    use super::{NovelSlot, PreparedEq};
    /// The novel-state data must stay out of the hot-loop struct: carrying it
    /// inline costs 16 bytes per equivalence class on every run, enabled or not.
    #[test]
    fn prepared_eq_does_not_carry_novel_state() {
        assert_eq!(
            std::mem::size_of::<PreparedEq>(),
            std::mem::size_of::<&[crate::util::oarfish_types::AlnInfo]>()
                + std::mem::size_of::<usize>()
                + std::mem::size_of::<u32>()
                + 4, // padding
        );
        assert_eq!(std::mem::size_of::<NovelSlot>(), 16);
    }
}
