//! Empirical joint endpoint likelihood for bulk long-read alignments.
//!
//! The model is deliberately trained only on reads with one retained candidate,
//! preventing ambiguous reads from manufacturing their own coverage support.

use crate::util::oarfish_types::{AlnInfo, InMemoryAlignmentStore, TranscriptInfo};
use bio_types::strand::Strand;
use tracing::info;

const GRID: usize = 20;
const CELLS: usize = GRID * GRID;
const PRIOR_MASS: f64 = 1_000.0;
const STRATA: usize = 4;

#[derive(Debug, Clone)]
struct EndpointModel {
    length_bounds: [usize; 3],
    counts: [[f64; CELLS]; STRATA],
    totals: [f64; STRATA],
}

pub(crate) struct EndpointProbabilities {
    pub probabilities: Vec<f64>,
    /// One gate per read, rather than one duplicated copy per alignment.
    pub support_gates: Vec<f32>,
    pub training_reads: usize,
}

pub(crate) struct AdaptiveEndpointProbabilities {
    pub probabilities: Vec<f64>,
    /// One gate per read, rather than one duplicated copy per alignment.
    pub support_gates: Vec<f32>,
    pub training_reads: usize,
    pub selected_prior_mass: f64,
    pub folds: usize,
    /// One gate per read; empty unless uncertainty gating is requested.
    pub uncertainty_gates: Vec<f32>,
    pub ambiguous_training_reads: usize,
}

impl EndpointModel {
    fn new(txps: &[TranscriptInfo]) -> Self {
        let mut lengths: Vec<usize> = txps.iter().map(|t| t.len.get()).collect();
        lengths.sort_unstable();
        let at = |q: usize| {
            lengths
                .get(q * lengths.len().saturating_sub(1) / 4)
                .copied()
                .unwrap_or(1)
        };
        Self {
            length_bounds: [at(1), at(2), at(3)],
            counts: [[0.0; CELLS]; STRATA],
            totals: [0.0; STRATA],
        }
    }

    fn stratum(&self, len: usize) -> usize {
        self.length_bounds.iter().take_while(|b| len > **b).count()
    }

    fn cell(aln: &AlnInfo, len: usize) -> usize {
        let (five_gap, three_gap) = endpoint_gaps(aln, len);
        let x = (five_gap * GRID as f64) as usize;
        let y = (three_gap * GRID as f64) as usize;
        x * GRID + y
    }
}

pub(crate) fn endpoint_gaps(aln: &AlnInfo, len: usize) -> (f64, f64) {
    let lenf = len.max(1) as f64;
    let left_gap = aln.start.saturating_sub(1) as f64 / lenf;
    let right_gap = len.saturating_sub(aln.end as usize) as f64 / lenf;
    let (five_gap, three_gap) = if aln.strand == Strand::Reverse {
        (right_gap, left_gap)
    } else {
        (left_gap, right_gap)
    };
    (
        five_gap.clamp(0.0, 1.0 - f64::EPSILON),
        three_gap.clamp(0.0, 1.0 - f64::EPSILON),
    )
}

impl EndpointModel {
    fn observe(&mut self, aln: &AlnInfo, len: usize) {
        let s = self.stratum(len);
        self.counts[s][Self::cell(aln, len)] += 1.0;
        self.totals[s] += 1.0;
    }

    fn probability(&self, aln: &AlnInfo, len: usize) -> f64 {
        let s = self.stratum(len);
        let prior = PRIOR_MASS / CELLS as f64;
        (self.counts[s][Self::cell(aln, len)] + prior) / (self.totals[s] + PRIOR_MASS)
    }

    fn support(&self, aln: &AlnInfo, len: usize) -> f64 {
        let s = self.stratum(len);
        self.counts[s][Self::cell(aln, len)]
    }
}

pub(crate) fn endpoint_probabilities(
    store: &InMemoryAlignmentStore,
    txps: &[TranscriptInfo],
    support_scale: f64,
) -> EndpointProbabilities {
    let mut model = EndpointModel::new(txps);
    let mut training_reads = 0usize;
    for (alns, _, _) in store.iter() {
        if let [aln] = alns {
            let len = txps[aln.ref_id as usize].len.get();
            model.observe(aln, len);
            training_reads += 1;
        }
    }

    let mut probabilities = Vec::with_capacity(store.total_len());
    let mut support_gates = Vec::with_capacity(store.len());
    for (alns, _, _) in store.iter() {
        let mut local: Vec<f64> = alns
            .iter()
            .map(|aln| model.probability(aln, txps[aln.ref_id as usize].len.get()))
            .collect();
        let sum: f64 = local.iter().sum();
        if sum.is_finite() && sum > 0.0 {
            local.iter_mut().for_each(|p| *p /= sum);
        } else if !local.is_empty() {
            let uniform = 1.0 / local.len() as f64;
            local.fill(uniform);
        }
        let mean_support = if alns.is_empty() {
            0.0
        } else {
            alns.iter()
                .map(|aln| model.support(aln, txps[aln.ref_id as usize].len.get()))
                .sum::<f64>()
                / alns.len() as f64
        };
        let gate = mean_support / (mean_support + support_scale);
        support_gates.push(gate as f32);
        probabilities.extend(local);
    }
    EndpointProbabilities {
        probabilities,
        support_gates,
        training_reads,
    }
}

pub fn apply_endpoint_probabilities(store: &mut InMemoryAlignmentStore, txps: &[TranscriptInfo]) {
    let result = endpoint_probabilities(store, txps, PRIOR_MASS / CELLS as f64);
    store.coverage_probabilities = result.probabilities;
    let training_reads = result.training_reads;
    info!(training_reads, "learned joint endpoint coverage model");
}

fn smooth_counts(counts: &[[f64; CELLS]; STRATA]) -> ([[f64; CELLS]; STRATA], [f64; STRATA]) {
    let mut smoothed = [[0.0; CELLS]; STRATA];
    let mut totals = [0.0; STRATA];
    for stratum in 0..STRATA {
        for x in 0..GRID {
            for y in 0..GRID {
                let mut weighted = 0.0;
                let mut weight_sum = 0.0;
                for nx in x.saturating_sub(1)..=(x + 1).min(GRID - 1) {
                    for ny in y.saturating_sub(1)..=(y + 1).min(GRID - 1) {
                        let weight =
                            if nx == x { 2.0 } else { 1.0 } * if ny == y { 2.0 } else { 1.0 };
                        weighted += weight * counts[stratum][nx * GRID + ny];
                        weight_sum += weight;
                    }
                }
                let cell = x * GRID + y;
                smoothed[stratum][cell] = weighted / weight_sum;
                totals[stratum] += smoothed[stratum][cell];
            }
        }
    }
    (smoothed, totals)
}

/// Cross-fitted endpoint probabilities with the symmetric Dirichlet prior mass
/// selected by held-out predictive likelihood of unambiguous reads.
pub(crate) fn adaptive_endpoint_probabilities(
    store: &InMemoryAlignmentStore,
    txps: &[TranscriptInfo],
    folds: usize,
    support_scale: f64,
    ablation: crate::prog_opts::CoverageAblation,
) -> AdaptiveEndpointProbabilities {
    const PRIOR_GRID: [f64; 7] = [10.0, 30.0, 100.0, 300.0, 1_000.0, 3_000.0, 10_000.0];
    let folds = folds.max(2);
    let template = EndpointModel::new(txps);
    let mut fold_counts = vec![[[0.0; CELLS]; STRATA]; folds];
    let mut total_counts = [[0.0; CELLS]; STRATA];
    let mut observations = Vec::new();
    let use_eq_classes = matches!(
        ablation,
        crate::prog_opts::CoverageAblation::EqClassTraining
            | crate::prog_opts::CoverageAblation::AllCandidates
    );
    let mut ambiguous_training_reads = 0usize;
    for (read_index, (alignments, score_probs, _)) in store.iter().enumerate() {
        if let [alignment] = alignments {
            let len = txps[alignment.ref_id as usize].len.get();
            let stratum = template.stratum(len);
            let cell = EndpointModel::cell(alignment, len);
            let fold = read_index % folds;
            fold_counts[fold][stratum][cell] += 1.0;
            total_counts[stratum][cell] += 1.0;
            observations.push((fold, stratum, cell));
        } else if use_eq_classes && !alignments.is_empty() {
            let total: f64 = score_probs.iter().map(|value| *value as f64).sum();
            if total > 0.0 {
                ambiguous_training_reads += 1;
                let fold = read_index % folds;
                for (alignment, score) in alignments.iter().zip(score_probs) {
                    let len = txps[alignment.ref_id as usize].len.get();
                    let stratum = template.stratum(len);
                    let cell = EndpointModel::cell(alignment, len);
                    let weight = *score as f64 / total;
                    fold_counts[fold][stratum][cell] += weight;
                    total_counts[stratum][cell] += weight;
                }
            }
        }
    }

    let training_for_fold = |fold: usize| {
        let mut counts = total_counts;
        for stratum in 0..STRATA {
            for cell in 0..CELLS {
                counts[stratum][cell] -= fold_counts[fold][stratum][cell];
            }
        }
        smooth_counts(&counts)
    };
    let fitted: Vec<_> = (0..folds).map(training_for_fold).collect();
    let selected_prior_mass = PRIOR_GRID
        .iter()
        .copied()
        .max_by(|left, right| {
            let score = |prior: f64| {
                observations
                    .iter()
                    .map(|&(fold, stratum, cell)| {
                        let (counts, totals) = &fitted[fold];
                        ((counts[stratum][cell] + prior / CELLS as f64) / (totals[stratum] + prior))
                            .max(f64::MIN_POSITIVE)
                            .ln()
                    })
                    .sum::<f64>()
            };
            score(*left).total_cmp(&score(*right))
        })
        .unwrap_or(PRIOR_MASS);

    // Held-out observations are needed only for selecting the prior. Releasing
    // them before allocating per-alignment output avoids overlapping two large
    // buffers on deep libraries.
    let training_reads = observations.len();
    drop(observations);
    drop(fold_counts);

    let mut probabilities = Vec::with_capacity(store.total_len());
    let mut support_gates = Vec::with_capacity(store.len());
    let use_uncertainty = matches!(
        ablation,
        crate::prog_opts::CoverageAblation::UncertaintyGate
            | crate::prog_opts::CoverageAblation::AllCandidates
    );
    let mut uncertainty_gates = if use_uncertainty {
        Vec::with_capacity(store.len())
    } else {
        Vec::new()
    };
    for (read_index, (alignments, _, _)) in store.iter().enumerate() {
        let fold = read_index % folds;
        let (counts, totals) = &fitted[fold];
        let mut local = Vec::with_capacity(alignments.len());
        let mut support = 0.0;
        for alignment in alignments {
            let len = txps[alignment.ref_id as usize].len.get();
            let stratum = template.stratum(len);
            let cell = EndpointModel::cell(alignment, len);
            local.push(
                (counts[stratum][cell] + selected_prior_mass / CELLS as f64)
                    / (totals[stratum] + selected_prior_mass),
            );
            support += counts[stratum][cell];
        }
        let uncertainty_gate = if use_uncertainty {
            let mut fold_disagreement = 0.0;
            for alignment in alignments {
                let len = txps[alignment.ref_id as usize].len.get();
                let stratum = template.stratum(len);
                let cell = EndpointModel::cell(alignment, len);
                let mut sum = 0.0;
                let mut sum_sq = 0.0;
                for (fold_counts, fold_totals) in &fitted {
                    let value = (fold_counts[stratum][cell] + selected_prior_mass / CELLS as f64)
                        / (fold_totals[stratum] + selected_prior_mass);
                    sum += value;
                    sum_sq += value * value;
                }
                let mean = sum / fitted.len() as f64;
                if mean > 0.0 {
                    fold_disagreement +=
                        (sum_sq / fitted.len() as f64 - mean * mean).max(0.0) / mean.powi(2);
                }
            }
            1.0 / (1.0 + fold_disagreement / alignments.len().max(1) as f64)
        } else {
            1.0
        };
        let sum: f64 = local.iter().sum();
        if sum > 0.0 && sum.is_finite() {
            local.iter_mut().for_each(|value| *value /= sum);
        } else if !local.is_empty() {
            local.fill(1.0 / alignments.len() as f64);
        }
        let mean_support = support / alignments.len().max(1) as f64;
        let gate = mean_support / (mean_support + support_scale);
        support_gates.push(gate as f32);
        if use_uncertainty {
            uncertainty_gates.push(uncertainty_gate as f32);
        }
        probabilities.extend(local);
    }

    AdaptiveEndpointProbabilities {
        probabilities,
        support_gates,
        training_reads,
        selected_prior_mass,
        folds,
        uncertainty_gates,
        ambiguous_training_reads,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn aln(start: u32, end: u32) -> AlnInfo {
        AlnInfo {
            ref_id: 0,
            start,
            end,
            prob: 0.0,
            strand: Strand::Forward,
            left_clip: 0,
            right_clip: 0,
        }
    }

    #[test]
    fn endpoint_cells_distinguish_opposite_truncation() {
        assert_ne!(
            EndpointModel::cell(&aln(1, 900), 1000),
            EndpointModel::cell(&aln(101, 1000), 1000)
        );
    }

    #[test]
    fn untrained_model_is_neutral() {
        let model = EndpointModel {
            length_bounds: [500, 1000, 2000],
            counts: [[0.0; CELLS]; STRATA],
            totals: [0.0; STRATA],
        };
        let p1 = model.probability(&aln(1, 900), 1000);
        let p2 = model.probability(&aln(101, 1000), 1000);
        assert!((p1 - p2).abs() < f64::EPSILON);
    }

    #[test]
    fn spatial_smoothing_borrows_from_neighboring_cells() {
        let mut counts = [[0.0; CELLS]; STRATA];
        counts[0][10 * GRID + 10] = 9.0;
        let (smoothed, _) = smooth_counts(&counts);
        assert!(smoothed[0][10 * GRID + 10] > 0.0);
        assert!(smoothed[0][10 * GRID + 9] > 0.0);
        assert_eq!(smoothed[0][0], 0.0);
    }
}
