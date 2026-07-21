//! Regularized log-linear combination of coverage likelihoods.

use crate::util::endpoint_probability::{AdaptiveEndpointProbabilities, EndpointProbabilities};
use crate::util::oarfish_types::InMemoryAlignmentStore;
use serde::Serialize;

const PROBABILITY_FLOOR: f64 = 1e-300;

#[derive(Debug, Serialize)]
pub(crate) struct AdaptiveCoverageDiagnostics {
    pub training_reads: usize,
    pub folds: usize,
    pub selected_prior_mass: f64,
    pub mean_reliability: f64,
    pub capped_reads: usize,
    pub scored_reads: usize,
    pub mean_quality_gate: f64,
    pub mean_uncertainty_gate: f64,
    pub ambiguous_training_reads: usize,
}

fn combine_read_into(
    logistic: &[f64],
    endpoint: &[f64],
    support_gate: f64,
    logistic_weight: f64,
    endpoint_weight: f64,
    output: &mut Vec<f64>,
) {
    if logistic.is_empty() {
        output.clear();
        return;
    }
    output.clear();
    output.extend(logistic.iter().zip(endpoint).map(|(&lp, &ep)| {
        logistic_weight * lp.max(PROBABILITY_FLOOR).ln()
            + endpoint_weight * support_gate.clamp(0.0, 1.0) * ep.max(PROBABILITY_FLOOR).ln()
    }));
    let maximum = output.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    output
        .iter_mut()
        .for_each(|value| *value = (*value - maximum).exp());
    let sum: f64 = output.iter().sum();
    if sum.is_finite() && sum > 0.0 {
        output.iter_mut().for_each(|value| *value /= sum);
    } else {
        let uniform = 1.0 / output.len() as f64;
        output.fill(uniform);
    }
}

#[cfg(test)]
fn combine_read(
    logistic: &[f64],
    endpoint: &[f64],
    support_gate: f64,
    logistic_weight: f64,
    endpoint_weight: f64,
) -> Vec<f64> {
    let mut output = Vec::with_capacity(logistic.len());
    combine_read_into(
        logistic,
        endpoint,
        support_gate,
        logistic_weight,
        endpoint_weight,
        &mut output,
    );
    output
}

pub(crate) fn apply_hybrid_probabilities(
    store: &mut InMemoryAlignmentStore,
    endpoint: EndpointProbabilities,
    logistic_weight: f64,
    endpoint_weight: f64,
) {
    assert_eq!(endpoint.probabilities.len(), store.total_len());
    assert_eq!(endpoint.support_gates.len(), store.len());
    let mut combined = Vec::new();
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        let gate = endpoint.support_gates[read_index] as f64;
        combine_read_into(
            &store.coverage_probabilities[start..end],
            &endpoint.probabilities[start..end],
            gate,
            logistic_weight,
            endpoint_weight,
            &mut combined,
        );
        store.coverage_probabilities[start..end].copy_from_slice(&combined);
    }
}

fn jensen_shannon(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(&l, &r)| {
            let middle = 0.5 * (l + r);
            let contribution = |value: f64| {
                if value > 0.0 {
                    value * (value / middle).ln()
                } else {
                    0.0
                }
            };
            0.5 * (contribution(l) + contribution(r))
        })
        .sum()
}

fn cap_bayes_factor(probabilities: &mut [f64], maximum_bayes_factor: f64) -> bool {
    if probabilities.len() < 2 {
        return false;
    }
    let maximum = probabilities.iter().copied().fold(0.0, f64::max);
    let minimum = probabilities.iter().copied().fold(f64::INFINITY, f64::min);
    if minimum > 0.0 && maximum / minimum <= maximum_bayes_factor {
        return false;
    }
    let floor = maximum / maximum_bayes_factor.max(1.0);
    probabilities
        .iter_mut()
        .for_each(|value| *value = value.max(floor));
    let sum: f64 = probabilities.iter().sum();
    probabilities.iter_mut().for_each(|value| *value /= sum);
    true
}

fn alignment_quality_gate(score_probabilities: &[f32]) -> f64 {
    if score_probabilities.len() < 2 {
        return 1.0;
    }
    let score_sum: f64 = score_probabilities.iter().map(|p| *p as f64).sum();
    if score_sum <= 0.0 {
        return 1.0;
    }
    let best = score_probabilities.iter().copied().fold(0.0_f32, f32::max) as f64;
    let chance = 1.0 / score_probabilities.len() as f64;
    ((best / score_sum - chance) / (1.0 - chance)).clamp(0.0, 1.0)
}

pub(crate) fn apply_adaptive_probabilities(
    store: &mut InMemoryAlignmentStore,
    endpoint: AdaptiveEndpointProbabilities,
    logistic_weight: f64,
    endpoint_weight: f64,
    maximum_bayes_factor: f64,
    ablation: crate::prog_opts::CoverageAblation,
) -> AdaptiveCoverageDiagnostics {
    assert_eq!(endpoint.probabilities.len(), store.total_len());
    assert_eq!(endpoint.support_gates.len(), store.len());
    let mut reliability_sum = 0.0;
    let mut capped_reads = 0usize;
    let mut scored_reads = 0usize;
    let mut quality_sum = 0.0;
    let mut uncertainty_sum = 0.0;
    let mut local = Vec::new();
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        let logistic_read = &store.coverage_probabilities[start..end];
        let endpoint_read = &endpoint.probabilities[start..end];
        let score_probs = &store.as_probabilities[start..end];
        let mut support = endpoint.support_gates[read_index] as f64;
        let disagreement = jensen_shannon(logistic_read, endpoint_read);
        let mut agreement = (1.0 - disagreement / std::f64::consts::LN_2).clamp(0.0, 1.0);
        if ablation == crate::prog_opts::CoverageAblation::NoSupportGate {
            support = 1.0;
        }
        if ablation == crate::prog_opts::CoverageAblation::NoAgreementGate {
            agreement = 1.0;
        }
        let quality = alignment_quality_gate(score_probs);
        let uncertainty = endpoint
            .uncertainty_gates
            .get(read_index)
            .copied()
            .unwrap_or(1.0) as f64;
        let use_quality = matches!(
            ablation,
            crate::prog_opts::CoverageAblation::QualityGate
                | crate::prog_opts::CoverageAblation::AllCandidates
        );
        let use_uncertainty = matches!(
            ablation,
            crate::prog_opts::CoverageAblation::UncertaintyGate
                | crate::prog_opts::CoverageAblation::AllCandidates
        );
        let reliability = support
            * agreement
            * if use_quality { quality } else { 1.0 }
            * if use_uncertainty { uncertainty } else { 1.0 };
        let lw = if ablation == crate::prog_opts::CoverageAblation::NoLogistic {
            0.0
        } else {
            logistic_weight
        };
        let ew = if ablation == crate::prog_opts::CoverageAblation::NoEndpoint {
            0.0
        } else {
            endpoint_weight
        };
        combine_read_into(logistic_read, endpoint_read, 1.0, lw, ew, &mut local);
        let uniform = 1.0 / local.len().max(1) as f64;
        local
            .iter_mut()
            .for_each(|value| *value = (1.0 - reliability) * uniform + reliability * *value);
        if ablation != crate::prog_opts::CoverageAblation::NoBayesCap {
            capped_reads += usize::from(cap_bayes_factor(&mut local, maximum_bayes_factor));
        }
        reliability_sum += reliability;
        scored_reads += 1;
        quality_sum += quality;
        uncertainty_sum += uncertainty;
        store.coverage_probabilities[start..end].copy_from_slice(&local);
    }
    AdaptiveCoverageDiagnostics {
        training_reads: endpoint.training_reads,
        folds: endpoint.folds,
        selected_prior_mass: endpoint.selected_prior_mass,
        mean_reliability: reliability_sum / scored_reads.max(1) as f64,
        capped_reads,
        scored_reads,
        mean_quality_gate: quality_sum / scored_reads.max(1) as f64,
        mean_uncertainty_gate: uncertainty_sum / scored_reads.max(1) as f64,
        ambiguous_training_reads: endpoint.ambiguous_training_reads,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn support_gate_can_make_endpoint_neutral() {
        let result = combine_read(&[0.8, 0.2], &[0.1, 0.9], 0.0, 1.0, 1.0);
        assert!((result[0] - 0.8).abs() < 1e-12);
        assert!((result[1] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn log_linear_combination_is_normalized() {
        let result = combine_read(&[0.8, 0.2], &[0.1, 0.9], 1.0, 0.5, 0.5);
        assert!((result.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(result.iter().all(|value| value.is_finite() && *value > 0.0));
    }

    #[test]
    fn bayes_factor_cap_limits_extreme_evidence() {
        let mut values = vec![0.999, 0.001];
        assert!(cap_bayes_factor(&mut values, 4.0));
        assert!(values[0] / values[1] <= 4.0 + 1e-12);
        assert!((values.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn disagreement_is_zero_only_for_matching_distributions() {
        assert!(jensen_shannon(&[0.8, 0.2], &[0.8, 0.2]).abs() < 1e-12);
        assert!(jensen_shannon(&[0.99, 0.01], &[0.01, 0.99]) > 0.5);
    }

    #[test]
    fn alignment_quality_gate_distinguishes_ties_from_separated_hits() {
        assert!(alignment_quality_gate(&[1.0, 1.0]).abs() < f64::EPSILON);
        assert!(alignment_quality_gate(&[1.0, 0.01]) > 0.98);
        assert_eq!(alignment_quality_gate(&[1.0]), 1.0);
    }
}
