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
}

fn combine_read(
    logistic: &[f64],
    endpoint: &[f64],
    support_gate: f64,
    logistic_weight: f64,
    endpoint_weight: f64,
) -> Vec<f64> {
    if logistic.is_empty() {
        return Vec::new();
    }
    let mut logs: Vec<f64> = logistic
        .iter()
        .zip(endpoint)
        .map(|(&lp, &ep)| {
            logistic_weight * lp.max(PROBABILITY_FLOOR).ln()
                + endpoint_weight * support_gate.clamp(0.0, 1.0) * ep.max(PROBABILITY_FLOOR).ln()
        })
        .collect();
    let maximum = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    logs.iter_mut()
        .for_each(|value| *value = (*value - maximum).exp());
    let sum: f64 = logs.iter().sum();
    if sum.is_finite() && sum > 0.0 {
        logs.iter_mut().for_each(|value| *value /= sum);
    } else {
        let uniform = 1.0 / logs.len() as f64;
        logs.fill(uniform);
    }
    logs
}

pub(crate) fn apply_hybrid_probabilities(
    store: &mut InMemoryAlignmentStore,
    logistic: &[f64],
    endpoint: EndpointProbabilities,
    logistic_weight: f64,
    endpoint_weight: f64,
) {
    assert_eq!(logistic.len(), store.total_len());
    assert_eq!(endpoint.probabilities.len(), store.total_len());
    assert_eq!(endpoint.support_gates.len(), store.total_len());
    let mut combined = Vec::with_capacity(store.total_len());
    let mut offset = 0usize;
    for (alignments, _, _) in store.iter() {
        let end = offset + alignments.len();
        let gate = endpoint.support_gates.get(offset).copied().unwrap_or(0.0);
        combined.extend(combine_read(
            &logistic[offset..end],
            &endpoint.probabilities[offset..end],
            gate,
            logistic_weight,
            endpoint_weight,
        ));
        offset = end;
    }
    store.coverage_probabilities = combined;
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

pub(crate) fn apply_adaptive_probabilities(
    store: &mut InMemoryAlignmentStore,
    logistic: &[f64],
    endpoint: AdaptiveEndpointProbabilities,
    logistic_weight: f64,
    endpoint_weight: f64,
    maximum_bayes_factor: f64,
) -> AdaptiveCoverageDiagnostics {
    assert_eq!(logistic.len(), store.total_len());
    assert_eq!(endpoint.probabilities.len(), store.total_len());
    let mut combined = Vec::with_capacity(store.total_len());
    let mut offset = 0usize;
    let mut reliability_sum = 0.0;
    let mut capped_reads = 0usize;
    let mut scored_reads = 0usize;
    for (alignments, _, _) in store.iter() {
        let end = offset + alignments.len();
        let logistic_read = &logistic[offset..end];
        let endpoint_read = &endpoint.probabilities[offset..end];
        let support = endpoint.support_gates.get(offset).copied().unwrap_or(0.0);
        let disagreement = jensen_shannon(logistic_read, endpoint_read);
        let agreement = (1.0 - disagreement / std::f64::consts::LN_2).clamp(0.0, 1.0);
        let reliability = support * agreement;
        let mut local = combine_read(
            logistic_read,
            endpoint_read,
            1.0,
            logistic_weight,
            endpoint_weight,
        );
        let uniform = 1.0 / local.len().max(1) as f64;
        local
            .iter_mut()
            .for_each(|value| *value = (1.0 - reliability) * uniform + reliability * *value);
        capped_reads += usize::from(cap_bayes_factor(&mut local, maximum_bayes_factor));
        reliability_sum += reliability;
        scored_reads += 1;
        combined.extend(local);
        offset = end;
    }
    store.coverage_probabilities = combined;
    AdaptiveCoverageDiagnostics {
        training_reads: endpoint.training_reads,
        folds: endpoint.folds,
        selected_prior_mass: endpoint.selected_prior_mass,
        mean_reliability: reliability_sum / scored_reads.max(1) as f64,
        capped_reads,
        scored_reads,
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
}
