//! Cross-fitted global degradation nuisance model for ONT direct RNA.

use crate::prog_opts::DegradationKernel;
use crate::util::endpoint_probability::{
    adaptive_endpoint_probabilities, endpoint_gaps, AdaptiveEndpointProbabilities,
};
use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use serde::Serialize;

const BINS: usize = 20;
const HAZARD_GRID: [f64; 8] = [0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 2.0, 4.0];

#[derive(Clone, Copy, Debug)]
struct DegradationParameters {
    intact_fraction: f64,
    technical_fraction: f64,
    hazard_per_kb: f64,
    shape: DegradationShape,
    correction_weight: f64,
}

#[derive(Clone, Copy, Debug)]
enum DegradationShape {
    Constant,
    Piecewise2 { rates: [f64; 2] },
}

#[derive(Debug, Serialize)]
pub(crate) struct DegradationDiagnostics {
    pub kernel: &'static str,
    pub training_reads: usize,
    pub folds: usize,
    pub selected_prior_mass: f64,
    pub intact_fraction: f64,
    pub technical_truncation_fraction: f64,
    pub degraded_fraction: f64,
    pub degradation_hazard_per_kb: f64,
    pub degradation_shape: String,
    pub held_out_log_predictive_density: f64,
    pub hazard_fold_sd: f64,
    pub mean_degradation_posterior: f64,
    pub degradation_correction_weight: f64,
    pub mean_degradation_strength: f64,
}

fn constant_hazard_bin_probability(bin: usize, transcript_kb: f64, hazard_per_kb: f64) -> f64 {
    let rate = (transcript_kb * hazard_per_kb).max(1e-9);
    // Sequenced fraction u = 1 - five_prime_gap. Conditional on a break,
    // fragment length follows a truncated exponential from the 3' anchor.
    let x0 = bin as f64 / BINS as f64;
    let x1 = (bin + 1) as f64 / BINS as f64;
    let u0 = 1.0 - x1;
    let u1 = 1.0 - x0;
    ((-rate * u0).exp() - (-rate * u1).exp()) / (1.0 - (-rate).exp())
}

fn cumulative_piecewise_hazard(u: f64, rates: [f64; 2]) -> f64 {
    let width = 0.5;
    let mut remaining = u;
    let mut total = 0.0;
    for rate in rates {
        let span = remaining.min(width).max(0.0);
        total += rate * span;
        remaining -= span;
    }
    total
}

fn degraded_bin_probability(bin: usize, transcript_kb: f64, params: DegradationParameters) -> f64 {
    match params.shape {
        DegradationShape::Constant => {
            constant_hazard_bin_probability(bin, transcript_kb, params.hazard_per_kb)
        }
        DegradationShape::Piecewise2 { rates } => {
            let x0 = bin as f64 / BINS as f64;
            let x1 = (bin + 1) as f64 / BINS as f64;
            let u0 = 1.0 - x1;
            let u1 = 1.0 - x0;
            let scale = transcript_kb * params.hazard_per_kb;
            let survival = |u: f64| (-scale * cumulative_piecewise_hazard(u, rates)).exp();
            (survival(u0) - survival(u1)) / (1.0 - survival(1.0)).max(f64::MIN_POSITIVE)
        }
    }
}

fn technical_bin_probability(bin: usize) -> f64 {
    // Length-invariant Beta(2, 1) retained fraction: ordinary sequencing
    // termination favors longer reads but does not intensify with transcript
    // length in the way a physical per-base degradation hazard does.
    let x0 = bin as f64 / BINS as f64;
    let x1 = (bin + 1) as f64 / BINS as f64;
    let u0 = 1.0 - x1;
    let u1 = 1.0 - x0;
    u1 * u1 - u0 * u0
}

fn degraded_fraction(params: DegradationParameters) -> f64 {
    (1.0 - params.intact_fraction - params.technical_fraction).max(0.0)
}

fn endpoint_probability(bin: usize, transcript_len: usize, params: DegradationParameters) -> f64 {
    let degraded = degraded_bin_probability(bin, transcript_len.max(1) as f64 / 1_000.0, params);
    (degraded_fraction(params) * degraded
        + params.technical_fraction * technical_bin_probability(bin)
        + if bin == 0 {
            params.intact_fraction
        } else {
            0.0
        })
    .max(1e-12)
}

fn degradation_posterior(bin: usize, transcript_len: usize, params: DegradationParameters) -> f64 {
    let degraded = degraded_fraction(params)
        * degraded_bin_probability(bin, transcript_len.max(1) as f64 / 1_000.0, params);
    degraded / endpoint_probability(bin, transcript_len, params)
}

fn degradation_strength(params: DegradationParameters) -> f64 {
    // Normalized Bernoulli variance: near zero when the library is almost
    // entirely intact, and maximal when both intact and degraded molecules are
    // sufficiently represented for degradation correction to be identifiable.
    let degraded = degraded_fraction(params);
    4.0 * degraded * (1.0 - degraded)
}

fn nondegraded_probability(bin: usize, params: DegradationParameters) -> f64 {
    let total = (params.intact_fraction + params.technical_fraction).max(1e-12);
    ((params.technical_fraction * technical_bin_probability(bin)
        + if bin == 0 {
            params.intact_fraction
        } else {
            0.0
        })
        / total)
        .max(1e-12)
}

fn degradation_likelihood_ratio(
    bin: usize,
    transcript_len: usize,
    params: DegradationParameters,
) -> f64 {
    endpoint_probability(bin, transcript_len, params) / nondegraded_probability(bin, params)
}

fn shape_candidates(kernel: DegradationKernel) -> Vec<DegradationShape> {
    match kernel {
        DegradationKernel::Constant => vec![DegradationShape::Constant],
        DegradationKernel::Piecewise2 => vec![
            DegradationShape::Constant,
            DegradationShape::Piecewise2 {
                rates: [0.67, 1.33],
            },
            DegradationShape::Piecewise2 {
                rates: [1.33, 0.67],
            },
        ],
    }
}

fn shape_penalty(shape: DegradationShape) -> f64 {
    match shape {
        DegradationShape::Constant => 0.0,
        DegradationShape::Piecewise2 { rates } => {
            5.0 * rates.iter().map(|rate| rate.ln().powi(2)).sum::<f64>()
        }
    }
}

fn shape_name(shape: DegradationShape) -> String {
    match shape {
        DegradationShape::Constant => "constant".to_string(),
        DegradationShape::Piecewise2 { rates } => format!(
            "piecewise2:{}",
            rates
                .iter()
                .map(|v| format!("{v:.2}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    }
}

fn fit_parameters(
    observations: &[(usize, usize)],
    kernel: DegradationKernel,
) -> DegradationParameters {
    // Weak technology-level priors prevent boundary estimates on small samples.
    // The Dirichlet prior favors intact molecules while retaining technical and
    // degraded components; the log-hazard prior centers at 0.4/kb.
    let score = |params: DegradationParameters| {
        let likelihood = observations
            .iter()
            .map(|&(bin, len)| endpoint_probability(bin, len, params).ln())
            .sum::<f64>();
        let pi_prior = 4.0 * params.intact_fraction.max(1e-12).ln()
            + params.technical_fraction.max(1e-12).ln()
            + degraded_fraction(params).max(1e-12).ln();
        let log_ratio = (params.hazard_per_kb / 0.4).ln();
        likelihood + pi_prior - 0.5 * log_ratio * log_ratio - shape_penalty(params.shape)
    };
    let candidates: Vec<_> = HAZARD_GRID
        .iter()
        .flat_map(|&hazard_per_kb| {
            shape_candidates(kernel).into_iter().map(move |shape| {
                // MAP-EM update with a weak Dirichlet(5,2,2) technology prior.
                // Component probabilities are invariant across EM iterations;
                // cache them so fitting does not repeatedly evaluate exponentials.
                let components: Vec<_> = observations
                    .iter()
                    .map(|&(bin, len)| {
                        (
                            bin == 0,
                            technical_bin_probability(bin),
                            degraded_bin_probability(
                                bin,
                                len.max(1) as f64 / 1_000.0,
                                DegradationParameters {
                                    intact_fraction: 0.7,
                                    technical_fraction: 0.15,
                                    hazard_per_kb,
                                    shape,
                                    correction_weight: 0.0,
                                },
                            ),
                        )
                    })
                    .collect();
                let mut intact_fraction: f64 = 0.7;
                let mut technical_fraction: f64 = 0.15;
                for _ in 0..50 {
                    let degraded = (1.0 - intact_fraction - technical_fraction).max(0.0);
                    let mut intact_responsibility = 0.0;
                    let mut technical_responsibility = 0.0;
                    for &(is_intact_bin, technical_probability, degraded_probability) in &components
                    {
                        let total = (if is_intact_bin { intact_fraction } else { 0.0 })
                            + technical_fraction * technical_probability
                            + degraded * degraded_probability;
                        if is_intact_bin {
                            intact_responsibility += intact_fraction / total;
                        }
                        technical_responsibility +=
                            technical_fraction * technical_probability / total;
                    }
                    let denominator = observations.len() as f64 + 6.0;
                    let updated_intact =
                        ((intact_responsibility + 4.0) / denominator).clamp(0.001, 0.998);
                    let updated_technical = ((technical_responsibility + 1.0) / denominator)
                        .clamp(0.001, 0.998 - updated_intact);
                    if (updated_intact - intact_fraction).abs() < 1e-8
                        && (updated_technical - technical_fraction).abs() < 1e-8
                    {
                        break;
                    }
                    intact_fraction = updated_intact;
                    technical_fraction = updated_technical;
                }
                DegradationParameters {
                    intact_fraction,
                    technical_fraction,
                    hazard_per_kb,
                    shape,
                    correction_weight: 0.0,
                }
            })
        })
        .collect();
    let mut best = candidates
        .iter()
        .copied()
        .max_by(|left, right| score(*left).total_cmp(&score(*right)))
        .unwrap_or(DegradationParameters {
            intact_fraction: 0.8,
            technical_fraction: 0.1,
            hazard_per_kb: 0.4,
            shape: DegradationShape::Constant,
            correction_weight: 0.0,
        });
    let null = candidates[0];
    let log_evidence_gain = (score(best) - score(null)).max(0.0);
    best.correction_weight = if best.hazard_per_kb == HAZARD_GRID[0] {
        0.0
    } else {
        1.0 - (-2.0 * log_evidence_gain).exp()
    };
    best
}

/// Remove sample-wide 3'-anchored degradation from endpoint evidence before
/// candidate comparison. Each read is corrected with parameters learned
/// without its fold. Correction and residual endpoint reliability are weighted
/// by sample-level model evidence, mixture identifiability, and the read-level
/// degradation posterior.
pub(crate) fn degradation_adjusted_endpoint_probabilities(
    store: &InMemoryAlignmentStore,
    txps: &[TranscriptInfo],
    folds: usize,
    support_scale: f64,
    kernel: DegradationKernel,
) -> (AdaptiveEndpointProbabilities, DegradationDiagnostics) {
    let folds = folds.max(2);
    let mut endpoint = adaptive_endpoint_probabilities(store, txps, folds, support_scale);
    let mut by_fold = vec![Vec::new(); folds];
    for (read_index, (alignments, _, _)) in store.iter().enumerate() {
        if let [alignment] = alignments {
            let len = txps[alignment.ref_id as usize].len.get();
            let (five_gap, three_gap) = endpoint_gaps(alignment, len);
            // Direct-RNA degradation is identifiable only for 3'-anchored reads.
            if three_gap <= 0.1 {
                let bin = (five_gap * BINS as f64) as usize;
                by_fold[read_index % folds].push((bin, len));
            }
        }
    }
    let parameters: Vec<_> = (0..folds)
        .map(|held_out| {
            let training: Vec<_> = by_fold
                .iter()
                .enumerate()
                .filter(|(fold, _)| *fold != held_out)
                .flat_map(|(_, observations)| observations.iter().copied())
                .collect();
            fit_parameters(&training, kernel)
        })
        .collect();
    let mut offset = 0usize;
    let mut posterior_sum = 0.0;
    let mut posterior_count = 0usize;
    for (read_index, (alignments, _, _)) in store.iter().enumerate() {
        let params = parameters[read_index % folds];
        let end = offset + alignments.len();
        let local = &mut endpoint.probabilities[offset..end];
        let mut read_posterior_sum = 0.0;
        let mut read_posterior_count = 0usize;
        for (probability, alignment) in local.iter_mut().zip(alignments) {
            let len = txps[alignment.ref_id as usize].len.get();
            let (five_gap, three_gap) = endpoint_gaps(alignment, len);
            if three_gap <= 0.1 {
                let bin = (five_gap * BINS as f64) as usize;
                let posterior = degradation_posterior(bin, len, params);
                *probability /= degradation_likelihood_ratio(bin, len, params).powf(
                    posterior * posterior * degradation_strength(params) * params.correction_weight,
                );
                posterior_sum += posterior;
                posterior_count += 1;
                read_posterior_sum += posterior;
                read_posterior_count += 1;
            }
        }
        let sum: f64 = local.iter().sum();
        if sum.is_finite() && sum > 0.0 {
            local.iter_mut().for_each(|value| *value /= sum);
        } else if !local.is_empty() {
            local.fill(1.0 / local.len() as f64);
        }
        let read_degradation = read_posterior_sum / read_posterior_count.max(1) as f64;
        endpoint.support_gates[offset..end]
            .iter_mut()
            .for_each(|gate| {
                *gate *=
                    1.0 - params.correction_weight * read_degradation * degradation_strength(params)
            });
        offset = end;
    }

    let mean_intact =
        parameters.iter().map(|p| p.intact_fraction).sum::<f64>() / parameters.len().max(1) as f64;
    let mean_hazard =
        parameters.iter().map(|p| p.hazard_per_kb).sum::<f64>() / parameters.len().max(1) as f64;
    let mean_technical = parameters.iter().map(|p| p.technical_fraction).sum::<f64>()
        / parameters.len().max(1) as f64;
    let mean_degraded = parameters
        .iter()
        .map(|&p| degraded_fraction(p))
        .sum::<f64>()
        / parameters.len().max(1) as f64;
    let mean_correction_weight = parameters.iter().map(|p| p.correction_weight).sum::<f64>()
        / parameters.len().max(1) as f64;
    let mean_degradation_strength = parameters
        .iter()
        .map(|&params| degradation_strength(params))
        .sum::<f64>()
        / parameters.len().max(1) as f64;
    let mut predictive_sum = 0.0;
    for (fold, observations) in by_fold.iter().enumerate() {
        for &(bin, len) in observations {
            predictive_sum += endpoint_probability(bin, len, parameters[fold]).ln();
        }
    }
    let mean_log_predictive =
        predictive_sum / by_fold.iter().map(Vec::len).sum::<usize>().max(1) as f64;
    let mean_hazard_square = parameters
        .iter()
        .map(|p| p.hazard_per_kb.powi(2))
        .sum::<f64>()
        / parameters.len().max(1) as f64;
    let diagnostics = DegradationDiagnostics {
        kernel: "ont_drna_3prime_competing_survival_mixture",
        training_reads: by_fold.iter().map(Vec::len).sum(),
        folds,
        selected_prior_mass: endpoint.selected_prior_mass,
        intact_fraction: mean_intact,
        technical_truncation_fraction: mean_technical,
        degraded_fraction: mean_degraded,
        degradation_hazard_per_kb: mean_hazard,
        degradation_shape: shape_name(
            parameters
                .first()
                .map(|p| p.shape)
                .unwrap_or(DegradationShape::Constant),
        ),
        held_out_log_predictive_density: mean_log_predictive,
        hazard_fold_sd: (mean_hazard_square - mean_hazard.powi(2)).max(0.0).sqrt(),
        mean_degradation_posterior: posterior_sum / posterior_count.max(1) as f64,
        degradation_correction_weight: mean_correction_weight,
        mean_degradation_strength,
    };
    (endpoint, diagnostics)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degraded_fragments_favor_short_sequenced_spans() {
        let params = DegradationParameters {
            intact_fraction: 0.6,
            technical_fraction: 0.2,
            hazard_per_kb: 1.0,
            shape: DegradationShape::Constant,
            correction_weight: 1.0,
        };
        let short_span = degraded_bin_probability(18, 2.0, params);
        let long_span = degraded_bin_probability(2, 2.0, params);
        assert!(short_span > long_span);
    }

    #[test]
    fn intact_component_concentrates_in_first_gap_bin() {
        let params = DegradationParameters {
            intact_fraction: 0.9,
            technical_fraction: 0.05,
            hazard_per_kb: 0.4,
            shape: DegradationShape::Constant,
            correction_weight: 1.0,
        };
        assert!(endpoint_probability(0, 2_000, params) > endpoint_probability(1, 2_000, params));
    }

    #[test]
    fn technical_component_is_normalized_and_length_invariant() {
        let total: f64 = (0..BINS).map(technical_bin_probability).sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn challenger_distributions_are_normalized() {
        for kernel in [DegradationKernel::Constant, DegradationKernel::Piecewise2] {
            for shape in shape_candidates(kernel) {
                let params = DegradationParameters {
                    intact_fraction: 0.6,
                    technical_fraction: 0.2,
                    hazard_per_kb: 0.7,
                    shape,
                    correction_weight: 1.0,
                };
                let total: f64 = (0..BINS)
                    .map(|bin| degraded_bin_probability(bin, 2.0, params))
                    .sum();
                assert!((total - 1.0).abs() < 1e-10, "{kernel:?} {shape:?}");
            }
        }
    }

    #[test]
    fn piecewise_challenger_embeds_constant_hazard() {
        assert!(shape_candidates(DegradationKernel::Piecewise2)
            .iter()
            .any(|shape| matches!(shape, DegradationShape::Constant)));
    }

    fn simulated_observations(shape: DegradationShape) -> Vec<(usize, usize)> {
        let mut observations = Vec::new();
        for len in [1_000usize, 5_000] {
            let params = DegradationParameters {
                intact_fraction: 0.1,
                technical_fraction: 0.1,
                hazard_per_kb: 0.7,
                shape,
                correction_weight: 1.0,
            };
            for bin in 0..BINS {
                let count = (2_000.0 * endpoint_probability(bin, len, params)).round() as usize;
                observations.extend(std::iter::repeat_n((bin, len), count));
            }
        }
        observations
    }

    #[test]
    fn regularized_piecewise_model_collapses_on_constant_data() {
        let fitted = fit_parameters(
            &simulated_observations(DegradationShape::Constant),
            DegradationKernel::Piecewise2,
        );
        assert!(matches!(fitted.shape, DegradationShape::Constant));
    }

    #[test]
    fn piecewise_model_recovers_directional_hazard_shape() {
        let shape = DegradationShape::Piecewise2 {
            rates: [0.67, 1.33],
        };
        let fitted = fit_parameters(
            &simulated_observations(shape),
            DegradationKernel::Piecewise2,
        );
        assert!(matches!(fitted.shape, DegradationShape::Piecewise2 { .. }));
    }

    #[test]
    fn degradation_correction_depends_on_physical_transcript_length() {
        let params = DegradationParameters {
            intact_fraction: 0.6,
            technical_fraction: 0.2,
            hazard_per_kb: 0.7,
            shape: DegradationShape::Constant,
            correction_weight: 1.0,
        };
        assert_ne!(
            degradation_likelihood_ratio(10, 1_000, params),
            degradation_likelihood_ratio(10, 5_000, params)
        );
    }

    #[test]
    fn fitted_intact_fraction_falls_for_truncated_observations() {
        let intact = fit_parameters(&vec![(0, 2_000); 100], DegradationKernel::Constant);
        let degraded = fit_parameters(&vec![(18, 2_000); 100], DegradationKernel::Constant);
        assert!(intact.intact_fraction > degraded.intact_fraction);
    }
}
