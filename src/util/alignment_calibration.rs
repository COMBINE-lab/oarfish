//! Agreement-calibrated alignment-score likelihoods.

use crate::util::censoring_probability::unexplained_clip;
use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub(crate) struct AlignmentCalibrationDiagnostics {
    pub ambiguous_reads: usize,
    pub calibrated_reads: usize,
    pub mean_exponent: f64,
    pub mean_candidate_count: f64,
}

#[inline]
fn agreement_exponent(top: f32, second: f32, candidates: usize) -> f64 {
    let margin = 1.0 - second as f64 / top as f64;
    let multiplicity = ((candidates - 1) as f64).sqrt();
    1.0 + 4.0 * (1.0 - margin).powi(2) / multiplicity
}

pub(crate) fn apply_agreement_calibration(
    store: &mut InMemoryAlignmentStore<'_>,
    transcripts: &[TranscriptInfo],
) -> AlignmentCalibrationDiagnostics {
    let mut ambiguous_reads = 0usize;
    let mut calibrated_reads = 0usize;
    let mut exponent_sum = 0.0;
    let mut candidate_sum = 0usize;

    for read in 0..store.len() {
        let start = store.boundaries[read];
        let end = store.boundaries[read + 1];
        let candidates = end - start;
        if candidates < 2 {
            continue;
        }
        ambiguous_reads += 1;
        candidate_sum += candidates;

        let probabilities = &store.as_probabilities[start..end];
        let mut top = 0usize;
        let mut second = 0.0f32;
        for (candidate, &probability) in probabilities.iter().enumerate() {
            if probability > probabilities[top] {
                second = probabilities[top];
                top = candidate;
            } else if candidate != top && probability > second {
                second = probability;
            }
        }
        let top_probability = probabilities[top];
        if top_probability <= 0.0 {
            continue;
        }

        let mut minimum_clip = u32::MAX;
        let mut maximum_span = 0u32;
        for alignment in &store.alignments[start..end] {
            minimum_clip = minimum_clip.min(unexplained_clip(
                alignment,
                transcripts[alignment.ref_id as usize].len.get(),
            ));
            maximum_span = maximum_span.max(alignment.alignment_span());
        }
        let top_alignment = &store.alignments[start + top];
        let top_clip = unexplained_clip(
            top_alignment,
            transcripts[top_alignment.ref_id as usize].len.get(),
        );
        if top_clip != minimum_clip || top_alignment.alignment_span() != maximum_span {
            continue;
        }

        // Sharpen only when the score winner agrees with two independent read
        // features. Near-ties receive more help; already decisive scores are
        // left close to their original temperature. Larger candidate sets are
        // shrunk toward the unmodified likelihood.
        let exponent = agreement_exponent(top_probability, second, candidates);
        if exponent > 1.0 + 1e-12 {
            for probability in &mut store.as_probabilities[start..end] {
                *probability = probability.powf(exponent as f32);
            }
            exponent_sum += exponent;
            calibrated_reads += 1;
        }
    }

    AlignmentCalibrationDiagnostics {
        ambiguous_reads,
        calibrated_reads,
        mean_exponent: exponent_sum / calibrated_reads.max(1) as f64,
        mean_candidate_count: candidate_sum as f64 / ambiguous_reads.max(1) as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::agreement_exponent;

    #[test]
    fn agreement_sharpens_ties_more_than_decisive_scores() {
        let tied = agreement_exponent(1.0, 0.9, 2);
        let decisive = agreement_exponent(1.0, 0.1, 2);
        assert!(tied > decisive);
        assert!(decisive > 1.0);
    }

    #[test]
    fn multiplicity_regularizes_calibration_strength() {
        assert!(agreement_exponent(1.0, 0.9, 2) > agreement_exponent(1.0, 0.9, 10));
    }
}
