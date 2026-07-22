//! Read-level likelihood for alignment-induced endpoint censoring.

use crate::util::oarfish_types::{AlnInfo, InMemoryAlignmentStore, TranscriptInfo};
use serde::Serialize;

const PRIOR_OBSERVATIONS: f64 = 100.0;
const PRIOR_SCALE_NT: f64 = 50.0;
const MAX_BAYES_FACTOR: f64 = 4.0;

#[derive(Debug, Serialize)]
pub(crate) struct CensoringDiagnostics {
    pub training_reads: usize,
    pub learned_scale_nt: f64,
    pub reliability: f64,
    pub ambiguous_reads: usize,
    pub informative_reads: usize,
    pub mean_unexplained_clip_nt: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio_types::strand::Strand;

    #[test]
    fn only_clipping_with_available_transcript_sequence_is_unexplained() {
        let alignment = AlnInfo {
            ref_id: 0,
            start: 101,
            end: 900,
            prob: 0.0,
            strand: Strand::Forward,
            left_clip: 40,
            right_clip: 150,
        };
        assert_eq!(unexplained_clip(&alignment, 1_000), 140);
    }
}

#[inline]
fn unexplained_clip(alignment: &AlnInfo, transcript_length: usize) -> u32 {
    let left_gap = alignment.start.saturating_sub(1);
    let right_gap = transcript_length.saturating_sub(alignment.end as usize) as u32;
    alignment.left_clip.min(left_gap) + alignment.right_clip.min(right_gap)
}

pub(crate) fn estimate_censor_scale(
    store: &InMemoryAlignmentStore,
    transcripts: &[TranscriptInfo],
) -> f64 {
    let mut training_reads = 0usize;
    let mut training_clip = 0u64;
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start == 1 {
            let alignment = &store.alignments[start];
            training_clip +=
                unexplained_clip(alignment, transcripts[alignment.ref_id as usize].len.get())
                    as u64;
            training_reads += 1;
        }
    }
    ((training_clip as f64 + PRIOR_OBSERVATIONS * PRIOR_SCALE_NT)
        / (training_reads as f64 + PRIOR_OBSERVATIONS))
        .clamp(25.0, 500.0)
}

pub(crate) fn apply_censoring_probabilities(
    store: &mut InMemoryAlignmentStore,
    transcripts: &[TranscriptInfo],
) -> CensoringDiagnostics {
    let mut training_reads = 0usize;
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start == 1 {
            training_reads += 1;
        }
    }
    let scale = estimate_censor_scale(store, transcripts);
    let reliability = training_reads as f64 / (training_reads as f64 + 1_000.0);
    let mut ambiguous_reads = 0;
    let mut informative_reads = 0;
    let mut total_clip = 0u64;
    let mut factors = Vec::new();
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start < 2 {
            continue;
        }
        ambiguous_reads += 1;
        factors.clear();
        factors.extend(store.alignments[start..end].iter().map(|alignment| {
            let clip =
                unexplained_clip(alignment, transcripts[alignment.ref_id as usize].len.get());
            total_clip += clip as u64;
            (-(clip as f64) / scale).exp().powf(reliability)
        }));
        let maximum = factors.iter().copied().fold(0.0, f64::max);
        let floor = maximum / MAX_BAYES_FACTOR;
        let mut informative = false;
        for (probability, factor) in store.as_probabilities[start..end].iter_mut().zip(&factors) {
            let factor = factor.max(floor);
            informative |= factor < maximum * (1.0 - 1e-12);
            *probability *= factor as f32;
        }
        informative_reads += usize::from(informative);
    }
    let candidates = store
        .boundaries
        .windows(2)
        .filter(|window| window[1] - window[0] > 1)
        .map(|window| window[1] - window[0])
        .sum::<usize>();
    CensoringDiagnostics {
        training_reads,
        learned_scale_nt: scale,
        reliability,
        ambiguous_reads,
        informative_reads,
        mean_unexplained_clip_nt: total_clip as f64 / candidates.max(1) as f64,
    }
}
