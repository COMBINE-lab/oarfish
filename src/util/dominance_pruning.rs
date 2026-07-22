//! Conservative Pareto-dominance pruning of read-to-transcript candidates.

use crate::util::oarfish_types::InMemoryAlignmentStore;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub(crate) struct DominancePruningDiagnostics {
    pub ambiguous_reads: usize,
    pub reads_pruned: usize,
    pub candidates_pruned: usize,
    pub candidates_before: usize,
    pub minimum_bayes_factor: f64,
}

#[inline]
fn dominates(winner: (f64, f64, u32), loser: (f64, f64, u32), minimum_bayes_factor: f64) -> bool {
    winner.0 >= loser.0
        && winner.1 >= loser.1
        && winner.2 >= loser.2
        && winner.0 * winner.1 >= minimum_bayes_factor * (loser.0 * loser.1).max(f64::MIN_POSITIVE)
}

pub(crate) fn apply_dominance_pruning(
    store: &mut InMemoryAlignmentStore,
    minimum_bayes_factor: f64,
) -> DominancePruningDiagnostics {
    let mut ambiguous_reads = 0;
    let mut reads_pruned = 0;
    let mut candidates_pruned = 0;
    let mut losers = Vec::new();
    let model_coverage = store.filter_opts.model_coverage;
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start < 2 {
            continue;
        }
        ambiguous_reads += 1;
        losers.clear();
        losers.resize(end - start, false);
        for loser in 0..end - start {
            let loser_score = store.as_probabilities[start + loser] as f64;
            let loser_coverage = if model_coverage {
                store.coverage_probabilities[start + loser]
            } else {
                1.0
            };
            let loser_span = store.alignments[start + loser].alignment_span();
            for winner in 0..end - start {
                if winner == loser {
                    continue;
                }
                let winner_score = store.as_probabilities[start + winner] as f64;
                let winner_coverage = if model_coverage {
                    store.coverage_probabilities[start + winner]
                } else {
                    1.0
                };
                let winner_span = store.alignments[start + winner].alignment_span();
                if dominates(
                    (winner_score, winner_coverage, winner_span),
                    (loser_score, loser_coverage, loser_span),
                    minimum_bayes_factor,
                ) {
                    losers[loser] = true;
                    break;
                }
            }
        }
        // The strict Bayes-factor condition should guarantee this, but retain
        // the best joint candidate defensively if numerical underflow marks all.
        if losers.iter().all(|value| *value) {
            let best = (0..end - start)
                .max_by(|&a, &b| {
                    let coverage = |offset| {
                        if model_coverage {
                            store.coverage_probabilities[start + offset]
                        } else {
                            1.0
                        }
                    };
                    let pa = store.as_probabilities[start + a] as f64 * coverage(a);
                    let pb = store.as_probabilities[start + b] as f64 * coverage(b);
                    pa.total_cmp(&pb)
                })
                .unwrap();
            losers[best] = false;
        }
        let local_pruned = losers.iter().filter(|value| **value).count();
        if local_pruned > 0 {
            reads_pruned += 1;
            candidates_pruned += local_pruned;
            for (offset, &loser) in losers.iter().enumerate() {
                if loser {
                    store.as_probabilities[start + offset] = 0.0;
                }
            }
        }
    }
    DominancePruningDiagnostics {
        ambiguous_reads,
        reads_pruned,
        candidates_pruned,
        candidates_before: store.total_len(),
        minimum_bayes_factor,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dominance_requires_every_signal_and_joint_margin() {
        assert!(dominates((1.0, 0.8, 1000), (0.4, 0.5, 900), 2.0));
        assert!(!dominates((1.0, 0.4, 1000), (0.4, 0.5, 900), 2.0));
        assert!(!dominates((1.0, 0.8, 800), (0.4, 0.5, 900), 2.0));
        assert!(!dominates((0.7, 0.8, 1000), (0.4, 0.5, 900), 4.0));
    }
}
