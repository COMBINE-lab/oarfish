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

/// Number of grid cells that can ever hold an observation.
///
/// `five_gap + three_gap = 1 - aligned_fraction <= 1`, and
/// `floor(a) + floor(b) <= floor(a + b)`, so a valid alignment always lands in a
/// cell with `x + y <= GRID - 1`. That is the triangle `GRID * (GRID + 1) / 2`
/// = 210 of 400 cells; the other 190 are structurally unreachable.
const FEASIBLE_CELLS: usize = GRID * (GRID + 1) / 2;

#[inline]
const fn is_feasible_cell(cell: usize) -> bool {
    cell / GRID + cell % GRID < GRID
}

/// Which endpoint-grid geometry corrections are active.
///
/// The default is the historical behavior, so every field being `false` must
/// reproduce the legacy model exactly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct GridConfig {
    /// Restrict the 3x3 smoother to reachable neighbours.
    pub feasible_smoothing: bool,
    /// Spread the Dirichlet prior over reachable cells only.
    pub feasible_prior: bool,
    /// Compare candidates as densities over integer gap pairs rather than as
    /// raw cell probabilities.
    pub nt_measure: bool,
}

impl GridConfig {
    pub(crate) fn from_ablation(ablation: crate::prog_opts::CoverageAblation) -> Self {
        use crate::prog_opts::CoverageAblation as A;
        match ablation {
            A::EndpointFeasibleSmoothing => Self {
                feasible_smoothing: true,
                ..Self::default()
            },
            A::EndpointFeasiblePrior => Self {
                feasible_prior: true,
                ..Self::default()
            },
            A::EndpointFeasibleGeometry => Self {
                feasible_smoothing: true,
                feasible_prior: true,
                ..Self::default()
            },
            A::EndpointNtMeasure => Self {
                nt_measure: true,
                ..Self::default()
            },
            A::EndpointNtMeasureGeometry => Self {
                feasible_smoothing: true,
                feasible_prior: true,
                nt_measure: true,
            },
            _ => Self::default(),
        }
    }

    /// Fraction of the total prior mass allocated to `cell`.
    #[inline]
    fn prior_share(&self, cell: usize) -> f64 {
        if self.feasible_prior {
            if is_feasible_cell(cell) {
                1.0 / FEASIBLE_CELLS as f64
            } else {
                0.0
            }
        } else {
            1.0 / CELLS as f64
        }
    }
}

/// Number of integer gap values that fall in `bin` on a transcript of `len`.
///
/// `EndpointModel::cell` assigns gap `g` to `floor(g / len * GRID)`, so bin `b`
/// holds the integers in `[ceil(b*len/GRID), ceil((b+1)*len/GRID))`.
#[inline]
fn axis_lattice_width(bin: usize, len: usize) -> f64 {
    let lo = (bin * len).div_ceil(GRID);
    let hi = if bin + 1 < GRID {
        ((bin + 1) * len).div_ceil(GRID)
    } else {
        len
    };
    hi.saturating_sub(lo).max(1) as f64
}

/// How many integer `(five_gap, three_gap)` pairs a cell covers.
///
/// A cell is a *probability*, but candidates of different length spread that
/// probability over different numbers of lattice points. Comparing candidates
/// therefore requires dividing by this area to obtain a common measure.
#[inline]
fn cell_lattice_area(cell: usize, len: usize) -> f64 {
    axis_lattice_width(cell / GRID, len) * axis_lattice_width(cell % GRID, len)
}

/// Smoothed-count-plus-prior estimate for a single cell.
///
/// With a default [`GridConfig`] this is exactly the historical expression
/// `(count + prior / CELLS) / (total + prior)`.
#[inline]
fn cell_probability(count: f64, total: f64, cell: usize, prior: f64, config: GridConfig) -> f64 {
    (count + prior * config.prior_share(cell)) / (total + prior)
}

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

fn smooth_counts(
    counts: &[[f64; CELLS]; STRATA],
    config: GridConfig,
) -> ([[f64; CELLS]; STRATA], [f64; STRATA]) {
    let mut smoothed = [[0.0; CELLS]; STRATA];
    let mut totals = [0.0; STRATA];
    for stratum in 0..STRATA {
        for x in 0..GRID {
            for y in 0..GRID {
                let cell = x * GRID + y;
                // An unreachable cell has no mass to redistribute and must not
                // borrow any; leaving it at zero keeps the smoother's support
                // equal to the model's support.
                if config.feasible_smoothing && !is_feasible_cell(cell) {
                    continue;
                }
                let mut weighted = 0.0;
                let mut weight_sum = 0.0;
                for nx in x.saturating_sub(1)..=(x + 1).min(GRID - 1) {
                    for ny in y.saturating_sub(1)..=(y + 1).min(GRID - 1) {
                        // Structurally-empty neighbours contribute a zero count
                        // but, on the legacy path, still inflate `weight_sum`.
                        // That deflates anti-diagonal cells by a cell-specific
                        // factor no scalar prior can undo.
                        if config.feasible_smoothing && !is_feasible_cell(nx * GRID + ny) {
                            continue;
                        }
                        let weight =
                            if nx == x { 2.0 } else { 1.0 } * if ny == y { 2.0 } else { 1.0 };
                        weighted += weight * counts[stratum][nx * GRID + ny];
                        weight_sum += weight;
                    }
                }
                smoothed[stratum][cell] = weighted / weight_sum;
                totals[stratum] += smoothed[stratum][cell];
            }
        }
    }
    if config.feasible_smoothing {
        // Renormalize so the correction changes only the *shape* of the
        // smoothed table. Without this the restricted weight sum also raises
        // total observed mass, which would confound the geometry change with a
        // change in the data-to-prior ratio.
        for stratum in 0..STRATA {
            let observed: f64 = counts[stratum].iter().sum();
            if totals[stratum] > 0.0 && observed > 0.0 {
                let scale = observed / totals[stratum];
                smoothed[stratum].iter_mut().for_each(|value| *value *= scale);
                totals[stratum] = observed;
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
    // The default grid's lower bound is a binding constraint on real ONT data:
    // the held-out selector picks 10.0 in every library measured. This extends
    // the search downward so the optimum can be interior.
    const WIDE_PRIOR_GRID: [f64; 9] = [
        1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1_000.0, 3_000.0, 10_000.0,
    ];
    let prior_grid: &[f64] = if ablation == crate::prog_opts::CoverageAblation::EndpointWidePriorGrid
    {
        &WIDE_PRIOR_GRID
    } else {
        &PRIOR_GRID
    };
    let folds = folds.max(2);
    let grid_config = GridConfig::from_ablation(ablation);
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
        smooth_counts(&counts, grid_config)
    };
    let fitted: Vec<_> = (0..folds).map(training_for_fold).collect();
    let selected_prior_mass = prior_grid
        .iter()
        .copied()
        .max_by(|left, right| {
            let score = |prior: f64| {
                observations
                    .iter()
                    .map(|&(fold, stratum, cell)| {
                        let (counts, totals) = &fitted[fold];
                        cell_probability(
                            counts[stratum][cell],
                            totals[stratum],
                            cell,
                            prior,
                            grid_config,
                        )
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
            let mut value = cell_probability(
                counts[stratum][cell],
                totals[stratum],
                cell,
                selected_prior_mass,
                grid_config,
            );
            if grid_config.nt_measure {
                // Candidates differ in length, so the same cell covers a
                // different number of integer gap pairs for each. Without this
                // the comparison is between probabilities of differently-sized
                // events.
                value /= cell_lattice_area(cell, len);
            }
            local.push(value);
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
                    let value = cell_probability(
                        fold_counts[stratum][cell],
                        fold_totals[stratum],
                        cell,
                        selected_prior_mass,
                        grid_config,
                    );
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
        let (smoothed, _) = smooth_counts(&counts, GridConfig::default());
        assert!(smoothed[0][10 * GRID + 10] > 0.0);
        assert!(smoothed[0][10 * GRID + 9] > 0.0);
        assert_eq!(smoothed[0][0], 0.0);
    }

    /// The reachable region is the triangle `x + y < GRID`, and a valid
    /// alignment can never land outside it.
    #[test]
    fn feasible_region_is_the_expected_triangle() {
        assert_eq!(FEASIBLE_CELLS, 210);
        assert_eq!((0..CELLS).filter(|c| is_feasible_cell(*c)).count(), 210);
        assert!(is_feasible_cell(19)); // x=0, y=19
        assert!(!is_feasible_cell(GRID + 19)); // x=1, y=19
    }

    #[test]
    fn valid_alignments_never_reach_an_infeasible_cell() {
        let len = 1_000usize;
        for start in [1u32, 7, 250, 500, 999] {
            for end in [start + 1, start + 50, 999, 1_000] {
                if end as usize > len {
                    continue;
                }
                let cell = EndpointModel::cell(&aln(start, end), len);
                assert!(is_feasible_cell(cell), "start={start} end={end}");
            }
        }
    }

    /// The smoother must not borrow mass across the feasibility boundary.
    #[test]
    fn feasible_smoothing_ignores_unreachable_neighbors() {
        let config = GridConfig {
            feasible_smoothing: true,
            ..GridConfig::default()
        };
        let mut counts = [[0.0; CELLS]; STRATA];
        // A cell on the anti-diagonal: x + y == GRID - 1.
        let boundary = 10 * GRID + 9;
        counts[0][boundary] = 16.0;
        let (legacy, _) = smooth_counts(&counts, GridConfig::default());
        let (fixed, _) = smooth_counts(&counts, config);
        // Legacy divides by a weight sum inflated by structurally-empty
        // neighbours, so the boundary cell is deflated relative to the fix.
        assert!(fixed[0][boundary] > legacy[0][boundary]);
        // Unreachable cells receive no mass at all under the fix.
        for cell in 0..CELLS {
            if !is_feasible_cell(cell) {
                assert_eq!(fixed[0][cell], 0.0, "cell {cell}");
            }
        }
    }

    /// The correction must change only the shape of the table, not its mass.
    #[test]
    fn feasible_smoothing_preserves_total_mass() {
        let config = GridConfig {
            feasible_smoothing: true,
            ..GridConfig::default()
        };
        let mut counts = [[0.0; CELLS]; STRATA];
        for cell in 0..CELLS {
            if is_feasible_cell(cell) {
                counts[0][cell] = (cell % 7) as f64;
            }
        }
        let observed: f64 = counts[0].iter().sum();
        let (_, totals) = smooth_counts(&counts, config);
        assert!((totals[0] - observed).abs() < 1e-9);
    }

    /// Deep-interior cells must be untouched by the geometry correction.
    #[test]
    fn feasible_smoothing_leaves_interior_cells_alone() {
        let config = GridConfig {
            feasible_smoothing: true,
            ..GridConfig::default()
        };
        let mut counts = [[0.0; CELLS]; STRATA];
        counts[0][2 * GRID + 3] = 12.0;
        let (legacy, _) = smooth_counts(&counts, GridConfig::default());
        let (fixed, _) = smooth_counts(&counts, config);
        // Mass is far from both the grid edge and the anti-diagonal, so the two
        // smoothers agree up to the mass-preserving rescale.
        for (x, y) in [(2usize, 3usize), (2, 4), (3, 3), (1, 2)] {
            let cell = x * GRID + y;
            assert!(
                (legacy[0][cell] - fixed[0][cell]).abs() < 1e-12,
                "cell ({x},{y}): {} vs {}",
                legacy[0][cell],
                fixed[0][cell]
            );
        }
    }

    /// A default config must reproduce the historical estimator exactly.
    #[test]
    fn default_grid_config_matches_legacy_probability() {
        let config = GridConfig::default();
        for cell in [0usize, 5, 199, 399] {
            let legacy = (3.0 + 1_000.0 / CELLS as f64) / (500.0 + 1_000.0);
            assert_eq!(cell_probability(3.0, 500.0, cell, 1_000.0, config), legacy);
        }
    }

    /// Restricting the prior to reachable cells raises the per-cell share by
    /// exactly `CELLS / FEASIBLE_CELLS` and zeroes the unreachable ones.
    #[test]
    fn feasible_prior_concentrates_mass_on_reachable_cells() {
        let config = GridConfig {
            feasible_prior: true,
            ..GridConfig::default()
        };
        let total: f64 = (0..CELLS).map(|c| config.prior_share(c)).sum();
        assert!((total - 1.0).abs() < 1e-12);
        assert_eq!(config.prior_share(GRID + 19), 0.0);
        let ratio = config.prior_share(0) / GridConfig::default().prior_share(0);
        assert!((ratio - CELLS as f64 / FEASIBLE_CELLS as f64).abs() < 1e-12);
    }
}
