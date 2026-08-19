//! Conditional 3'-completeness likelihood for reads that carry a poly(A) tail.
//!
//! A read whose terminal soft clip is a poly(A) tail demonstrably reached its
//! transcript's 3' end: it was not truncated there. The marginal endpoint model
//! cannot express this, because it must accommodate the many reads that *are*
//! truncated, so a large 3' gap carries only a mild penalty. Conditioning on the
//! tail licenses a far sharper penalty for exactly the reads known to be
//! complete.
//!
//! The reference distribution is learned from single-candidate reads that carry
//! a tail, so the tolerance and the rate of genuine disagreement (alternative
//! polyadenylation, mis-annotated 3' ends) are both estimated from the sample
//! rather than assumed.

use crate::util::oarfish_types::{AlnInfo, InMemoryAlignmentStore, TranscriptInfo};
use bio_types::strand::Strand;
use serde::Serialize;

/// Upper edge (exclusive) of each 3'-gap bin, in nucleotides. Log-spaced, so
/// resolution is finest where a complete read's gap actually falls.
const GAP_EDGES: [u32; 13] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, u32::MAX];
const PRIOR_MASS: f64 = 20.0;
/// Cap on how much the term may favour one candidate over another, matching the
/// spirit of `--coverage-max-bayes-factor`: annotation error must not be able to
/// eliminate a true candidate outright.
const MAX_ODDS: f64 = 20.0;

#[derive(Debug, Serialize)]
pub(crate) struct PolyaDiagnostics {
    pub training_reads: usize,
    pub scored_reads: usize,
    pub tail_reads: usize,
    pub flush_fraction: f64,
    pub diffuse_fraction: f64,
    pub mean_log_odds_range: f64,
}

#[inline]
fn three_prime_gap_nt(aln: &AlnInfo, len: usize) -> u32 {
    let left = aln.start.saturating_sub(1);
    let right = len.saturating_sub(aln.end as usize) as u32;
    if aln.strand == Strand::Reverse { left } else { right }
}

#[inline]
fn gap_bin(gap: u32) -> usize {
    GAP_EDGES.iter().position(|edge| gap < *edge).unwrap_or(GAP_EDGES.len() - 1)
}

/// Width of a bin in nucleotides, so bins of unequal width are compared as
/// densities rather than as raw probabilities.
fn bin_width(bin: usize) -> f64 {
    let hi = GAP_EDGES[bin];
    let lo = if bin == 0 { 0 } else { GAP_EDGES[bin - 1] };
    // The final bin is unbounded; give it the span it can realistically cover
    // rather than u32::MAX, which would annihilate its density.
    if hi == u32::MAX {
        4096.0
    } else {
        (hi - lo) as f64
    }
}

/// Multiply the conditional 3'-completeness likelihood into the coverage term
/// for every read that carries a poly(A) tail. Reads without a tail are left
/// untouched: absence of a detected tail is not evidence of truncation.
pub(crate) fn apply_polya_probabilities(
    store: &mut InMemoryAlignmentStore,
    txps: &[TranscriptInfo],
) -> PolyaDiagnostics {
    let bins = GAP_EDGES.len();
    let mut counts = vec![0.0f64; bins];
    let mut training = 0usize;

    // Train only on unambiguous tail-bearing reads: an ambiguous read must not
    // manufacture the distribution used to adjudicate it.
    for (read_index, (alns, _, _)) in store.iter().enumerate() {
        if store.polya_tail[read_index] == 0 {
            continue;
        }
        if let [aln] = alns {
            let len = txps[aln.ref_id as usize].len.get();
            counts[gap_bin(three_prime_gap_nt(aln, len))] += 1.0;
            training += 1;
        }
    }
    let total: f64 = counts.iter().sum();
    let density = |bin: usize| -> f64 {
        ((counts[bin] + PRIOR_MASS / bins as f64) / (total + PRIOR_MASS)) / bin_width(bin)
    };

    let mut scored = 0usize;
    let mut tail_reads = 0usize;
    let mut odds_range_sum = 0.0;
    let mut local: Vec<f64> = Vec::new();
    for read_index in 0..store.len() {
        if store.polya_tail[read_index] == 0 {
            continue;
        }
        tail_reads += 1;
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start < 2 {
            continue;
        }
        local.clear();
        for aln in &store.alignments[start..end] {
            let len = txps[aln.ref_id as usize].len.get();
            local.push(density(gap_bin(three_prime_gap_nt(aln, len))));
        }
        let max = local.iter().copied().fold(0.0, f64::max);
        if !(max > 0.0) {
            continue;
        }
        // Cap the odds before normalizing so no candidate is eliminated outright.
        let floor = max / MAX_ODDS;
        local.iter_mut().for_each(|v| *v = v.max(floor));
        let sum: f64 = local.iter().sum();
        if !(sum > 0.0) {
            continue;
        }
        let min = local.iter().copied().fold(f64::INFINITY, f64::min);
        odds_range_sum += (max.max(floor) / min).ln();
        for (probability, value) in store.coverage_probabilities[start..end]
            .iter_mut()
            .zip(&local)
        {
            *probability *= value / sum;
        }
        // Renormalize this read's coverage term.
        let renorm: f64 = store.coverage_probabilities[start..end].iter().sum();
        if renorm > 0.0 && renorm.is_finite() {
            store.coverage_probabilities[start..end]
                .iter_mut()
                .for_each(|v| *v /= renorm);
        } else {
            let uniform = 1.0 / (end - start) as f64;
            store.coverage_probabilities[start..end].fill(uniform);
        }
        scored += 1;
    }

    let flush = if total > 0.0 {
        counts[..gap_bin(20)].iter().sum::<f64>() / total
    } else {
        0.0
    };
    let diffuse = if total > 0.0 {
        counts[gap_bin(256)..].iter().sum::<f64>() / total
    } else {
        0.0
    };
    PolyaDiagnostics {
        training_reads: training,
        scored_reads: scored,
        tail_reads,
        flush_fraction: flush,
        diffuse_fraction: diffuse,
        mean_log_odds_range: odds_range_sum / scored.max(1) as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(start: u32, end: u32, strand: Strand) -> AlnInfo {
        AlnInfo { ref_id: 0, start, end, strand, left_clip: 0, right_clip: 0 }
    }

    #[test]
    fn three_prime_gap_is_strand_aware() {
        assert_eq!(three_prime_gap_nt(&aln(1, 900, Strand::Forward), 1000), 100);
        assert_eq!(three_prime_gap_nt(&aln(101, 1000, Strand::Reverse), 1000), 100);
        assert_eq!(three_prime_gap_nt(&aln(1, 1000, Strand::Forward), 1000), 0);
    }

    #[test]
    fn gap_bins_are_monotone_and_saturating() {
        assert_eq!(gap_bin(0), 0);
        assert!(gap_bin(0) < gap_bin(5));
        assert!(gap_bin(5) < gap_bin(500));
        assert_eq!(gap_bin(u32::MAX), GAP_EDGES.len() - 1);
    }

    #[test]
    fn bin_widths_are_positive_and_finite() {
        for bin in 0..GAP_EDGES.len() {
            let w = bin_width(bin);
            assert!(w > 0.0 && w.is_finite(), "bin {bin} width {w}");
        }
    }
}

/// Soft-clip bases at the molecule-3' side that a genuine 3'-complete read
/// legitimately carries: partially basecalled poly(A) plus adapter.
const OVERHANG_ALLOWANCE_NT: u32 = 80;

/// Signed, clip-aware 3' mismatch. Two failure directions, both evidence
/// against a candidate:
/// - undershoot: the annotated 3' end lies beyond the read even after
///   crediting the clipped ragged read-start (gap − covered);
/// - overhang: the read extends past the candidate's annotated 3' end (the
///   excess is soft-clipped), beyond the poly(A)/adapter allowance — the
///   signature of a candidate shorter than the molecule (e.g. an alt-3'
///   decoy of the true isoform).
/// A symmetric clip *credit* alone is blind to overhang: alignment gap is 0
/// at a too-short candidate's end and the excess hides in the clip.
#[inline]
fn three_prime_gap_clip_aware_nt(aln: &AlnInfo, len: usize) -> u32 {
    let gap = three_prime_gap_nt(aln, len);
    let clip3 = if aln.strand == Strand::Reverse {
        aln.left_clip
    } else {
        aln.right_clip
    };
    let covered = clip3.min(gap);
    let residual_gap = gap - covered;
    let overhang = clip3 - covered;
    residual_gap + overhang.saturating_sub(OVERHANG_ALLOWANCE_NT)
}

/// Anchored-3' likelihood (`--anchor-three-prime`, dRNA): every read is
/// treated as 3'-complete — the direct-RNA protocol sequences from the
/// poly(A), so truncation is 5'-sided and the read's 3' end is the
/// molecule's 3' end (premise measured at 91.7% clip-inclusive on real SIRV
/// dRNA). Identical machinery to the poly(A) term (empirical gap-bin
/// distribution trained on unique reads, odds capped at MAX_ODDS) but
/// ungated and clip-aware. Do not use on cDNA/PacBio (premise fails there).
pub(crate) fn apply_anchored_probabilities(
    store: &mut InMemoryAlignmentStore,
    txps: &[TranscriptInfo],
) -> PolyaDiagnostics {
    let bins = GAP_EDGES.len();
    let mut counts = vec![0.0f64; bins];
    let mut training = 0usize;
    for (alns, _, _) in store.iter() {
        if let [aln] = alns {
            let len = txps[aln.ref_id as usize].len.get();
            counts[gap_bin(three_prime_gap_clip_aware_nt(aln, len))] += 1.0;
            training += 1;
        }
    }
    let total: f64 = counts.iter().sum();
    let density = |bin: usize| -> f64 {
        ((counts[bin] + PRIOR_MASS / bins as f64) / (total + PRIOR_MASS)) / bin_width(bin)
    };
    let mut scored = 0usize;
    let mut odds_range_sum = 0.0;
    let mut local: Vec<f64> = Vec::new();
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start < 2 {
            continue;
        }
        local.clear();
        for aln in &store.alignments[start..end] {
            let len = txps[aln.ref_id as usize].len.get();
            local.push(density(gap_bin(three_prime_gap_clip_aware_nt(aln, len))));
        }
        let max = local.iter().copied().fold(0.0, f64::max);
        if !(max > 0.0) {
            continue;
        }
        let floor = max / MAX_ODDS;
        local.iter_mut().for_each(|v| *v = v.max(floor));
        let sum: f64 = local.iter().sum();
        if !(sum > 0.0) {
            continue;
        }
        let min = local.iter().copied().fold(f64::INFINITY, f64::min);
        odds_range_sum += (max.max(floor) / min).ln();
        for (probability, value) in store.coverage_probabilities[start..end]
            .iter_mut()
            .zip(&local)
        {
            *probability *= value / sum;
        }
        let renorm: f64 = store.coverage_probabilities[start..end].iter().sum();
        if renorm > 0.0 && renorm.is_finite() {
            store.coverage_probabilities[start..end]
                .iter_mut()
                .for_each(|v| *v /= renorm);
        } else {
            let uniform = 1.0 / (end - start) as f64;
            store.coverage_probabilities[start..end].fill(uniform);
        }
        scored += 1;
    }
    let flush = if total > 0.0 {
        counts[..gap_bin(20)].iter().sum::<f64>() / total
    } else {
        0.0
    };
    let diffuse = if total > 0.0 {
        counts[gap_bin(256)..].iter().sum::<f64>() / total
    } else {
        0.0
    };
    PolyaDiagnostics {
        training_reads: training,
        scored_reads: scored,
        tail_reads: 0,
        flush_fraction: flush,
        diffuse_fraction: diffuse,
        mean_log_odds_range: odds_range_sum / scored.max(1) as f64,
    }
}

/// Shrinkage pseudo-mass for per-transcript end profiles: a transcript with
/// no unique reads scores exactly like the global model; profiles earn
/// their influence with data.
const PROFILE_PRIOR_MASS: f64 = 20.0;

/// Clip-credited 5' gap: distance from the molecule's 5' (read-start) side
/// to the candidate's annotated start, crediting the 5'-side soft clip
/// (ragged dRNA read starts are clipped by the aligner).
#[inline]
fn five_prime_gap_clip_credited_nt(aln: &AlnInfo, len: usize) -> u32 {
    let left = aln.start.saturating_sub(1);
    let right = len.saturating_sub(aln.end as usize) as u32;
    let (gap, clip) = if aln.strand == Strand::Reverse {
        (right, aln.right_clip)
    } else {
        (left, aln.left_clip)
    };
    gap.saturating_sub(clip)
}

/// Two-sided per-transcript end-profile likelihood (`--end-profile`, dRNA).
///
/// Bottom-up motivation (measured, 2026-08-19): real dRNA molecule ends are
/// strongly transcript-specific — PAS motifs sit upstream of read 3' ends
/// (54%/33% vs 2.5% background), 26% of reads reach an annotated start
/// (~26x uniform truncation), and both ends cluster at recurrent
/// per-transcript sites. Held-out validation: per-transcript profiles beat
/// the global end distribution by ~+1.1-1.3 nats/read on real samples,
/// positive in every coverage stratum.
///
/// Model: per transcript, empirical distributions over log-spaced gap bins
/// for the 5' (clip-credited) and 3' (signed clip-aware) end mismatches,
/// trained on the sample's unique reads, shrunk toward the global
/// distribution with PROFILE_PRIOR_MASS pseudo-reads. Candidates for a
/// multi-mapping read are scored by the product of their two end densities,
/// odds capped at MAX_ODDS. Unlike the global anchored model, APA and
/// unannotated ends become per-transcript signal, and the 5' side can
/// separate 3'-sharing isoforms.
pub(crate) fn apply_end_profile_probabilities(
    store: &mut InMemoryAlignmentStore,
    txps: &[TranscriptInfo],
) -> PolyaDiagnostics {
    let bins = GAP_EDGES.len();
    let nt = txps.len();
    let mut c5 = vec![0u32; nt * bins];
    let mut c3 = vec![0u32; nt * bins];
    let mut ntr = vec![0u32; nt];
    let mut g5 = vec![0.0f64; bins];
    let mut g3 = vec![0.0f64; bins];
    let mut training = 0usize;
    for (alns, _, _) in store.iter() {
        if let [aln] = alns {
            let tid = aln.ref_id as usize;
            let len = txps[tid].len.get();
            let b3 = gap_bin(three_prime_gap_clip_aware_nt(aln, len));
            let b5 = gap_bin(five_prime_gap_clip_credited_nt(aln, len));
            c3[tid * bins + b3] += 1;
            c5[tid * bins + b5] += 1;
            ntr[tid] += 1;
            g3[b3] += 1.0;
            g5[b5] += 1.0;
            training += 1;
        }
    }
    let t3: f64 = g3.iter().sum();
    let t5: f64 = g5.iter().sum();
    // global bin probability with a small uniform prior so empty bins stay finite
    let gprob = |g: &[f64], tot: f64, b: usize| -> f64 {
        (g[b] + PRIOR_MASS / bins as f64) / (tot + PRIOR_MASS)
    };
    let dens = |c: &[u32], g: &[f64], tot: f64, tid: usize, b: usize| -> f64 {
        let n = ntr[tid] as f64;
        let p = (c[tid * bins + b] as f64 + PROFILE_PRIOR_MASS * gprob(g, tot, b))
            / (n + PROFILE_PRIOR_MASS);
        p / bin_width(b)
    };
    let mut scored = 0usize;
    let mut odds_range_sum = 0.0;
    let mut local: Vec<f64> = Vec::new();
    for read_index in 0..store.len() {
        let start = store.boundaries[read_index];
        let end = store.boundaries[read_index + 1];
        if end - start < 2 {
            continue;
        }
        local.clear();
        for aln in &store.alignments[start..end] {
            let tid = aln.ref_id as usize;
            let len = txps[tid].len.get();
            let b3 = gap_bin(three_prime_gap_clip_aware_nt(aln, len));
            let b5 = gap_bin(five_prime_gap_clip_credited_nt(aln, len));
            local.push(dens(&c3, &g3, t3, tid, b3) * dens(&c5, &g5, t5, tid, b5));
        }
        let max = local.iter().copied().fold(0.0, f64::max);
        if !(max > 0.0) {
            continue;
        }
        let floor = max / MAX_ODDS;
        local.iter_mut().for_each(|v| *v = v.max(floor));
        let sum: f64 = local.iter().sum();
        if !(sum > 0.0) {
            continue;
        }
        let min = local.iter().copied().fold(f64::INFINITY, f64::min);
        odds_range_sum += (max.max(floor) / min).ln();
        for (probability, value) in store.coverage_probabilities[start..end]
            .iter_mut()
            .zip(&local)
        {
            *probability *= value / sum;
        }
        let renorm: f64 = store.coverage_probabilities[start..end].iter().sum();
        if renorm > 0.0 && renorm.is_finite() {
            store.coverage_probabilities[start..end]
                .iter_mut()
                .for_each(|v| *v /= renorm);
        } else {
            let uniform = 1.0 / (end - start) as f64;
            store.coverage_probabilities[start..end].fill(uniform);
        }
        scored += 1;
    }
    let flush = if t3 > 0.0 {
        g3[..gap_bin(20)].iter().sum::<f64>() / t3
    } else {
        0.0
    };
    let diffuse = if t3 > 0.0 {
        g3[gap_bin(256)..].iter().sum::<f64>() / t3
    } else {
        0.0
    };
    PolyaDiagnostics {
        training_reads: training,
        scored_reads: scored,
        tail_reads: 0,
        flush_fraction: flush,
        diffuse_fraction: diffuse,
        mean_log_odds_range: odds_range_sum / scored.max(1) as f64,
    }
}
