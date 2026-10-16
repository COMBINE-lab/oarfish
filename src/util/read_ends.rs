//! Backend-neutral evidence at a read's ends, and the score adjustments built on it.
//!
//! Mapping backends differ in what they can say about a read's ends: a SAM record carries a
//! CIGAR (terminal clips) and, on the record that holds SEQ, the read's bases (its poly(A)
//! tail); rammap reports query spans and is handed the read's sequence; a mapper that reports
//! spans between its outermost anchors cannot say how much of the read lies beyond them. Each
//! backend therefore reports an [`EndEvidence`] per alignment, with `None` wherever it cannot
//! tell, and the models here consume only that. Anything that cannot be established is left
//! unadjusted rather than guessed.

/// What an alignment says about the read's ends, in reference orientation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EndEvidence {
    /// Read bases beyond the alignment at the reference-left and reference-right ends.
    pub clip_left: Option<u32>,
    pub clip_right: Option<u32>,
    /// Length of the read's poly(A) tail, which in reference orientation sits at the right end
    /// (`Some(0)`: the read has none; `None`: this record cannot tell).
    pub tail: Option<u32>,
}

impl EndEvidence {
    pub const UNKNOWN: EndEvidence = EndEvidence {
        clip_left: None,
        clip_right: None,
        tail: None,
    };
}

/// Shortest A-run counted as a poly(A) tail.
const MIN_TAIL: u32 = 8;

/// Length of the run of `base` at the start of `seq` (pass the read backward to measure a run
/// at its end), tolerating up to one other base per ten; runs shorter than [`MIN_TAIL`] count as
/// no tail.
pub fn terminal_run(seq: impl Iterator<Item = u8>, base: u8) -> u32 {
    let (mut n, mut other, mut best) = (0_u32, 0_u32, 0_u32);
    for b in seq {
        n += 1;
        if b.to_ascii_uppercase() == base {
            best = n;
        } else {
            other += 1;
            if other * 10 > n.max(10) {
                break;
            }
        }
    }
    if best >= MIN_TAIL { best } else { 0 }
}

/// Removes the score credit a read's poly(A) tail earns where it aligns to a transcript, which
/// happens when the reference sequence itself ends in an A-run (an annotated tail). Tail bases
/// are not evidence for any transcript; left in, they let a transcript with an annotated tail
/// take reads from siblings without one.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct TailNeutral {
    /// The backend's score per matched base.
    pub match_score: i32,
}

impl TailNeutral {
    /// The read's tail length and the strand of the record it was measured on, from the first
    /// alignment whose backend can tell; `None` if none can or the read has no tail.
    pub fn read_tail(evidence: impl Iterator<Item = (EndEvidence, bool)>) -> Option<(u32, bool)> {
        evidence
            .filter_map(|(ev, rc)| ev.tail.map(|n| (n, rc)))
            .next()
            .filter(|&(n, _)| n > 0)
    }

    /// `score` without the credit for tail bases this alignment covers: the tail bases not
    /// clipped at its right end. Alignments on the other strand from the one the tail was
    /// measured on, or whose right clip is unknown, are left unchanged.
    pub fn adjust(&self, score: i32, tail: Option<(u32, bool)>, ev: &EndEvidence, rc: bool) -> i32 {
        match (tail, ev.clip_right) {
            (Some((n, tail_rc)), Some(clip)) if tail_rc == rc && score > i32::MIN => {
                score - self.match_score * n.saturating_sub(clip) as i32
            }
            _ => score,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tail detection: an A-run at the end counts (one stray base per ten tolerated); short runs
    /// and A-poor ends do not.
    #[test]
    fn terminal_run_detects_tails() {
        let run = |s: &str| terminal_run(s.bytes().rev(), b'A');
        assert_eq!(run("CGTACGAAAAAAAAAAAA"), 12);
        assert_eq!(run("CGTACGAAAAAGAAAAAAAAA"), 15);
        assert_eq!(run("CGTACGTAAAAA"), 0);
        assert_eq!(run("AAAAAAAAAACGTACGTCGT"), 0);
        assert_eq!(terminal_run("TTTTTTTTTTGCA".bytes(), b'T'), 10);
    }

    /// A transcript that aligns the tail loses its credit; one that clips it is unchanged; so do
    /// alignments whose clip is unknown or on the other strand.
    #[test]
    fn tail_credit_is_removed_only_where_the_tail_aligned() {
        let tn = TailNeutral { match_score: 2 };
        let tail = Some((12, false));
        let ev = |clip: Option<u32>| EndEvidence {
            clip_left: Some(0),
            clip_right: clip,
            tail: None,
        };
        assert_eq!(tn.adjust(100, tail, &ev(Some(0)), false), 76);
        assert_eq!(tn.adjust(100, tail, &ev(Some(12)), false), 100);
        assert_eq!(tn.adjust(100, tail, &ev(Some(5)), false), 86);
        assert_eq!(tn.adjust(100, tail, &ev(None), false), 100);
        assert_eq!(tn.adjust(100, tail, &ev(Some(0)), true), 100);
        assert_eq!(tn.adjust(100, None, &ev(Some(0)), false), 100);
    }
}
