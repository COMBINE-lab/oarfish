//! Mapping-backend abstraction.
//!
//! oarfish maps raw reads with the pure-Rust
//! [`rammap`](https://github.com/jwanglab/rammap) backend.
//!
//! The backend is exposed to the rest of the crate through a small, uniform
//! surface so the hot mapping loops in [`crate::bulk`] don't need to care about
//! its internals:
//!
//! * [`Mapper`] — the built aligner type. rammap's `Aligner` is `Send + Sync`,
//!   so a single shared instance is held behind an [`Arc`] across worker threads.
//! * [`OMapping`] — the per-alignment record type yielded for a read. It
//!   implements [`AlnRecordLike`](crate::util::oarfish_types::AlnRecordLike)
//!   (transcriptome filter) and is convertible to a bramble `GenomicAlignment`
//!   via [`crate::util::projection::mapping_to_genomic_alignment`] (genome path).
//! * [`map_read`] — map one read, returning `anyhow::Result<Vec<OMapping>>`.

// ---------------------------------------------------------------------------
// Mapping backend: rammap (pure Rust)
// ---------------------------------------------------------------------------

mod backend {
    use crate::util::oarfish_types::AlnRecordLike;
    use noodles_sam::header::Header;
    use std::sync::Arc;

    /// The built aligner. rammap's `Aligner` is not `Clone` but is `Send + Sync`,
    /// so we share a single instance across worker threads behind an `Arc`
    /// (each `map_seq` call allocates its own lightweight per-call buffers).
    pub(crate) type Mapper = Arc<rammap::api::Aligner>;

    /// A rammap alignment plus the query length (rammap's `Mapping` does not carry
    /// the query length, but the transcriptome filter needs it via
    /// [`AlnRecordLike::opt_sequence_len`]).
    #[derive(Debug, Clone)]
    pub(crate) struct RammapMapping {
        pub(crate) m: rammap::api::Mapping,
        pub(crate) query_len: usize,
        /// poly(A) run at the reference-right end of the read in this mapping's orientation.
        pub(crate) tail: u32,
    }

    /// The per-alignment record yielded for a read.
    pub(crate) type OMapping = RammapMapping;

    /// Map a single read against the index, returning all of its alignments.
    ///
    /// rammap's `map_seq` is infallible and operates on a single read, which fits
    /// oarfish's existing per-thread mapping loop directly. We wrap the result in
    /// `Ok(..)` so the call sites are uniform across the mapping surface.
    #[inline]
    pub(crate) fn map_read(
        aligner: &Mapper,
        name: &[u8],
        seq: &[u8],
    ) -> anyhow::Result<Vec<OMapping>> {
        let qname = std::str::from_utf8(name).unwrap_or("read");
        let res = aligner.map_seq(qname, seq);
        // In reference orientation the tail is an A-run at the right end: the read's 3' A-run
        // when it maps forward, its 5' T-run when it maps reverse.
        let tail_fw = crate::util::read_ends::terminal_run(seq.iter().rev().copied(), b'A');
        let tail_rc = crate::util::read_ends::terminal_run(seq.iter().copied(), b'T');
        Ok(res
            .mappings
            .into_iter()
            .map(|m| {
                let tail = if m.strand == rammap::api::Strand::Reverse {
                    tail_rc
                } else {
                    tail_fw
                };
                RammapMapping {
                    m,
                    query_len: seq.len(),
                    tail,
                }
            })
            .collect())
    }

    impl AlnRecordLike for RammapMapping {
        fn opt_sequence_len(&self) -> Option<usize> {
            Some(self.query_len)
        }

        fn is_reverse_complemented(&self) -> bool {
            self.m.strand == rammap::api::Strand::Reverse
        }

        fn is_unmapped(&self) -> bool {
            // rammap only returns mapped alignments.
            false
        }

        fn ref_id(&self, _header: &Header) -> anyhow::Result<usize> {
            // rammap's `target_id` is the 0-based index into the aligner's
            // reference list, which is the same order the SAM header is built in
            // (see `crate::util::aligner`), so it equals the header reference id.
            Ok(self.m.target_id)
        }

        fn aln_span(&self) -> Option<usize> {
            // reference-consuming span of the alignment.
            Some(self.m.target_end - self.m.target_start)
        }

        fn aln_score(&self) -> Option<i64> {
            Some(self.m.score as i64)
        }

        fn aln_start(&self) -> u32 {
            self.m.target_start as u32
        }

        fn aln_end(&self) -> u32 {
            self.m.target_end as u32
        }

        fn is_supp(&self) -> bool {
            self.m.is_supplementary
        }

        fn end_evidence(&self) -> crate::util::read_ends::EndEvidence {
            // Query spans are on the read as given; in reference orientation a reverse mapping
            // swaps which end each clip is at.
            let (qs, qe) = (self.m.query_start as u32, self.m.query_end as u32);
            let ql = self.query_len as u32;
            let (clip_left, clip_right) = if self.m.strand == rammap::api::Strand::Reverse {
                (ql.saturating_sub(qe), qs)
            } else {
                (qs, ql.saturating_sub(qe))
            };
            crate::util::read_ends::EndEvidence {
                clip_left: Some(clip_left),
                clip_right: Some(clip_right),
                tail: Some(self.tail),
            }
        }

        fn name(&self) -> Option<String> {
            None
        }
    }

    /// The rammap alignment score (`AS`).
    #[inline]
    pub(crate) fn alignment_score(m: &OMapping) -> i32 {
        m.m.score
    }
}

use crate::util::oarfish_types::AlnRecordLike;
pub(crate) use backend::{OMapping, alignment_score, map_read};

/// The rammap aligner, as used by the genome-projection path.
pub(crate) type RammapMapper = std::sync::Arc<rammap::api::Aligner>;

/// Transcriptome read-mode mapping backend: rammap (default) or spearfish.
#[derive(Clone)]
pub(crate) enum Mapper {
    Rammap(RammapMapper),
    Spearfish(std::sync::Arc<spearfish::Mapper>),
}

/// spearfish scores are real-valued; oarfish's filter takes integers. Scores are scaled by
/// this factor.
pub(crate) const SPEARFISH_SCORE_SCALE: f32 = 10.0;

/// Default score→probability denominator for spearfish hits (scaled scores): a hit Δ below
/// the read's best (spearfish units) weighs exp(-Δ). Swept on simulated reads (ONT 95%, ONT 90%,
/// HiFi) and spike-ins, as `--score-prob-denom` was for alignments (5 -> 3): exp(-Δ) was best or
/// tied on every set; exp(-Δ / 2), the previous default, spread reads over near-tied isoforms.
pub(crate) const SPEARFISH_SCORE_DENOM: f32 = SPEARFISH_SCORE_SCALE;

/// A spearfish hit plus the read length (for the aligned-fraction filter).
#[derive(Debug, Clone)]
pub(crate) struct SpearfishMapping {
    pub(crate) hit: spearfish::Hit,
    pub(crate) query_len: usize,
}

/// One transcriptome alignment from either backend.
#[derive(Debug, Clone)]
pub(crate) enum TxMapping {
    Rammap(OMapping),
    Spearfish(SpearfishMapping),
}

/// Map one read in transcriptome read mode with whichever backend is configured.
#[inline]
pub(crate) fn map_read_txome(
    mapper: &Mapper,
    name: &[u8],
    seq: &[u8],
) -> anyhow::Result<Vec<TxMapping>> {
    match mapper {
        Mapper::Rammap(a) => Ok(map_read(a, name, seq)?
            .into_iter()
            .map(TxMapping::Rammap)
            .collect()),
        Mapper::Spearfish(m) => {
            let r = m.map(seq);
            Ok(r.hits
                .into_iter()
                .map(|hit| {
                    TxMapping::Spearfish(SpearfishMapping {
                        hit,
                        query_len: seq.len(),
                    })
                })
                .collect())
        }
    }
}

impl AlnRecordLike for TxMapping {
    fn opt_sequence_len(&self) -> Option<usize> {
        match self {
            TxMapping::Rammap(m) => m.opt_sequence_len(),
            TxMapping::Spearfish(m) => Some(m.query_len),
        }
    }
    fn is_reverse_complemented(&self) -> bool {
        match self {
            TxMapping::Rammap(m) => m.is_reverse_complemented(),
            TxMapping::Spearfish(m) => m.hit.reverse,
        }
    }
    fn is_unmapped(&self) -> bool {
        match self {
            TxMapping::Rammap(m) => m.is_unmapped(),
            TxMapping::Spearfish(_) => false,
        }
    }
    fn ref_id(&self, header: &noodles_sam::header::Header) -> anyhow::Result<usize> {
        match self {
            TxMapping::Rammap(m) => m.ref_id(header),
            // spearfish transcript ids follow index order, which is the header order.
            TxMapping::Spearfish(m) => Ok(m.hit.tx as usize),
        }
    }
    fn aln_span(&self) -> Option<usize> {
        match self {
            TxMapping::Rammap(m) => m.aln_span(),
            TxMapping::Spearfish(m) => Some((m.hit.target_end - m.hit.target_start) as usize),
        }
    }
    fn aln_score(&self) -> Option<i64> {
        match self {
            TxMapping::Rammap(m) => m.aln_score(),
            TxMapping::Spearfish(m) => Some((m.hit.score * SPEARFISH_SCORE_SCALE).round() as i64),
        }
    }
    fn aln_start(&self) -> u32 {
        match self {
            TxMapping::Rammap(m) => m.aln_start(),
            TxMapping::Spearfish(m) => m.hit.target_start,
        }
    }
    fn aln_end(&self) -> u32 {
        match self {
            TxMapping::Rammap(m) => m.aln_end(),
            TxMapping::Spearfish(m) => m.hit.target_end,
        }
    }
    fn is_supp(&self) -> bool {
        match self {
            TxMapping::Rammap(m) => m.is_supp(),
            TxMapping::Spearfish(_) => false,
        }
    }
    fn name(&self) -> Option<String> {
        None
    }
}
