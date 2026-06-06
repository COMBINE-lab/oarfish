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
        Ok(res
            .mappings
            .into_iter()
            .map(|m| RammapMapping {
                m,
                query_len: seq.len(),
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

pub(crate) use backend::{Mapper, OMapping, alignment_score, map_read};
