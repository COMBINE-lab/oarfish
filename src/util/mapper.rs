//! Mapping-backend abstraction.
//!
//! oarfish can map raw reads with one of two backends:
//!
//! * the default [`minimap2`]/`minimap2-rs` backend (matches the `dev` branch), or
//! * the pure-Rust [`rammap`](https://github.com/jwanglab/rammap) backend, enabled
//!   with the `rammap` Cargo feature.
//!
//! Both backends are exposed to the rest of the crate through a small, uniform
//! surface so the hot mapping loops in [`crate::bulk`] don't need to care which
//! one is in use:
//!
//! * [`Mapper`] — the built aligner type (cheaply `clone`-able per worker thread;
//!   `minimap2::Aligner<Built>` clones its `Arc`-shared index, and the rammap
//!   variant is an [`Arc`] around a single shared aligner).
//! * [`OMapping`] — the per-alignment record type yielded for a read. It always
//!   implements [`AlnRecordLike`](crate::util::oarfish_types::AlnRecordLike)
//!   (transcriptome filter) and is convertible to a bramble `GenomicAlignment`
//!   via [`crate::util::projection::mapping_to_genomic_alignment`] (genome path).
//! * [`map_read`] — map one read, returning `anyhow::Result<Vec<OMapping>>`.

// ---------------------------------------------------------------------------
// Default backend: minimap2-rs
// ---------------------------------------------------------------------------

#[cfg(not(feature = "rammap"))]
mod backend {
    /// The built aligner. `minimap2::Aligner<Built>` is `Clone` (it shares the
    /// underlying index behind an `Arc`), so each worker thread can hold its own
    /// handle with per-thread mapping buffers.
    pub(crate) type Mapper = minimap2::Aligner<minimap2::Built>;

    /// The per-alignment record yielded for a read.
    pub(crate) type OMapping = minimap2::Mapping;

    /// Map a single read against the index, returning all of its alignments.
    #[inline]
    pub(crate) fn map_read(
        aligner: &Mapper,
        name: &[u8],
        seq: &[u8],
    ) -> anyhow::Result<Vec<OMapping>> {
        // (seq, cigar = true, md = false, max_frag_len = None, extra_flags = None, name)
        aligner
            .map(seq, true, false, None, None, Some(name))
            .map_err(anyhow::Error::msg)
    }

    /// The minimap2 alignment score (`AS`), or 0 if absent.
    #[inline]
    pub(crate) fn alignment_score(m: &OMapping) -> i32 {
        m.alignment
            .as_ref()
            .and_then(|a| a.alignment_score)
            .unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// Alternative backend: rammap (pure Rust)
// ---------------------------------------------------------------------------

#[cfg(feature = "rammap")]
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
    /// `Ok(..)` so the call sites are identical to the minimap2 backend.
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

// `OMapping` is only named by other modules in the rammap cfg path.
#[allow(unused_imports)]
pub(crate) use backend::{Mapper, OMapping, alignment_score, map_read};
