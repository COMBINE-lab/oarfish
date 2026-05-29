//! Bridge between bramble-rs (genome→transcriptome projection) and oarfish's
//! quantification machinery.
//!
//! In oarfish's "genome-space" modes the *input* references are chromosomes,
//! but the quantification *targets* are annotated transcripts. This module
//! builds the transcript set (header + [`TranscriptInfo`]) from a bramble
//! [`G2TTree`] rather than from the input BAM/aligner header, and converts
//! bramble [`ProjectedAlignment`]s into oarfish [`ProjectedAlnRecord`]s that
//! feed the existing alignment filters and EM.
//!
//! The key invariant: bramble assigns each transcript a dense 0-based id in
//! annotation/build order, and we build the SAM header in that same order, so
//! `ProjectedAlignment::transcript_id` directly indexes `txps`/`txps_name` and
//! becomes `AlnInfo::ref_id` with no remapping.

use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, anyhow};
use tracing::{info, warn};

use noodles_sam::header::Header;
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;

use bramble_rs::ProjectedAlignment;
use bramble_rs::ProjectionConfig;
use bramble_rs::annotation::load_transcripts;
use bramble_rs::fasta::FastaDb;
use bramble_rs::g2t::{G2TTree, build_g2t_from_refnames};

use crate::prog_opts::Args;
use crate::util::oarfish_types::{ProjectedAlnRecord, TranscriptInfo};

/// Load a GTF/GFF annotation and build the genome→transcriptome index.
///
/// `refnames[i]` is the chromosome whose 0-based reference id is `i`; it must
/// match the `ref_id` placed on every `GenomicAlignment` handed to
/// [`bramble_rs::project_group`] (i.e. the input BAM header order, or the
/// minimap2 target order). When `genome_fasta` is provided, exon sequences are
/// loaded so bramble can perform soft-clip rescue.
pub fn load_g2t(
    annotation: &Path,
    refnames: &[String],
    genome_fasta: Option<&Path>,
) -> anyhow::Result<G2TTree> {
    info!("loading annotation from {}", annotation.display());
    let transcripts = load_transcripts(annotation)
        .with_context(|| format!("failed to load annotation {}", annotation.display()))?;
    info!("loaded {} transcripts from annotation", transcripts.len());

    // The FastaDb only needs to outlive `build_g2t_from_refnames`; bramble
    // copies the per-exon sequences it needs into the index.
    let fasta = match genome_fasta {
        Some(p) => {
            info!("loading genome FASTA for soft-clip rescue from {}", p.display());
            Some(FastaDb::load(p).with_context(|| format!("failed to load genome FASTA {}", p.display()))?)
        }
        None => None,
    };

    let g2t = build_g2t_from_refnames(&transcripts, refnames, fasta.as_ref())
        .context("failed to build genome-to-transcriptome index")?;
    info!("built g2t index over {} transcripts", g2t.num_transcripts());
    Ok(g2t)
}

/// Build a transcriptome SAM header plus the parallel `TranscriptInfo` and name
/// vectors from a bramble [`G2TTree`].
///
/// The header references are the transcripts, in bramble's dense `tid` order,
/// so a returned `txps[tid]` / `names[tid]` lines up with
/// `ProjectedAlignment::transcript_id`.
pub fn build_transcriptome_header_and_info(
    g2t: &G2TTree,
    args: &Args,
) -> anyhow::Result<(Header, Vec<TranscriptInfo>, Vec<String>)> {
    let n = g2t.num_transcripts();
    let mut builder = Header::builder();
    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(n);
    let mut names: Vec<String> = Vec::with_capacity(n);

    let mut n_zero_len = 0usize;
    for tid in 0..n {
        let tid = tid as u32;
        let name = g2t
            .transcript_name(tid)
            .ok_or_else(|| anyhow!("g2t index missing a name for transcript id {tid}"))?
            .to_string();
        let raw_len = g2t
            .transcript_len(tid)
            .ok_or_else(|| anyhow!("g2t index missing a length for transcript id {tid}"))?;
        // A transcript with no usable exons would have length 0; substitute 1 so
        // the header/`TranscriptInfo` stay index-aligned with bramble's tids. No
        // alignment will project to such a transcript.
        if raw_len == 0 {
            n_zero_len += 1;
        }
        let len = NonZeroUsize::new(raw_len.max(1) as usize).expect("len >= 1");

        builder = builder.add_reference_sequence(
            name.clone(),
            HeaderMap::<header_val::map::ReferenceSequence>::new(len),
        );

        let info = if args.model_coverage {
            TranscriptInfo::with_len_and_bin_width(len, args.bin_width)
        } else {
            TranscriptInfo::with_len(len)
        };
        txps.push(info);
        names.push(name);
    }

    builder = builder.add_program(
        "bramble",
        HeaderMap::<header_val::map::Program>::default(),
    );

    if n_zero_len > 0 {
        warn!(
            "{} transcript(s) in the annotation had zero exonic length; they were retained with length 1.",
            n_zero_len
        );
    }
    info!(
        "constructed transcriptome header with {} transcripts.",
        names.len()
    );

    Ok((builder.build(), txps, names))
}

/// Convert bramble's per-alignment projection results into oarfish's neutral
/// [`ProjectedAlnRecord`] hand-off type (consumed by
/// [`crate::util::oarfish_types::AlignmentFilters::filter_projected`]).
pub fn projected_to_records(projected: &[ProjectedAlignment]) -> Vec<ProjectedAlnRecord> {
    projected
        .iter()
        .map(|p| ProjectedAlnRecord {
            ref_id: p.transcript_id,
            start: p.transcript_start,
            end: p.transcript_end,
            aligned_len: p.aligned_len,
            query_aligned_len: p.query_aligned_len,
            is_reverse: p.is_reverse,
            similarity: p.similarity_score,
        })
        .collect()
}

/// The projection configuration to use for the run. oarfish is a long-read
/// quantifier, so genome-space projection always uses long-read tolerances;
/// soft-clip rescue is enabled when a genome FASTA is available.
pub fn projection_config(args: &Args) -> ProjectionConfig {
    ProjectionConfig {
        long_reads: true,
        use_fasta: args.genome_fasta.is_some(),
    }
}
