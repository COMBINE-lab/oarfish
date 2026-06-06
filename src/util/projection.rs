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

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use anyhow::{Context, anyhow};
use tracing::{info, warn};

use noodles_sam::header::Header;
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;

use bramble_rs::GenomicAlignment;
use bramble_rs::ProjectedAlignment;
use bramble_rs::ProjectionConfig;
use bramble_rs::annotation::{Transcript, load_transcripts};
use bramble_rs::fasta::FastaDb;
use bramble_rs::g2t::{G2TTree, build_g2t_from_refnames};

use crate::prog_opts::Args;
use crate::util::oarfish_types::{ProjectedAlnRecord, TranscriptInfo};

/// Load a GTF/GFF annotation and build the genome→transcriptome index.
///
/// `refnames[i]` is the chromosome whose 0-based reference id is `i`; it must
/// match the `ref_id` placed on every `GenomicAlignment` handed to
/// [`bramble_rs::project_group`] (i.e. the input BAM header order, or the
/// aligner target order). When `genome_fasta` is provided, exon sequences are
/// loaded so bramble can perform soft-clip rescue.
pub fn load_g2t(
    annotation: &Path,
    refnames: &[String],
    rescue_fasta: Option<FastaDb>,
) -> anyhow::Result<G2TTree> {
    info!("loading annotation from {}", annotation.display());
    let transcripts = load_transcripts(annotation)
        .with_context(|| format!("failed to load annotation {}", annotation.display()))?;
    info!("loaded {} transcripts from annotation", transcripts.len());
    build_g2t_from_transcripts(&transcripts, refnames, rescue_fasta)
}

/// Build the genome→transcriptome index from already-parsed transcripts. Lets
/// callers (e.g. genome-read mode) parse the annotation once and reuse the
/// transcripts for both junction extraction and index building.
///
/// `rescue_fasta` is the (already-resolved) reference sequence for bramble's
/// soft-clip rescue, sourced by the caller either from a FASTA file or from the
/// loaded aligner index; `None` disables rescue. bramble copies the per-exon
/// sequences it needs into the index, so the `FastaDb` is dropped on return.
pub fn build_g2t_from_transcripts(
    transcripts: &[Transcript],
    refnames: &[String],
    rescue_fasta: Option<FastaDb>,
) -> anyhow::Result<G2TTree> {
    let g2t = build_g2t_from_refnames(transcripts, refnames, rescue_fasta.as_ref())
        .context("failed to build genome-to-transcriptome index")?;
    info!("built g2t index over {} transcripts", g2t.num_transcripts());
    Ok(g2t)
}

/// Write a BED12 of transcript models derived from `transcripts`, suitable for
/// the aligner's `--junc-bed`-equivalent junction loader. The aligner derives
/// the splice junctions from the exon-block structure and uses them to
/// guide/score spliced alignment, which non-trivially improves alignment
/// accuracy.
///
/// Only multi-exon transcripts (which actually define junctions) are written.
/// Returns the number of transcript models written. Coordinates: bramble `Exon`
/// is 1-based `start` with `end == 1-based-inclusive-end + 1`, so the 0-based
/// half-open BED exon is `[start-1, end-1)` with size `end-start`.
pub fn write_annotation_junction_bed(
    transcripts: &[Transcript],
    path: &Path,
) -> anyhow::Result<usize> {
    let f = File::create(path)
        .with_context(|| format!("failed to create junction BED {}", path.display()))?;
    let mut w = BufWriter::new(f);
    let mut n = 0usize;
    for tx in transcripts {
        if tx.exons.len() < 2 {
            continue; // single-exon transcripts contribute no junctions
        }
        let mut exons: Vec<(u32, u32)> = tx
            .exons
            .iter()
            .map(|e| (e.start.saturating_sub(1), e.end.saturating_sub(1)))
            .collect();
        exons.sort_by_key(|e| e.0);
        let chrom_start = exons[0].0;
        let chrom_end = exons.last().unwrap().1;
        let strand = match tx.strand {
            '+' => '+',
            '-' => '-',
            _ => '.',
        };
        let mut sizes = String::new();
        let mut starts = String::new();
        for (s, e) in &exons {
            sizes.push_str(&(e - s).to_string());
            sizes.push(',');
            starts.push_str(&(s - chrom_start).to_string());
            starts.push(',');
        }
        // BED12: chrom start end name score strand thickStart thickEnd rgb blockCount blockSizes blockStarts
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tx.seqname,
            chrom_start,
            chrom_end,
            tx.id,
            1000,
            strand,
            chrom_start,
            chrom_end,
            0,
            exons.len(),
            sizes,
            starts,
        )?;
        n += 1;
    }
    w.flush()?;
    Ok(n)
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
///
/// `src_scores[i]` is the alignment score of the i-th source genomic alignment
/// (the slice passed to `project_group`); each `ProjectedAlignment` carries its
/// `input_index` into that slice, so projections from the same genomic locus
/// share a score. An empty `src_scores` yields `aln_score = 0` for all.
pub fn projected_to_records(
    projected: &[ProjectedAlignment],
    src_scores: &[i32],
) -> Vec<ProjectedAlnRecord> {
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
            aln_score: src_scores.get(p.input_index).copied().unwrap_or(0),
        })
        .collect()
}

/// The projection configuration to use for the run. oarfish is a long-read
/// quantifier, so genome-space projection always uses long-read tolerances;
/// soft-clip rescue is enabled when a genome FASTA is available.
/// Resolve the genome FASTA used for bramble's soft-clip rescue, honoring the
/// rescue-on-by-default policy and `--no-rescue`. Priority: explicit
/// `--genome-fasta`; else `--genome` itself when it is a FASTA (genome read
/// mode). Returns `None` when rescue is disabled, or when no sequence source is
/// available — e.g. `--genome` is a prebuilt index and no `--genome-fasta` was
/// given (index-sourced rescue is a separate path).
pub fn rescue_fasta_path(args: &Args) -> Option<std::path::PathBuf> {
    if args.no_rescue {
        return None;
    }
    if let Some(p) = &args.genome_fasta {
        return Some(p.clone());
    }
    if let Some(g) = &args.genome
        && crate::util::file_utils::is_fasta(g).unwrap_or(false)
    {
        return Some(g.clone());
    }
    None
}

/// Load the soft-clip-rescue reference from a FASTA path, if one is resolved by
/// [`rescue_fasta_path`] (explicit `--genome-fasta`, or a FASTA `--genome`).
/// Returns `None` when rescue is disabled or no FASTA source applies (e.g.
/// `--genome` is a prebuilt index — sequences are then sourced from the index).
pub fn load_rescue_fasta(args: &Args) -> anyhow::Result<Option<FastaDb>> {
    match rescue_fasta_path(args) {
        Some(p) => {
            info!("loading genome FASTA for soft-clip rescue from {}", p.display());
            Ok(Some(FastaDb::load(&p).with_context(|| {
                format!("failed to load genome FASTA {}", p.display())
            })?))
        }
        None => Ok(None),
    }
}

/// Long-read projection config. `use_fasta` reflects whether a rescue reference
/// was actually resolved (from a FASTA or the aligner index) — soft-clip rescue
/// is on by default and disabled by `--no-rescue` or the absence of any source.
pub fn projection_config(args: &Args, use_fasta: bool) -> ProjectionConfig {
    ProjectionConfig {
        long_reads: true,
        use_fasta,
        junc_miss_discount: args.junc_miss_discount,
    }
}

/// Reverse-complement a nucleotide sequence into a fresh uppercase buffer.
///
/// bramble's soft-clip rescue expects the read sequence in *forward-reference*
/// orientation — the orientation in which a BAM stores `SEQ` (reverse-complemented
/// relative to the original read when the alignment is on the reverse strand), and
/// the orientation the CIGAR is written in. The genome-BAM path gets this for free
/// from the BAM. In genome *read* mode we hold the read in its original FASTQ
/// orientation, so a reverse-strand mapping must be reverse-complemented before it
/// is handed to bramble; otherwise the clip-rescue slices the wrong end of the read.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'U' | b'u' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Convert a rammap spliced genome alignment (wrapped as `OMapping`) into a
/// bramble [`GenomicAlignment`]. rammap's `target_id` is the 0-based index into
/// the aligner's reference list (the same order the `G2TTree` is built in; see
/// [`crate::util::aligner::get_genome_aligner_from_args`]), so it is used
/// directly as the bramble `ref_id`. The query sequence is not attached here
/// (soft-clip rescue is handled by the caller); `read_len` carries the query
/// length needed for NH / aligned-fraction. Returns `None` for mappings lacking
/// a CIGAR (the projection needs one).
pub fn mapping_to_genomic_alignment(
    m: &crate::util::mapper::OMapping,
    query_name: &str,
    read_len: usize,
) -> Option<GenomicAlignment> {
    let cigar_ops = m.m.cigar_ops.as_ref()?;
    if cigar_ops.is_empty() {
        return None;
    }
    // rammap CigarOp uses BAM op codes (0=M,1=I,2=D,3=N,4=S,5=H,...), matching
    // the `(len, op)` representation bramble expects. BUT rammap's `cigar_ops`
    // cover only the *aligned* span (M/N/I/D); the soft-clips are implicit in
    // `query_start`/`query_end` (rammap only materializes them when writing SAM).
    // bramble's soft-clip rescue keys on explicit leading/trailing `S` ops, so we
    // reconstruct them here in reference-forward orientation, otherwise rescue
    // never fires (and isoform accuracy silently drops to the no-rescue baseline).
    let is_reverse = matches!(m.m.strand, rammap::api::Strand::Reverse);
    let clip5 = if is_reverse {
        read_len.saturating_sub(m.m.query_end)
    } else {
        m.m.query_start
    };
    let clip3 = if is_reverse {
        m.m.query_start
    } else {
        read_len.saturating_sub(m.m.query_end)
    };
    let mut cigar: Vec<(u32, u8)> = Vec::with_capacity(cigar_ops.len() + 2);
    if clip5 > 0 {
        cigar.push((clip5 as u32, 4)); // 4 = S (soft clip)
    }
    cigar.extend(cigar_ops.iter().map(|c| (c.len, c.op)));
    if clip3 > 0 {
        cigar.push((clip3 as u32, 4));
    }
    let ts_strand = m.m.trans_strand.map(|s| match s {
        rammap::api::Strand::Forward => '+',
        rammap::api::Strand::Reverse => '-',
    });
    Some(GenomicAlignment {
        query_name: query_name.to_string(),
        ref_id: m.m.target_id as i32,
        ref_start: (m.m.target_start as i64) + 1, // rammap target_start is 0-based; SAM POS is 1-based
        is_reverse,
        cigar,
        sequence: None,
        is_paired: false,
        is_first_in_pair: false,
        xs_strand: None,
        ts_strand,
        hit_index: 0,
        mate_ref_id: None,
        mate_ref_start: None,
        mate_is_unmapped: false,
        read_len,
    })
}
