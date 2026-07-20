use clap::Parser;
use noodles_sam::Header;
use std::num::NonZeroUsize;
use util::oarfish_types::NamedDigestVec;

// Genome-projection mapping churns many small allocations across worker threads;
// glibc malloc fragments/retains that into large peak RSS. mimalloc keeps the
// peak down and returns memory to the OS. (Disable with `--no-default-features`.)
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use num_format::{Locale, ToFormattedString};
use std::{fs::File, io};
use tracing::info;
use tracing_subscriber::{EnvFilter, filter::LevelFilter, fmt, prelude::*};

use noodles_bam as bam;
use noodles_bgzf as bgzf;

mod alignment_parser;
mod bootstrap;
mod bulk;
mod em;
mod em_accel;
mod prog_opts;
mod single_cell;
mod util;

use crate::prog_opts::{Args, FilterGroup};
use crate::util::aligner::{get_aligner_from_args, get_aligner_from_fastas};
use crate::util::digest_utils;
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{AlignmentFilters, TranscriptInfo};
use crate::util::projection;
use crate::util::{
    binomial_probability::binomial_continuous_prob, kde_utils, logistic_probability::logistic_prob,
};

fn get_txp_info_from_header(
    header: &Header,
    args: &Args,
) -> anyhow::Result<(Vec<TranscriptInfo>, Vec<String>)> {
    let num_ref_seqs = header.reference_sequences().len();
    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(num_ref_seqs);
    let mut txps_name: Vec<String> = Vec::with_capacity(num_ref_seqs);

    // loop over the transcripts in the header and fill in the relevant
    // information here.
    if args.model_coverage {
        for (rseq, rmap) in header.reference_sequences().iter() {
            txps.push(TranscriptInfo::with_len_and_bin_width(
                rmap.length(),
                args.bin_width,
            ));
            txps_name.push(rseq.to_string());
        }
    } else {
        for (rseq, rmap) in header.reference_sequences().iter() {
            txps.push(TranscriptInfo::with_len(rmap.length()));
            txps_name.push(rseq.to_string());
        }
    }
    info!(
        "parsed reference information for {} transcripts.",
        txps.len().to_formatted_string(&Locale::en)
    );
    Ok((txps, txps_name))
}

fn get_filter_opts(args: &Args) -> anyhow::Result<AlignmentFilters> {
    // `--score-prob-denom` tunes the transcriptome-mode score→probability
    // conversion. In genome (projection) mode a read's projected alignments are
    // weighted by bramble's exonic-coverage similarity (see `--projected-prob-source`),
    // not that conversion, so the flag has no effect there — reject it explicitly
    // rather than silently ignoring it.
    if args.score_prob_denom.is_some()
        && (args.genome.is_some() || args.genome_alignments.is_some())
    {
        anyhow::bail!(
            "--score-prob-denom does not apply in genome (projection) mode: a read's \
             projected alignments are weighted by bramble's exonic-coverage similarity \
             (see --projected-prob-source), not the score→probability conversion that \
             --score-prob-denom controls. Please omit --score-prob-denom in genome mode."
        );
    }

    // set all of the filter options that the user
    // wants to apply.
    match args.filter_group {
        Some(FilterGroup::NoFilters) => {
            info!("disabling alignment filters.");
            // override individual parameters if the user passed them in explicitly
            let fpc = args
                .five_prime_clip
                .provided_or_u32("overriding 5' clip with user-provided value", u32::MAX);
            let tpc = args
                .three_prime_clip
                .provided_or_i64("overriding 3' clip with user-provided value", i64::MAX);
            let st = args
                .score_threshold
                .provided_or_f32("overriding score threshold with user-provided value", 0_f32);
            let maf = args.min_aligned_fraction.provided_or_f32(
                "overriding min aligned fraction with user-provided value",
                0_f32,
            );
            let mal = args.min_aligned_len.provided_or_u32(
                "overriding min aligned length with user-provided value",
                1_u32,
            );

            Ok(AlignmentFilters::builder()
                .five_prime_clip(fpc)
                .three_prime_clip(tpc)
                .score_threshold(st)
                .min_aligned_fraction(maf)
                .min_aligned_len(mal)
                .which_strand(args.strand_filter)
                .model_coverage(args.model_coverage)
                .logistic_growth_rate(args.growth_rate)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .score_prob_denom(args.score_prob_denom.unwrap_or(5.0))
                .build())
        }
        Some(FilterGroup::NanocountFilters) => {
            info!("setting filters to nanocount defaults.");
            // override individual parameters if the user passed them in explicitly
            let fpc = args
                .five_prime_clip
                .provided_or_u32("overriding 5' clip with user-provided value", u32::MAX);
            let tpc = args
                .three_prime_clip
                .provided_or_i64("overriding 3' clip with user-provided value", 50_i64);
            let st = args.score_threshold.provided_or_f32(
                "overriding score threshold with user-provided value",
                0.95_f32,
            );
            let maf = args.min_aligned_fraction.provided_or_f32(
                "overriding min aligned fraction with user-provided value",
                0.5_f32,
            );
            let mal = args.min_aligned_len.provided_or_u32(
                "overriding min aligned length with user-provided value",
                50_u32,
            );

            Ok(AlignmentFilters::builder()
                .five_prime_clip(fpc)
                .three_prime_clip(tpc)
                .score_threshold(st)
                .min_aligned_fraction(maf)
                .min_aligned_len(mal)
                .which_strand(bio_types::strand::Strand::Forward)
                .model_coverage(args.model_coverage)
                .logistic_growth_rate(args.growth_rate)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .score_prob_denom(args.score_prob_denom.unwrap_or(5.0))
                .build())
        }
        None => {
            info!("setting user-provided filter parameters.");
            Ok(AlignmentFilters::builder()
                .five_prime_clip(args.five_prime_clip.try_as_u32()?)
                .three_prime_clip(args.three_prime_clip.try_as_i64()?)
                .score_threshold(args.score_threshold.try_as_f32()?)
                .min_aligned_fraction(args.min_aligned_fraction.try_as_f32()?)
                .min_aligned_len(args.min_aligned_len.try_as_u32()?)
                .which_strand(args.strand_filter)
                .model_coverage(args.model_coverage)
                .logistic_growth_rate(args.growth_rate)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .score_prob_denom(args.score_prob_denom.unwrap_or(5.0))
                .build())
        }
    }
}

/// Genome-BAM mode: project a spliced, genome-aligned (name-collated) BAM onto
/// the transcripts in `--annotation` and quantify.
fn run_genome_alignments(args: &Args, filter_opts: AlignmentFilters) -> anyhow::Result<()> {
    let bam_path = args
        .genome_alignments
        .clone()
        .expect("genome_alignments present");
    let annotation = args
        .annotation
        .clone()
        .expect("--annotation is required with --genome-alignments");

    info!("oarfish is operating in genome-alignment (projection) mode");

    // open the genome BAM with a multithreaded bgzf decoder
    let afile = File::open(&bam_path)?;
    let decomp_threads = 1.max(args.threads.saturating_sub(1));
    let worker_count = NonZeroUsize::new(decomp_threads).expect("decompression threads >= 1");
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, afile);
    let mut reader = bam::io::Reader::from(decoder);

    // the genome header gives us the chromosome -> RefId mapping; also enforce
    // that the BAM is name-collated (required for per-read-name projection).
    let genome_header = alignment_parser::read_and_verify_genome_header(&mut reader, &bam_path)?;
    let refnames: Vec<String> = genome_header
        .reference_sequences()
        .iter()
        .map(|(rseq, _)| rseq.to_string())
        .collect();

    // build the genome->transcriptome index and the transcriptome header/info.
    // genome-alignments (BAM) mode has no aligner index, so the rescue reference
    // comes from a FASTA (`--genome-fasta`) only.
    let rescue_db = projection::load_rescue_fasta(args)?;
    let use_fasta = rescue_db.is_some();
    let g2t = projection::load_g2t(&annotation, &refnames, rescue_db)?;
    let (txp_header, mut txps, txps_name) =
        projection::build_transcriptome_header_and_info(&g2t, args)?;
    let proj_config = projection::projection_config(args, use_fasta);

    let seqcol_digest = digest_utils::digest_from_header(&txp_header)?;
    let digest: NamedDigestVec = vec![("transcriptome_digest".to_string(), seqcol_digest)].into();

    bulk::quantify_genome_alignments_from_bam(
        &genome_header,
        &txp_header,
        &g2t,
        &proj_config,
        filter_opts,
        &mut reader,
        &mut txps,
        &txps_name,
        args,
        digest,
    )
}

/// Genome-read mode: spliced-align raw reads to the genome with the rammap
/// aligner, then project onto the transcripts in `--annotation` and quantify.
fn run_genome_reads(args: &mut Args, filter_opts: AlignmentFilters) -> anyhow::Result<()> {
    let annotation = args
        .annotation
        .clone()
        .expect("--annotation is required with --genome");
    let read_paths = args
        .reads
        .clone()
        .expect("--reads is required in genome read mode");

    // parse the annotation once; reused for both junction extraction and the
    // genome->transcriptome index.
    info!("loading annotation from {}", annotation.display());
    let transcripts = bramble_rs::annotation::load_transcripts(&annotation)?;
    info!("loaded {} transcripts from annotation", transcripts.len());

    // Determine the splice-junction BED to bias spliced alignment toward
    // annotated junctions (rammap `load_junctions_bed`).
    // Prefer a user-supplied BED; otherwise derive a BED12 of transcript models
    // from the annotation (default on). Computed before the aligner is built so
    // it can be loaded into the index by the aligner builder.
    let junc_bed: Option<std::path::PathBuf> = if let Some(b) = args.junctions.clone() {
        Some(b)
    } else if !args.ignore_annotation_junctions {
        let mut bed = args.output.clone().expect("output prefix required");
        let fname = format!(
            "{}.annot_junctions.bed",
            bed.file_name()
                .map(|f| f.to_string_lossy().into_owned())
                .unwrap_or_default()
        );
        bed.set_file_name(fname);
        let n = projection::write_annotation_junction_bed(&transcripts, &bed)?;
        info!(
            "derived {} spliced transcript models from the annotation for the splice-junction BED",
            n
        );
        Some(bed)
    } else {
        info!("not using annotated splice junctions (--ignore-annotation-junctions)");
        None
    };

    // build the spliced genome aligner (loading the junction BED into the index)
    // and obtain its reference (chromosome) names in target-id order (== the g2t
    // RefId order).
    let (aligner, refnames) =
        crate::util::aligner::get_genome_aligner_from_args(args, junc_bed.as_deref())?;

    // build the genome->transcriptome index and the transcriptome header/info.
    // Rescue is on by default: source its reference from a FASTA (`--genome-fasta`,
    // or a FASTA `--genome`), else from the reference sequences resident in the
    // loaded genome index (`--genome` is a prebuilt index). `--no-rescue` / no
    // available source => off.
    let rescue_db = match projection::load_rescue_fasta(args)? {
        Some(db) => Some(db),
        None if !args.no_rescue => crate::util::aligner::rescue_db_from_genome_index(&aligner),
        None => None,
    };
    let use_fasta = rescue_db.is_some();
    let g2t = projection::build_g2t_from_transcripts(&transcripts, &refnames, rescue_db)?;
    let (txp_header, mut txps, txps_name) =
        projection::build_transcriptome_header_and_info(&g2t, args)?;
    let proj_config = projection::projection_config(args, use_fasta);

    let seqcol_digest = digest_utils::digest_from_header(&txp_header)?;
    let digest: NamedDigestVec = vec![("transcriptome_digest".to_string(), seqcol_digest)].into();

    bulk::quantify_genome_raw_reads(
        &txp_header,
        aligner,
        &g2t,
        &proj_config,
        filter_opts,
        &read_paths,
        &mut txps,
        &txps_name,
        args,
        digest,
    )
}

fn main() -> anyhow::Result<()> {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();
    let (filtered_layer, reload_handle) = tracing_subscriber::reload::Layer::new(env_filter);

    // set up the logging.  Here we will take the
    // logging level from the environment variable if
    // it is set.  Otherwise, we'll set the default
    tracing_subscriber::registry()
        // log level to INFO.
        .with(fmt::layer().with_writer(io::stderr))
        .with(filtered_layer)
        .init();

    let mut args = Args::parse();

    if args.model_coverage && !args.single_cell {
        args.coverage_model = crate::prog_opts::CoverageModel::Logistic;
    } else if matches!(
        args.coverage_model,
        crate::prog_opts::CoverageModel::Logistic
            | crate::prog_opts::CoverageModel::Hybrid
            | crate::prog_opts::CoverageModel::Adaptive
            | crate::prog_opts::CoverageModel::Degradation
            | crate::prog_opts::CoverageModel::Auto
    ) {
        // Only the historical model needs per-transcript coverage bins.  The
        // endpoint model learns directly from retained alignment coordinates.
        args.model_coverage = true;
    }
    if matches!(
        args.coverage_model,
        crate::prog_opts::CoverageModel::Hybrid
            | crate::prog_opts::CoverageModel::Adaptive
            | crate::prog_opts::CoverageModel::Degradation
            | crate::prog_opts::CoverageModel::Auto
    ) && args.logistic_weight == 0.0
        && args.endpoint_weight == 0.0
    {
        anyhow::bail!(
            "hybrid/adaptive/degradation coverage requires a nonzero logistic or endpoint weight"
        );
    }
    if args.coverage_model == crate::prog_opts::CoverageModel::Degradation
        && args.seq_tech != Some(crate::prog_opts::SequencingTech::OntDRNA)
    {
        anyhow::bail!("--coverage-model degradation currently requires --seq-tech ont-drna");
    }
    if args.coverage_model == crate::prog_opts::CoverageModel::Auto && args.seq_tech.is_none() {
        anyhow::bail!("--coverage-model auto requires --seq-tech so it can select a technology kernel");
    }
    if args.degradation_kernel != crate::prog_opts::DegradationKernel::Constant
        && (!matches!(
            args.coverage_model,
            crate::prog_opts::CoverageModel::Degradation
                | crate::prog_opts::CoverageModel::Auto
        ) || args.seq_tech != Some(crate::prog_opts::SequencingTech::OntDRNA))
    {
        anyhow::bail!(
            "non-constant --degradation-kernel values require ONT direct-RNA auto/degradation coverage"
        );
    }
    if args.single_cell && args.coverage_model != crate::prog_opts::CoverageModel::None {
        anyhow::bail!(
            "--coverage-model is currently supported only for bulk quantification; use the legacy --model-coverage single-cell path"
        );
    }

    if args.single_cell && args.em_accel != crate::prog_opts::EmAccel::None {
        anyhow::bail!("--em-accel is currently supported only for bulk quantification");
    }

    // change the logging filter if the user specified quiet or
    // verbose.
    if args.quiet {
        reload_handle.modify(|filter| *filter = EnvFilter::new("WARN"))?;
    }
    if args.verbose {
        reload_handle.modify(|filter| *filter = EnvFilter::new("TRACE"))?;
    }

    // Apply an explicit DP scratch-cache cap if the user set one; otherwise the
    // rammap mapper uses its default (128 MB) or the RAMMAP_DP_CACHE_CAP_MB env var.
    if let Some(mb) = args.dp_cache_cap_mb {
        rammap::set_dp_cache_cap_mb(mb);
    }

    // if we are just indexing, don't bother with anything else
    if args.only_index {
        let (header, _reader, _aligner, digest) = get_aligner_from_fastas(&mut args)?;
        info!(
            "indexing completed; index over {} references written to {}",
            header.reference_sequences().len(),
            &args.index_out.expect("nonempty").display()
        );
        info!(
            "reference digest = {}",
            serde_json::to_string_pretty(&digest.to_json())?
        );
        return Ok(());
    }

    // if we are quantifying, handle this below
    let filter_opts = get_filter_opts(&args)?;

    // Genome-space modes: rather than quantifying against the input BAM/aligner
    // references directly, we project genomic alignments onto the transcripts in
    // `--annotation` (via bramble) and quantify those. These paths build the
    // transcriptome header / `TranscriptInfo` from the annotation, not from the
    // input header, so they are handled separately and return early.
    if args.genome_alignments.is_some() {
        run_genome_alignments(&args, filter_opts)?;
        info!("oarfish completed successfully.");
        return Ok(());
    } else if args.genome.is_some() {
        run_genome_reads(&mut args, filter_opts)?;
        info!("oarfish completed successfully.");
        return Ok(());
    }

    let (header, reader, aligner, digest) = if args.alignments.is_none() {
        get_aligner_from_args(&mut args)?
    } else {
        let alignments = args.alignments.clone().unwrap();
        let afile = File::open(&alignments)?;

        let decomp_threads = if args.single_cell {
            // we will overlap quantification with parsing, so don't try to use too many
            // parser threads, and adjust the worker threads accordingly.

            // is there a better heuristic than this?
            // <= 6 threads, use only 1 for decompression
            // 6-8 threads, use 2 for decompression
            // > 8 threads, use 3 for decompression
            match args.threads {
                1..=6 => 1,
                7 | 8 => 2,
                _ => 3,
            }
        } else {
            // try to use all but 1 thread, and assume we have at least 2.
            1.max(args.threads.saturating_sub(1))
        };

        let worker_count = NonZeroUsize::new(decomp_threads).expect("decompression threads >= 1");
        if args.single_cell {
            args.threads = 1.max(args.threads.saturating_sub(decomp_threads));
        }

        let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, afile);
        let mut reader = bam::io::Reader::from(decoder);
        // parse the header, and ensure that the reads were mapped with a known long-read
        // aligner (as far as we can tell).
        let header = alignment_parser::read_and_verify_header(&mut reader, &alignments)?;
        let seqcol_digest = digest_utils::digest_from_header(&header)?;
        (
            header,
            Some(reader),
            None,
            vec![("bam_digest".to_string(), seqcol_digest)].into(),
        )
    };

    // where we'll write down the per-transcript information we need
    // to track.
    let (mut txps, txps_name) = get_txp_info_from_header(&header, &args)?;

    if args.single_cell {
        // TODO: do this better (quiet the EM during single-cell quant)
        reload_handle.modify(|filter| {
            *filter = if args.quiet {
                EnvFilter::new("WARN")
                    .add_directive("oarfish=warn".parse().unwrap())
                    .add_directive("oarfish::single_cell=warn".parse().unwrap())
            } else if args.verbose {
                EnvFilter::new("TRACE")
                    .add_directive("oarfish=info".parse().unwrap())
                    .add_directive("oarfish::single_cell=trace".parse().unwrap())
            } else {
                // be quiet about normal things in single-cell mode
                // e.g. EM iterations, and only print out info for
                // oarfish::single_cell events.
                EnvFilter::new("INFO")
                    .add_directive("oarfish=warn".parse().unwrap())
                    .add_directive("oarfish::single_cell=info".parse().unwrap())
            }
        })?;

        single_cell::quantify_single_cell_from_collated_bam(
            &header,
            &filter_opts,
            &mut reader.unwrap(),
            &mut txps,
            &args,
            digest,
        )?;
    } else if args.alignments.is_some() {
        bulk::quantify_bulk_alignments_from_bam(
            &header,
            filter_opts,
            &mut reader.unwrap(),
            &mut txps,
            &txps_name,
            &args,
            digest,
        )?;
    } else {
        bulk::quantify_bulk_alignments_raw_reads(
            &header,
            aligner.expect("need valid alinger to align reads"),
            filter_opts,
            &args.reads.clone().expect("expected read file(s)"),
            &mut txps,
            &txps_name,
            &args,
            digest,
        )?;
    }

    info!("oarfish completed successfully.");
    Ok(())
}
