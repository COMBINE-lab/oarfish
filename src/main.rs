use clap::Parser;
use noodles_sam::Header;
use prog_opts::FragmentEndModel;
use std::num::NonZeroUsize;
use util::oarfish_types::{FragmentEndFalloffDist, NamedDigestVec};

// Or now
// use minimap2::ffi as mm_ffi;
//use minimap2_temp as minimap2;
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
mod prog_opts;
mod single_cell;
mod util;

use crate::prog_opts::{Args, FilterGroup};
use crate::util::aligner::{get_aligner_from_args, get_aligner_from_fastas};
use crate::util::digest_utils;
use crate::util::kde_utils;
use crate::util::oarfish_types::{AlignmentFilters, TranscriptInfo};

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
    // if the user disabled the fragment end distance modeling, then
    // we will pass None to the filter group for this parameter, otherwise
    // we will fill it out with the provided standard deviation and threshold.
    let frag_dist = match &args.fragment_end_model {
        FragmentEndModel::None => {
            info!("disabled fragment end-distance modeling.");
            None
        }
        x => {
            let mname = match x {
                FragmentEndModel::FivePrimeStart(d) => format!("5' end start {}", d),
                FragmentEndModel::ThreePrimeStart(d) => format!("3' end start {}", d),
                FragmentEndModel::EitherStart(d) => format!("either end start {}", d),
                FragmentEndModel::None => {
                    anyhow::bail!("Must have a fragment end model at this point")
                }
            };
            info!(
                "applying {} fragment end-distance modeling with standard deviation = {} and threshold = {}.",
                mname, args.end_dist_std_dev, args.end_dist_thresh
            );

            let frag_min_prob = args.frag_min_prob.ln();

            Some(FragmentEndFalloffDist::new(
                0.,
                args.end_dist_std_dev,
                args.end_dist_thresh,
                args.fragment_end_model,
                frag_min_prob,
            ))
        }
    };

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
                .falloff_dist(frag_dist)
                .alignment_score_denom(args.alignment_score_denom)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
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
                .falloff_dist(frag_dist)
                .alignment_score_denom(args.alignment_score_denom)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
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
                .falloff_dist(frag_dist)
                .alignment_score_denom(args.alignment_score_denom)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .build())
        }
    }
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

    // change the logging filter if the user specified quiet or
    // verbose.
    if args.quiet {
        reload_handle.modify(|filter| *filter = EnvFilter::new("WARN"))?;
    }
    if args.verbose {
        reload_handle.modify(|filter| *filter = EnvFilter::new("TRACE"))?;
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
        // parse the header, and ensure that the reads were mapped with minimap2 (as far as we
        // can tell).
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
