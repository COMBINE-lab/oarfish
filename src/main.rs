use clap::Parser;
use std::num::NonZeroUsize;

use anyhow::Context;

use std::{fs::File, io};

use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

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
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{AlignmentFilters, TranscriptInfo};
use crate::util::{binomial_probability::binomial_continuous_prob, kde_utils};

fn get_filter_opts(args: &Args) -> AlignmentFilters {
    // set all of the filter options that the user
    // wants to apply.
    match args.filter_group {
        Some(FilterGroup::NoFilters) => {
            info!("disabling alignment filters.");
            AlignmentFilters::builder()
                .five_prime_clip(u32::MAX)
                .three_prime_clip(i64::MAX)
                .score_threshold(0_f32)
                .min_aligned_fraction(0_f32)
                .min_aligned_len(1_u32)
                .which_strand(args.strand_filter)
                .model_coverage(args.model_coverage)
                .build()
        }
        Some(FilterGroup::NanocountFilters) => {
            info!("setting filters to nanocount defaults.");
            AlignmentFilters::builder()
                .five_prime_clip(u32::MAX)
                .three_prime_clip(50_i64)
                .score_threshold(0.95_f32)
                .min_aligned_fraction(0.5_f32)
                .min_aligned_len(50_u32)
                .which_strand(bio_types::strand::Strand::Forward)
                .model_coverage(args.model_coverage)
                .build()
        }
        None => {
            info!("setting user-provided filter parameters.");
            AlignmentFilters::builder()
                .five_prime_clip(args.five_prime_clip)
                .three_prime_clip(args.three_prime_clip)
                .score_threshold(args.score_threshold)
                .min_aligned_fraction(args.min_aligned_fraction)
                .min_aligned_len(args.min_aligned_len)
                .which_strand(args.strand_filter)
                .model_coverage(args.model_coverage)
                .build()
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

    let args = Args::parse();

    // change the logging filter if the user specified quiet or
    // verbose.
    if args.quiet {
        reload_handle.modify(|filter| *filter = EnvFilter::new("WARN"))?;
    }
    if args.verbose {
        reload_handle.modify(|filter| *filter = EnvFilter::new("TRACE"))?;
    }

    let filter_opts = get_filter_opts(&args);

    let afile = File::open(&args.alignments)?;

    let worker_count = NonZeroUsize::new(1.max(args.threads.saturating_sub(1)))
        .expect("decompression threads >= 1");
    let decoder = bgzf::MultithreadedReader::with_worker_count(worker_count, afile);
    let mut reader = bam::io::Reader::from(decoder);
    // parse the header, and ensure that the reads were mapped with minimap2 (as far as we
    // can tell).
    let header = alignment_parser::read_and_verify_header(&mut reader, &args.alignments)?;
    let seqcol_digest = {
        info!("calculating seqcol digest");
        let sc = seqcol_rs::SeqCol::from_sam_header(
            header
                .reference_sequences()
                .iter()
                .map(|(k, v)| (k.as_slice(), v.length().into())),
        );
        let d = sc.digest(seqcol_rs::DigestConfig::default()).context(
            "failed to compute the seqcol digest for the information from the alignment header",
        )?;
        info!("done calculating seqcol digest");
        d
    };
    let num_ref_seqs = header.reference_sequences().len();

    // where we'll write down the per-transcript information we need
    // to track.
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
        txps.len()
    );

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
                EnvFilter::new("INFO")
                    .add_directive("oarfish=warn".parse().unwrap())
                    .add_directive("oarfish::single_cell=info".parse().unwrap())
            }
        })?;

        single_cell::quantify_single_cell_from_collated_bam(
            &header,
            &filter_opts,
            &mut reader,
            &mut txps,
            &args,
            seqcol_digest,
        )?;
    } else {
        bulk::quantify_bulk_alignments_from_bam(
            &header,
            filter_opts,
            &mut reader,
            &mut txps,
            &txps_name,
            &args,
            seqcol_digest,
        )?;
    }

    Ok(())
}
