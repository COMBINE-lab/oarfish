use clap::Parser;
use core::str;
use std::num::NonZeroUsize;

use anyhow::Context;
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};

use std::{fs::File, io};

use num_format::{Locale, ToFormattedString};
use serde_json::json;
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

use noodles_bam as bam;
use noodles_bgzf as bgzf;

mod alignment_parser;
mod bootstrap;
mod em;
mod prog_opts;
mod single_cell;
mod util;

use crate::prog_opts::{Args, FilterGroup};
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_output};
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

/// Produce a [serde_json::Value] that encodes the relevant arguments and
/// parameters of the run that we wish to record to file. Ultimately, this
/// will be written to the corresponding `meta_info.json` file for this run.
fn get_json_info(args: &Args, emi: &EMInfo, seqcol_digest: &str) -> serde_json::Value {
    let prob = if args.model_coverage {
        "scaled_binomial"
    } else {
        "no_coverage"
    };

    json!({
        "prob_model" : prob,
        "bin_width" : args.bin_width,
        "filter_options" : &emi.eq_map.filter_opts,
        "discard_table" : &emi.eq_map.discard_table,
        "alignments": &args.alignments,
        "output": &args.output,
        "verbose": &args.verbose,
        "single_cell": &args.single_cell,
        "quiet": &args.quiet,
        "em_max_iter": &args.max_em_iter,
        "em_convergence_thresh": &args.convergence_thresh,
        "threads": &args.threads,
        "filter_group": &args.filter_group,
        "short_quant": &args.short_quant,
        "num_bootstraps": &args.num_bootstraps,
        "seqcol_digest": seqcol_digest
    })
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
        reload_handle.modify(|filter| *filter = EnvFilter::new("WARN"))?;
        single_cell::quantify_single_cell_from_collated_bam(
            &header,
            &filter_opts,
            &mut reader,
            &mut txps,
            &args,
            seqcol_digest,
        )?;
    } else {
        // now parse the actual alignments for the reads and store the results
        // in our in-memory stor
        let mut store = InMemoryAlignmentStore::new(filter_opts, &header);
        alignment_parser::parse_alignments(&mut store, &header, &mut reader, &mut txps)?;

        // print discard table information in which the user might be interested.
        info!("\ndiscard_table: \n{}\n", store.discard_table.to_table());

        // no longer need the reader
        drop(reader);

        // if we are using the KDE, create that here.
        let kde_opt: Option<kders::kde::KDEModel> = if args.use_kde {
            Some(kde_utils::get_kde_model(&txps, &store)?)
        } else {
            None
        };

        if store.filter_opts.model_coverage {
            //obtaining the Cumulative Distribution Function (CDF) for each transcript
            binomial_continuous_prob(&mut txps, &args.bin_width, args.threads);
            //Normalize the probabilities for the records of each read
            normalize_read_probs(&mut store, &txps, &args.bin_width);
        }

        info!(
            "Total number of alignment records : {}",
            store.total_len().to_formatted_string(&Locale::en)
        );
        info!(
            "number of aligned reads : {}",
            store.num_aligned_reads().to_formatted_string(&Locale::en)
        );
        info!(
            "number of unique alignments : {}",
            store.unique_alignments().to_formatted_string(&Locale::en)
        );

        // if we are seeding the quantification estimates with short read
        // abundances, then read those in here.
        let init_abundances = args.short_quant.as_ref().map(|sr_path| {
            read_short_quant_vec(sr_path, &txps_name).unwrap_or_else(|e| panic!("{}", e))
        });

        // wrap up all of the relevant information we need for estimation
        // in an EMInfo struct and then call the EM algorithm.
        let emi = EMInfo {
            eq_map: &store,
            txp_info: &txps,
            max_iter: args.max_em_iter,
            convergence_thresh: args.convergence_thresh,
            init_abundances,
            kde_model: kde_opt,
        };

        if args.use_kde {
            /*
            // run EM for model train iterations
            let orig_iter = emi.max_iter;
            emi.max_iter = 10;
            let counts = em::em(&emi, args.threads);
            // relearn the kde
            let new_model =
            kde_utils::refresh_kde_model(&txps, &store, &emi.kde_model.unwrap(), &counts);
            info!("refreshed KDE model");
            emi.kde_model = Some(new_model?);
            emi.max_iter = orig_iter;
            */
        }

        let counts = if args.threads > 4 {
            em::em_par(&emi, args.threads)
        } else {
            em::em(&emi, args.threads)
        };

        let aux_txp_counts = crate::util::aux_counts::get_aux_counts(&store, &txps)?;

        // prepare the JSON object we'll write
        // to meta_info.json
        let json_info = get_json_info(&args, &emi, &seqcol_digest);

        // write the output
        write_output(&args.output, json_info, &header, &counts, &aux_txp_counts)?;

        // if the user requested bootstrap replicates,
        // compute and write those out now.
        if args.num_bootstraps > 0 {
            let breps = em::bootstrap(&emi, args.num_bootstraps, args.threads);

            let mut new_arrays = vec![];
            let mut bs_fields = vec![];
            for (i, b) in breps.into_iter().enumerate() {
                let bs_array = Float64Array::from_vec(b);
                bs_fields.push(Field::new(
                    format!("bootstrap.{}", i),
                    bs_array.data_type().clone(),
                    false,
                ));
                new_arrays.push(bs_array.boxed());
            }
            let chunk = Chunk::new(new_arrays);
            write_infrep_file(&args.output, bs_fields, chunk)?;
        }
    }

    Ok(())
}
