use clap::Parser;

use anyhow::Context;
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};

use std::{
    fs::File,
    io::{self, BufReader},
    path::PathBuf,
};

use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

use noodles_bam as bam;

mod alignment_parser;
mod bootstrap;
mod em;
mod util;

use crate::util::binomial_probability::binomial_continuous_prob;
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_output};

/// These represent different "meta-options", specific settings
/// for all of the different filters that should be applied in
/// different cases.
#[derive(Clone, Debug, clap::ValueEnum, Serialize)]
enum FilterGroup {
    NoFilters,
    NanocountFilters,
}

/// accurate transcript quantification from long-read RNA-seq data
#[derive(Parser, Debug, Serialize)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// be quiet (i.e. don't output log messages that aren't at least warnings)
    #[arg(long, conflicts_with = "verbose")]
    quiet: bool,

    /// be verbose (i.e. output all non-developer logging messages)
    #[arg(long)]
    verbose: bool,

    /// path to the file containing the input alignments
    #[arg(short, long, required = true)]
    alignments: PathBuf,
    /// location where output quantification file should be written
    #[arg(short, long, required = true)]
    output: PathBuf,

    #[arg(long, help_heading = "filters", value_enum)]
    filter_group: Option<FilterGroup>,

    /// maximum allowable distance of the right-most end of an alignment from the 3' transcript end
    #[arg(short, long, conflicts_with = "filter-group", help_heading="filters", default_value_t = u32::MAX as i64)]
    three_prime_clip: i64,
    /// maximum allowable distance of the left-most end of an alignment from the 5' transcript end
    #[arg(short, long, conflicts_with = "filter-group", help_heading="filters", default_value_t = u32::MAX)]
    five_prime_clip: u32,
    /// fraction of the best possible alignment score that a secondary alignment must have for
    /// consideration
    #[arg(
        short,
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = 0.95
    )]
    score_threshold: f32,
    /// fraction of a query that must be mapped within an alignemnt to consider the alignemnt
    /// valid
    #[arg(
        short,
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = 0.5
    )]
    min_aligned_fraction: f32,
    /// minimum number of nucleotides in the aligned portion of a read
    #[arg(
        short = 'l',
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = 50
    )]
    min_aligned_len: u32,
    /// allow both forward-strand and reverse-complement alignments
    #[arg(
        short = 'n',
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        value_parser
    )]
    allow_negative_strand: bool,
    /// apply the coverage model
    #[arg(long, help_heading = "coverage model", value_parser)]
    model_coverage: bool,
    /// maximum number of iterations for which to run the EM algorithm
    #[arg(long, help_heading = "EM", default_value_t = 1000)]
    max_em_iter: u32,
    /// maximum number of iterations for which to run the EM algorithm
    #[arg(long, help_heading = "EM", default_value_t = 1e-3)]
    convergence_thresh: f64,
    /// maximum number of cores that the oarfish can use to obtain binomial probability
    #[arg(short = 'j', long, default_value_t = 1)]
    threads: usize,
    /// location of short read quantification (if provided)
    #[arg(short = 'q', long, help_heading = "EM")]
    short_quant: Option<String>,
    /// number of bootstrap replicates to produce to assess quantification uncertainty
    #[arg(long, default_value_t = 0)]
    num_bootstraps: u32,
    /// number of bins to use in coverage model
    #[arg(short, long, help_heading = "coverage model", default_value_t = 10)]
    bins: u32,
}

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
                .allow_rc(true)
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
                .allow_rc(false)
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
                .allow_rc(args.allow_negative_strand)
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
        "binomial"
    } else {
        "no_coverage"
    };

    json!({
        "prob_model" : prob,
        "num_bins" : args.bins,
        "filter_options" : &emi.eq_map.filter_opts,
        "discard_table" : &emi.eq_map.discard_table,
        "alignments": &args.alignments,
        "output": &args.output,
        "verbose": &args.verbose,
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

    let mut reader = File::open(&args.alignments)
        .map(BufReader::new)
        .map(bam::io::Reader::new)?;

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
            txps.push(TranscriptInfo::with_len_and_bins(rmap.length(), args.bins));
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

    // now parse the actual alignments for the reads and store the results
    // in our in-memory stor
    let mut store = InMemoryAlignmentStore::new(filter_opts, &header);
    alignment_parser::parse_alignments(&mut store, &header, &mut reader, &mut txps)?;

    // print discard table information in which the user might be interested.
    info!("\ndiscard_table: \n{}\n", store.discard_table.to_table());

    if store.filter_opts.model_coverage {
        //obtaining the Cumulative Distribution Function (CDF) for each transcript
        binomial_continuous_prob(&mut txps, &args.bins, args.threads);
        //Normalize the probabilities for the records of each read
        normalize_read_probs(&mut store, &txps, &args.bins);
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
        txp_info: &mut txps,
        max_iter: args.max_em_iter,
        convergence_thresh: args.convergence_thresh,
        init_abundances,
    };

    let counts = if args.threads > 4 {
        em::em_par(&emi, args.threads)
    } else {
        em::em(&emi, args.threads)
    };

    // prepare the JSON object we'll write
    // to meta_info.json
    let json_info = get_json_info(&args, &emi, &seqcol_digest);

    // write the output
    write_output(&args.output, json_info, &header, &counts)?;

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

    Ok(())
}
