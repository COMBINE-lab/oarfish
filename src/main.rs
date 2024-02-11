use clap::Parser;

use std::{
    fs::File,
    io::{self, BufReader},
    path::PathBuf,
};

use num_format::{Locale, ToFormattedString};
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

use noodles_bam as bam;

mod util;
mod em;
use crate::util::binomial_probability::binomial_continuous_prob;
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::write_output;


/// These represent different "meta-options", specific settings 
/// for all of the different filters that should be applied in 
/// different cases.
#[derive(Clone, Debug, clap::ValueEnum)]
enum FilterGroup {
    NoFilters,
    NanocountFilters,
}

/// accurate transcript quantification from long-read RNA-seq data
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// path to the file containing the input alignments
    #[arg(short, long, required = true)]
    alignments: PathBuf,
    /// location where output quantification file should be written
    #[arg(short, long, required = true)]
    output: String,

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
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    /// location of short read quantification (if provided)
    #[arg(short = 'q', long, help_heading = "EM")]
    short_quant: Option<String>,
    /// number of bins to use in coverage model
    #[arg(short, long, help_heading = "coverage model", default_value_t = 10)]
    bins: u32,
}

fn main() -> io::Result<()> {
    // set up the logging.  Here we will take the
    // logging level from the environment variable if
    // it is set.  Otherwise, we'll set the default
    tracing_subscriber::registry()
        // log level to INFO.
        .with(fmt::layer().with_writer(io::stderr))
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let args = Args::parse();

    // set all of the filter options that the user
    // wants to apply.
    let filter_opts = match args.filter_group {
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
    };

    let mut reader = File::open(&args.alignments)
        .map(BufReader::new)
        .map(bam::io::Reader::new)?;

    // read the bam file header, print out some basic info
    let header = reader.read_header()?;
    info!(
        "read header from BAM file {}, contains {} reference sequences.",
        args.alignments.display(),
        header
            .reference_sequences()
            .len()
            .to_formatted_string(&Locale::en)
    );

    let mut saw_minimap2 = false;
    let mut progs = vec![];
    // explicitly check that alignment was done with a supported
    // aligner (right now, just minimap2).
    for (prog, _pmap) in header.programs().iter() {
        if prog == "minimap2" {
            saw_minimap2 = true;
            break;
        } else {
            progs.push(prog);
        }
    }
    assert!(
        saw_minimap2,
        "Currently, only minimap2 is supported as an aligner. The bam file listed {:?}.",
        progs
    );
    info!("saw minimap2 as a program in the header; proceeding.");

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
    info!("parsed reference information for {} transcripts.", txps.len());

    // we'll need these to keep track of which alignments belong
    // to which reads.
    let mut prev_read = String::new();
    let mut num_unmapped = 0_u64;
    let mut records_for_read = vec![];
    let mut store = InMemoryAlignmentStore::new(filter_opts, &header);

    // Parse the input alignemnt file, gathering the alignments aggregated
    // by their source read. **Note**: this requires that we have a
    // name-sorted input bam file (currently, aligned against the transcriptome).
    //
    // *NOTE*: this had to be changed from `records` to `record_bufs` or
    // critical information was missing from the records. This happened when
    // moving to the new version of noodles. Track `https://github.com/zaeleus/noodles/issues/230`
    // to see if it's clear why this is the case
    for (_i, result) in reader.record_bufs(&header).enumerate() {
        let record = result?;
        // unmapped reads don't contribute to quantification
        // but we track them.
        if record.flags().is_unmapped() {
            num_unmapped += 1;
            continue;
        }
        let record_copy = record.clone();
        if let Some(rname) = record.name() {
            let rstring: String = String::from_utf8_lossy(rname.as_ref()).into_owned();
            // if this is an alignment for the same read, then
            // push it onto our temporary vector.
            if prev_read == rstring {
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            } else {
                // otherwise, record the alignment range for the
                // previous read record.
                if !prev_read.is_empty() {
                    store.add_group(&mut txps, &mut records_for_read);
                    records_for_read.clear();
                }
                // the new "prev_read" name is the current read name
                // so it becomes the first on the new alignment range
                // vector.
                prev_read = rstring;
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            }
        }
    }
    // if we end with a non-empty alignment range vector, then
    // add that group.
    if !records_for_read.is_empty() {
        store.add_group(&mut txps, &mut records_for_read);
        records_for_read.clear();
    }

    info!(
        "the alignment file contained {} unmapped read records.",
        num_unmapped.to_formatted_string(&Locale::en)
    );
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

    let init_abundances = args.short_quant.map(|sr_path| {
        read_short_quant_vec(&sr_path, &txps_name).unwrap_or_else(|e| panic!("{}", e))
    });

    // wrap up all of the relevant information we need for estimation
    // in an EMInfo struct and then call the EM algorithm.
    let mut emi = EMInfo {
        eq_map: &store,
        txp_info: &mut txps,
        max_iter: args.max_em_iter,
        convergence_thresh: args.convergence_thresh,
        init_abundances,
    };

    let counts = em::em(&mut emi);

    // write the output
    write_output(
        &args.output,
        &args.model_coverage,
        &args.bins,
        &emi,
        &header,
        &counts,
    )?;

    Ok(())
}
