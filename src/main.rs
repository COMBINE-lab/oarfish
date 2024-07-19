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

mod alignment_parser;
mod em;
mod util;
use crate::util::binomial_probability::binomial_continuous_prob;
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo, entropy_function, Fragment, construct_equivalence_classes,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_output, write_cdf, write_info, write_EquivalenceClass};
use crate::util::kde_function;
use itertools::izip;
use num::integer::gcd;

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
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    /// location of short read quantification (if provided)
    #[arg(short = 'q', long, help_heading = "EM")]
    short_quant: Option<String>,
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
    let num_ref_seqs = header.reference_sequences().len();

    // where we'll write down the per-transcript information we need
    // to track.
    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(num_ref_seqs);
    let mut txps_name: Vec<String> = Vec::with_capacity(num_ref_seqs);
    let mut txp_len_gcd: Option<usize> = None;
    let mut txp_len_min: Option<usize> = None;
    let mut final_len_gcd: usize = 1;

    if args.model_coverage {
        for (_, rmap) in header.reference_sequences().iter() {
            let length = rmap.length().get() as usize;

            // Compute the running GCD of transcript lengths
            txp_len_gcd = Some(match txp_len_gcd {
                Some(gcd_value) => gcd(gcd_value, length),
                None => length,
            });

            // Determine the minimum length of transcripts
            txp_len_min = Some(match txp_len_min {
                Some(current_min) => current_min.min(length),
                None => length,
            });
        }

        if let Some(gcd_val) = txp_len_gcd {
            if gcd_val <= 2 {
                if let Some(min_val) = txp_len_min {
                    final_len_gcd = min_val;
                }
            } else {
                final_len_gcd = gcd_val;
            }
        }

        eprintln!("the bin length is: {:?}", final_len_gcd);
    }

    // loop over the transcripts in the header and fill in the relevant
    // information here.
    if args.model_coverage {
        for (rseq, rmap) in header.reference_sequences().iter() {
            //txps.push(TranscriptInfo::with_len_and_bins(rmap.length(), args.bins));
            txps.push(TranscriptInfo::with_len_and_bins(rmap.length(), final_len_gcd));
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

    //count the number of reads aligned to each transcripts
    for (alns, _as_probs, _coverage_prob, _as_val, _read_name, _kdeprobs) in store.iter() {
        //iterate over the alignments of a read
        for a in alns.iter() {
            let target_id: usize = a.ref_id as usize;
            let start_aln = a.start as usize;
            let end_aln = a.end as usize;
            (start_aln..=end_aln).for_each(|i| txps[target_id].txp_counts[i - 1] += 1.0);
            txps[target_id].num_read += 1.0;
        }
    }

    if store.filter_opts.model_coverage {
        //obtaining the Cumulative Distribution Function (CDF) for each transcript
        binomial_continuous_prob(&mut txps, args.threads, &txps_name, final_len_gcd);
        //Normalize the probabilities for the records of each read
        normalize_read_probs(&mut store, &txps, final_len_gcd);
    }

    info!(
        "Total number of alignment records : {}",
        store.total_len().to_formatted_string(&Locale::en)
    );
    info!(
        "number of aligned reads : {}",
        store.num_aligned_reads().to_formatted_string(&Locale::en)
    );

    // if we are seeding the quantification estimates with short read
    // abundances, then read those in here.
    let init_abundances = args.short_quant.map(|sr_path| {
        read_short_quant_vec(&sr_path, &txps_name).unwrap_or_else(|e| panic!("{}", e))
    });

    //computing the entropy of the transcripts based on the counts of each point
    for element in txps.iter_mut() {
        // Compute the sum of counts
        let sum: f64 = element.txp_counts.iter().map(|&x| x as f64).sum();
        if sum > 0.0 {
            // Normalize the counts by dividing each element by the sum
            let p: Vec<f64> = element.txp_counts.iter().map(|&count| count as f64 / sum).collect();
            (element.entropy, element.entropy_max, element.entropy_ratio, element.entropy_ratio_1)= entropy_function(&p);
            element.tin = 100.0 * (element.entropy.exp() / (element.lenf as f64))
        } else {
            element.entropy_max = (element.lenf as f64).log2();
            element.entropy_ratio_1 = 1.0;
            element.tin = 0.0;
        }
    }

    let mut fragments: Vec<Fragment> = vec![];
    for (alns, _as_probs, _coverage_prob, _as_val, read_name, _kde_probs) in store.iter() {
        let mut eq_vec = vec![];
        let mut rname: String = String::new();
        //iterate over the alignments of a read
        for (a, read_n) in izip!(alns, read_name) {
            let target_id: usize = a.ref_id as usize;
            eq_vec.push(target_id);
            rname = read_n.to_string();
        }
        fragments.push(Fragment { name: rname, mappings: eq_vec.iter().cloned().collect()});
    }

    let equivalence_classes = construct_equivalence_classes(fragments);
    //println!("equivalence class is: {:?}", equivalence_classes);
    // wrap up all of the relevant information we need for kde estimation
    
    let init2 = init_abundances.clone();
    let mut emi_kde = EMInfo {
        eq_map: &store,
        txp_info: &mut txps,
        max_iter: 10,
        convergence_thresh: args.convergence_thresh,
        init_abundances: init2,
    };
    let weight_count = em::em_par(&mut emi_kde, args.threads);

    let mut alignment_data = vec![];
    let mut alignment_weight = vec![];
    for (alns, as_probs, coverage_probs, _as_val, _read_name, _kde_probs) in store.iter() 
    {
        let mut denom = 0.0_f64;
        for (a, p, cp) in izip!(alns, as_probs, coverage_probs)
        {
            let target_id: usize = a.ref_id as usize;
            alignment_data.push(txps[target_id].lenf);
            alignment_data.push(a.alignment_span() as f64);
            alignment_weight.push(weight_count[target_id]);
            denom += weight_count[target_id] * (*p as f64) * (if args.model_coverage { *cp } else { 1.0 });
            //denom += weight_count[target_id] * (*p as f64) * (1.0);
        }
//
        if denom != 0.0 {
            // Loop over all possible assignment locations and proportionally
            // allocate the read according to our model and current parameter
            // estimates.
            for (a, p, cp) in izip!(alns, as_probs, coverage_probs)
            {
                let target_id: usize = a.ref_id as usize;
                let inc = (weight_count[target_id] * (*p as f64) * (if args.model_coverage { *cp } else { 1.0 })) / denom;
                //let inc = (weight_count[target_id] * (*p as f64) * (1.0 )) / denom;
                alignment_weight.push(inc);
            }
        } else {
            for (a, p, cp) in izip!(alns, as_probs, coverage_probs)
            {
                alignment_weight.push(0.0);
            }
        }
//
    }
//
    kde_function::kde_computation(&alignment_data, &alignment_weight, &mut store);

    // in an EMInfo struct and then call the EM algorithm.
    let mut emi = EMInfo {
        eq_map: &store,
        txp_info: &mut txps,
        max_iter: args.max_em_iter,
        convergence_thresh: args.convergence_thresh,
        init_abundances,
    };

    let counts = em::em_par(&mut emi, args.threads);

    // write the output
    write_output(
        &args.output,
        &args.model_coverage,
        &args.bins,
        &emi,
        &header,
        &counts,
    )?;

    //if args.model_coverage {
    //    write_cdf(
    //        &args.output,
    //        &emi,
    //        &header,
    //    )?;
////
    //    write_EquivalenceClass(
    //        &args.output,
    //        &args.model_coverage,
    //        &equivalence_classes,
    //    )?;
    //}
////
    //write_info(
    //    &args.output,
    //    &args.model_coverage,
    //    &emi,
    //    &txps_name,
    //)?;

    Ok(())
}
