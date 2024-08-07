use clap::Parser;
use core::str;
use std::num::NonZeroUsize;
use util::write_function;

use anyhow::Context;
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};

use std::{
    fs::{create_dir_all, File},
    io::{self, Write},
    path::{Path, PathBuf},
};

use num_format::{Locale, ToFormattedString};
use serde::Serialize;
use serde_json::json;
use tracing::{error, info};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

use path_tools::WithAdditionalExtension;

use noodles_bam as bam;
use noodles_bgzf as bgzf;

mod alignment_parser;
mod bootstrap;
mod em;
mod util;

use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_output};
use crate::util::{binomial_probability::binomial_continuous_prob, kde_utils};

/// These represent different "meta-options", specific settings
/// for all of the different filters that should be applied in
/// different cases.
#[derive(Clone, Debug, clap::ValueEnum, Serialize)]
enum FilterGroup {
    NoFilters,
    NanocountFilters,
}

fn parse_strand(arg: &str) -> anyhow::Result<bio_types::strand::Strand> {
    match arg {
        "+" | "fw" | "FW" | "f" | "F" => Ok(bio_types::strand::Strand::Forward),
        "-" | "rc" | "RC" | "r" | "R" => Ok(bio_types::strand::Strand::Reverse),
        "." | "both" | "either" => Ok(bio_types::strand::Strand::Unknown),
        _ => anyhow::bail!("Cannot parse {} as a valid strand type", arg),
    }
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
    /// only alignments to this strand will be allowed; options are (fw /+, rc/-, or both/.)
    #[arg(
        short = 'd',
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = bio_types::strand::Strand::Unknown,
        value_parser = parse_strand
    )]
    strand_filter: bio_types::strand::Strand,
    /// input is assumed to be a single-cell BAM and to have the `CB:z` tag for all read records
    #[arg(long)]
    single_cell: bool,
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
    /// width of the bins used in the coverage model
    #[arg(short, long, help_heading = "coverage model", default_value_t = 100)]
    bin_width: u32,
    /// use a KDE model of the observed fragment length distribution
    #[arg(short, long, hide = true)]
    use_kde: bool,
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

/// Produce a [serde_json::Value] that encodes the relevant arguments and
/// parameters of the run that we wish to record to file. Ultimately, this
/// will be written to the corresponding `meta_info.json` file for this run.
fn get_single_cell_json_info(args: &Args, seqcol_digest: &str) -> serde_json::Value {
    let prob = if args.model_coverage {
        "scaled_binomial"
    } else {
        "no_coverage"
    };

    json!({
        "prob_model" : prob,
        "bin_width" : args.bin_width,
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
        "seqcol_digest": seqcol_digest
    })
}

use noodles_sam::alignment::record::data::field::Value;

pub fn get_while<'a, R: io::BufRead>(
    filter_opts: &AlignmentFilters,
    header: &'a noodles_sam::Header,
    txps: &mut [TranscriptInfo],
    records_for_read: &mut Vec<noodles_sam::alignment::record_buf::RecordBuf>,
    iter: &mut core::iter::Peekable<noodles_bam::io::reader::RecordBufs<R>>,
    barcode: &[u8],
) -> anyhow::Result<InMemoryAlignmentStore<'a>> {
    let mut astore = InMemoryAlignmentStore::new(filter_opts.clone(), header);
    let num_rec_processed = alignment_parser::parse_alignments_for_barcode(
        &mut astore,
        txps,
        iter,
        barcode,
        records_for_read,
    )?;
    Ok(astore)
}

use crossbeam::queue::ArrayQueue;
use std::sync::Arc;

struct QuantOutputInfo {
    barcode_file: std::io::BufWriter<File>,
    row_ids: Vec<u32>,
    col_ids: Vec<u32>,
    vals: Vec<f32>,
    row_index: usize,
}

pub fn quantify_single_cell_from_collated_bam<R: io::BufRead>(
    header: &noodles_sam::Header,
    filter_opts: &AlignmentFilters,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    args: &Args,
    seqcol_digest: String,
) -> anyhow::Result<()> {
    // if there is a parent directory
    if let Some(p) = args.output.parent() {
        // unless this was a relative path with one component,
        // which we should treat as the file prefix, then grab
        // the non-empty parent and create it.
        if p != Path::new("") {
            create_dir_all(p)?;
        }
    }

    std::thread::scope(|s| {
        let bc_path = args.output.with_additional_extension(".barcodes.txt");
        let bc_file = File::create(bc_path)?;
        let bc_writer = Arc::new(std::sync::Mutex::new(QuantOutputInfo {
            barcode_file: std::io::BufWriter::new(bc_file),
            row_ids: Vec::new(),
            col_ids: Vec::new(),
            vals: Vec::new(),
            row_index: 0usize,
        }));

        // get the data for the next cell
        let mut records_for_read: Vec<noodles_sam::alignment::RecordBuf> = Vec::new();
        let mut peekable_bam_iter = reader.record_bufs(header).peekable();
        let nthreads = args.threads;
        const CB_TAG: [u8; 2] = [b'C', b'B'];

        let q: Arc<ArrayQueue<Arc<(InMemoryAlignmentStore, Vec<u8>)>>> =
            Arc::new(ArrayQueue::new(4 * nthreads));
        let done_parsing = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let mut thread_handles: Vec<std::thread::ScopedJoinHandle<'_, anyhow::Result<usize>>> =
            Vec::with_capacity(nthreads);

        for _worker_id in 0..nthreads {
            let in_q = q.clone();
            let done_parsing = done_parsing.clone();
            let mut new_txps = Vec::with_capacity(txps.len());
            new_txps.extend_from_slice(txps);
            let bc_out = bc_writer.clone();

            let handle = s.spawn(move || {
                let mut col_ids = Vec::with_capacity(new_txps.len());
                let mut row_ids = Vec::with_capacity(new_txps.len());
                let mut vals = Vec::with_capacity(new_txps.len());
                let mut num_cells = 0_usize;

                while !done_parsing.load(std::sync::atomic::Ordering::SeqCst) {
                    while let Some(astore) = in_q.pop() {
                        //println!("num_rec = {}", astore.len());
                        // wrap up all of the relevant information we need for estimation
                        // in an EMInfo struct and then call the EM algorithm.
                        let emi = EMInfo {
                            eq_map: &astore.0,
                            txp_info: &new_txps,
                            max_iter: 1000,
                            convergence_thresh: 0.001,
                            init_abundances: None,
                            kde_model: None,
                        };
                        let counts = em::em(&emi, 1);
                        col_ids.clear();
                        vals.clear();
                        for (col_idx, v) in counts.iter().enumerate() {
                            if *v > 0.0 {
                                col_ids.push(col_idx as u32);
                                vals.push((*v) as f32);
                            }
                        }
                        row_ids.resize(col_ids.len(), 0_u32);
                        num_cells += 1;

                        let row_index: usize;
                        {
                            let writer_deref = bc_out.lock();
                            let writer = &mut *writer_deref.unwrap();
                            writeln!(&mut writer.barcode_file, "{}", unsafe {
                                std::str::from_utf8_unchecked(&astore.1)
                            })?;

                            // get the row index and then increment it
                            row_index = writer.row_index;
                            writer.row_index += 1;
                            row_ids.fill(row_index as u32);

                            writer.col_ids.extend_from_slice(&col_ids);
                            writer.row_ids.extend_from_slice(&row_ids);
                            writer.vals.extend_from_slice(&vals);
                        }
                    }
                }
                Ok(num_cells)
            });
            thread_handles.push(handle);
        }

        // parser thread
        while let Some(next_res) = peekable_bam_iter.peek() {
            let rec = next_res.as_ref().unwrap();
            let barcode = match rec.data().get(&CB_TAG) {
                None => anyhow::bail!("could not get CB tag value"),
                Some(v) => match v {
                    noodles_sam::alignment::record_buf::data::field::Value::String(x) => {
                        x.to_ascii_uppercase()
                    }
                    _ => anyhow::bail!("CB tag value had unexpected type!"),
                },
            };

            let mut astore = std::sync::Arc::new((
                get_while(
                    filter_opts,
                    header,
                    txps,
                    &mut records_for_read,
                    &mut peekable_bam_iter,
                    &barcode,
                )?,
                barcode,
            ));

            // push the store on to the work queue
            while let Err(store) = q.push(astore) {
                astore = store;
                while q.is_full() {}
            }
        }

        done_parsing.store(true, std::sync::atomic::Ordering::SeqCst);

        let mut total_cells = 0_usize;
        for h in thread_handles {
            match h.join() {
                Ok(Ok(nc)) => {
                    total_cells += nc;
                }
                Ok(Err(e)) => {
                    error!("error result from thread {:?}", e);
                }
                Err(_e) => {
                    error!("thread panicked");
                }
            }
        }

        let trimat = {
            let writer_deref = bc_writer.lock();
            let writer = &mut *writer_deref.unwrap();
            let num_rows = total_cells;
            sprs::TriMatI::<f32, u32>::from_triplets(
                (num_rows, txps.len()),
                writer.row_ids.clone(),
                writer.col_ids.clone(),
                writer.vals.clone(),
            )
        };
        let info = get_single_cell_json_info(&args, &seqcol_digest);
        write_function::write_single_cell_output(&args.output, info, &header, &trimat)?;
        Ok(())
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
        quantify_single_cell_from_collated_bam(
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
