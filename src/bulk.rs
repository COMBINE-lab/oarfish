use crate::alignment_parser;
use crate::em;
use crate::kde_utils;
use crate::prog_opts::Args;
use crate::util::constants::EMPTY_READ_NAME;
use crate::util::oarfish_types::AlnInfo;
use crate::util::oarfish_types::DiscardTable;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, InputSourceType, ReadChunkWithNames,
    ReadSource, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_out_prob, write_output};
use crate::{binomial_continuous_prob, normalize_read_probs};
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};
use crossbeam::channel::bounded;
use crossbeam::channel::Receiver;
use crossbeam::channel::Sender;
#[allow(unused_imports)]
use minimap2_sys as mm_ffi;
//use minimap2_temp as minimap2;

use needletail::parse_fastx_file;
use noodles_bam as bam;
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use std::io::BufRead;
use swapvec::{SwapVec, SwapVecConfig};
use tracing::{info, warn};

/// Produce a [serde_json::Value] that encodes the relevant arguments and
/// parameters of the run that we wish to record to file. Ultimately, this
/// will be written to the corresponding `meta_info.json` file for this run.
fn get_json_info(args: &Args, emi: &EMInfo, seqcol_digest: &str) -> serde_json::Value {
    let prob = if args.model_coverage {
        "scaled_binomial"
    } else {
        "no_coverage"
    };

    let source = if args.alignments.is_some() {
        "from_bam"
    } else {
        "from_raw_reads"
    };

    json!({
        "prob_model" : prob,
        "alignment_source" : source,
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
        "write_assignment_probs": &emi.eq_map.filter_opts.write_assignment_probs_type,
        "short_quant": &args.short_quant,
        "num_bootstraps": &args.num_bootstraps,
        "seqcol_digest": seqcol_digest
    })
}

fn perform_inference_and_write_output(
    header: &noodles_sam::header::Header,
    store: &mut InMemoryAlignmentStore,
    name_vec: Option<SwapVec<String>>,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    seqcol_digest: String,
    args: &Args,
) -> anyhow::Result<()> {
    // print discard table information in which the user might be interested.
    info!("\ndiscard_table: \n{}\n", store.discard_table.to_table());

    // if we are using the KDE, create that here.
    let kde_opt: Option<kders::kde::KDEModel> = if args.use_kde {
        Some(kde_utils::get_kde_model(txps, store)?)
    } else {
        None
    };

    if store.filter_opts.model_coverage {
        //obtaining the Cumulative Distribution Function (CDF) for each transcript
        binomial_continuous_prob(txps, &args.bin_width, args.threads);
        //Normalize the probabilities for the records of each read
        normalize_read_probs(store, txps, &args.bin_width);
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
        read_short_quant_vec(sr_path, txps_name).unwrap_or_else(|e| panic!("{}", e))
    });

    // wrap up all of the relevant information we need for estimation
    // in an EMInfo struct and then call the EM algorithm.
    let emi = EMInfo {
        eq_map: store,
        txp_info: txps,
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

    let aux_txp_counts = crate::util::aux_counts::get_aux_counts(store, txps)?;

    // prepare the JSON object we'll write
    // to meta_info.json
    let json_info = get_json_info(args, &emi, &seqcol_digest);

    // write the output
    write_output(&args.output, json_info, header, &counts, &aux_txp_counts)?;

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

    if args.write_assignment_probs.is_some() {
        let name_vec = name_vec
            .expect("cannot write assignment probabilities without valid vector of read names");
        write_out_prob(&args.output, &emi, name_vec, txps_name)?;
    }

    Ok(())
}

pub fn quantify_bulk_alignments_from_bam<R: BufRead>(
    header: &noodles_sam::Header,
    filter_opts: AlignmentFilters,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: String,
) -> anyhow::Result<()> {
    let mut name_vec = if filter_opts.write_assignment_probs {
        Some(SwapVec::<String>::with_config(SwapVecConfig {
            swap_after: Default::default(),
            batch_size: Default::default(),
            compression: Some(swapvec::Compression::Lz4),
        }))
    } else {
        None
    };
    // now parse the actual alignments for the reads and store the results
    // in our in-memory stor
    let mut store = InMemoryAlignmentStore::new(filter_opts, header);
    alignment_parser::parse_alignments(
        &mut store,
        &mut name_vec,
        header,
        reader,
        txps,
        args.sort_check_num,
        args.quiet,
    )?;
    perform_inference_and_write_output(
        header,
        &mut store,
        name_vec,
        txps,
        txps_name,
        seqcol_digest,
        args,
    )
}

fn get_source_type(pb: &std::path::Path) -> InputSourceType {
    let faq_endings = vec![
        ".fasta",
        ".fastq",
        ".FASTA",
        ".FASTQ",
        ".fa",
        ".fq",
        ".FA",
        ".FQ",
        ".fasta.gz",
        ".fastq.gz",
        ".FASTA.GZ",
        ".FASTQ.GZ",
        ".fa.gz",
        ".fq.gz",
        ".FA.GZ",
        ".FQ.GZ",
    ];
    let ubam_endings = vec![".bam", ".BAM", ".ubam", ".UBAM"];
    if let Some(ps) = pb.to_str() {
        for fe in faq_endings {
            if ps.ends_with(fe) {
                return InputSourceType::Fastx;
            }
        }
        for be in ubam_endings {
            if ps.ends_with(be) {
                return InputSourceType::Ubam;
            }
        }
    }

    InputSourceType::Unknown
}

#[allow(clippy::too_many_arguments)]
pub fn quantify_bulk_alignments_raw_reads(
    header: &noodles_sam::Header,
    aligner: minimap2::Aligner,
    filter_opts: AlignmentFilters,
    read_paths: &[std::path::PathBuf],
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: String,
) -> anyhow::Result<()> {
    // now parse the actual alignments for the reads and store the results
    // in our in-memory stor

    // we will take only shared refs to this
    let mut txp_info_view: Vec<TranscriptInfo> = Vec::with_capacity(txps.len());
    for ti in txps.iter() {
        txp_info_view.push(ti.clone());
    }

    // at least one mapping thread, otherwise everything but the fastx parser
    // and the in memory alignment store populator
    let map_threads = args.threads.saturating_sub(2).max(1);

    type ReadGroup = ReadChunkWithNames;
    type AlignmentGroupInfo = (Vec<AlnInfo>, Vec<f32>, Vec<usize>, Option<Vec<String>>);

    let (read_sender, read_receiver): (Sender<ReadGroup>, Receiver<ReadGroup>) =
        bounded(args.threads * 10);

    const READ_CHUNK_SIZE: usize = 200;
    let mut rpaths = vec![];
    read_paths.clone_into(&mut rpaths);

    // Producer thread: reads sequences and sends them to the channel
    let producer = std::thread::spawn(move || {
        let mut ctr = 0_usize;
        let mut chunk_size = 0_usize;
        let mut read_chunk = ReadChunkWithNames::new();

        // work shared between the two different
        // source types
        let mark_chunk = |chunk_size: &mut usize,
                          ctr: &mut usize,
                          read_chunk: &mut ReadGroup,
                          read_sender: &Sender<ReadGroup>| {
            *chunk_size += 1;
            *ctr += 1;
            if *chunk_size >= READ_CHUNK_SIZE {
                read_sender
                    .send(read_chunk.clone())
                    .expect("Error sending sequence");
                // prepare for the next chunk
                read_chunk.clear();
                *chunk_size = 0;
            }
        };

        // read from either a UBAM or (possibly compressed) FASTX file
        for read_path in rpaths {
            match get_source_type(&read_path) {
                InputSourceType::Ubam => {
                    let mut reader = std::fs::File::open(read_path)
                        .map(bam::io::Reader::new)
                        .expect("could not create BAM reader");
                    let header = reader.read_header().expect("could not read BAM header");
                    for result in reader.record_bufs(&header) {
                        let record = result.expect("Error reading ubam record");
                        record.add_to_read_group(&mut read_chunk);
                        mark_chunk(&mut chunk_size, &mut ctr, &mut read_chunk, &read_sender);
                    }
                }
                s @ (InputSourceType::Fastx | InputSourceType::Unknown) => {
                    if matches!(s, InputSourceType::Unknown) {
                        warn!("could not determine input file type for {} from suffix; assuming (possibly gzipped) fastx", &read_path.display());
                    }
                    let mut reader =
                        parse_fastx_file(read_path).expect("valid path/file to read sequences");
                    while let Some(result) = reader.next() {
                        let record = result.expect("Error reading record");
                        record.add_to_read_group(&mut read_chunk);
                        mark_chunk(&mut chunk_size, &mut ctr, &mut read_chunk, &read_sender);
                    }
                }
            }
        }
        // if any reads remain, send them off
        if chunk_size > 0 {
            read_sender
                .send(read_chunk)
                .expect("Error sending sequence");
        }
        ctr
    });

    // we need the scope here so we can borrow the relevant non-'static data
    let (mut store, name_vec) = std::thread::scope(|s| {
        const ALN_GROUP_CHUNK_LIMIT: usize = 100;

        let (aln_group_sender, aln_group_receiver): (
            Sender<AlignmentGroupInfo>,
            Receiver<AlignmentGroupInfo>,
        ) = bounded(args.threads * 100);

        // Consumer threads: receive sequences and perform alignment
        let write_assignment_probs: bool = args.write_assignment_probs.is_some();
        let consumers: Vec<_> = (0..map_threads)
            .map(|_| {
                let receiver = read_receiver.clone();
                let mut filter = filter_opts.clone();
                let loc_aligner = aligner.clone();

                let my_txp_info_view = &txp_info_view;
                let aln_group_sender = aln_group_sender.clone();
                s.spawn(move || {
                    let mut discard_table = DiscardTable::new();

                    let mut chunk_size = 0_usize;
                    let mut aln_group_alns: Vec<AlnInfo> = Vec::new();
                    let mut aln_group_probs: Vec<f32> = Vec::new();
                    let mut aln_group_boundaries: Vec<usize> = Vec::new();
                    let mut aln_group_read_names = write_assignment_probs.then(Vec::new);
                    aln_group_boundaries.push(0);

                    // get the next chunk of reads
                    for read_chunk in receiver {
                        // iterate over every read
                        for (name, seq) in read_chunk.iter() {
                            // map the next read, with cigar string
                            let map_res_opt =
                                loc_aligner.map(seq, true, false, None, None, Some(name));
                            if let Ok(mut mappings) = map_res_opt {
                                let (ag, aprobs) = filter.filter(
                                    &mut discard_table,
                                    header,
                                    my_txp_info_view,
                                    &mut mappings,
                                );

                                if !ag.is_empty() {
                                    aln_group_alns.extend_from_slice(&ag);
                                    aln_group_probs.extend_from_slice(&aprobs);
                                    aln_group_boundaries.push(aln_group_alns.len());
                                    // if we are storing read names
                                    if let Some(ref mut names_vec) = aln_group_read_names {
                                        let name_str = String::from_utf8_lossy(name).into_owned();
                                        names_vec.push(name_str);
                                    }
                                    chunk_size += 1;
                                }
                                if chunk_size >= ALN_GROUP_CHUNK_LIMIT {
                                    aln_group_sender
                                        .send((
                                            aln_group_alns.clone(),
                                            aln_group_probs.clone(),
                                            aln_group_boundaries.clone(),
                                            aln_group_read_names,
                                        ))
                                        .expect("Error sending alignment group");
                                    aln_group_alns.clear();
                                    aln_group_probs.clear();
                                    aln_group_boundaries.clear();
                                    aln_group_boundaries.push(0);
                                    aln_group_read_names = write_assignment_probs.then(Vec::new);
                                    chunk_size = 0;
                                }
                            } else {
                                warn!(
                                    "Error encountered mappread_ing read : {}",
                                    map_res_opt.unwrap_err()
                                );
                            }
                        }
                    }
                    if chunk_size > 0 {
                        aln_group_sender
                            .send((
                                aln_group_alns,
                                aln_group_probs,
                                aln_group_boundaries,
                                aln_group_read_names,
                            ))
                            .expect("Error sending alignment group");
                    }
                    discard_table
                })
            })
            .collect();

        #[allow(clippy::useless_asref)]
        let txps_mut = txps.as_mut();
        let filter_opts_store = filter_opts.clone();
        let aln_group_consumer = s.spawn(move || {
            let mut name_vec = if filter_opts_store.write_assignment_probs {
                Some(SwapVec::<String>::with_config(SwapVecConfig {
                    swap_after: Default::default(),
                    batch_size: Default::default(),
                    compression: Some(swapvec::Compression::Lz4),
                }))
            } else {
                None
            };

            let mut store = InMemoryAlignmentStore::new(filter_opts_store, header);

            let pb = if args.quiet {
                indicatif::ProgressBar::hidden()
            } else {
                indicatif::ProgressBar::new_spinner().with_message("Number of reads mapped")
            };

            pb.set_style(
                indicatif::ProgressStyle::with_template(
                    "[{elapsed_precise}] {spinner:4.green/blue} {msg} {human_pos:>12}",
                )
                .unwrap()
                .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
            );
            pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

            for (ags, aprobs, aln_boundaries, read_names) in aln_group_receiver {
                // if we are getting read names out then we are going to "reverse" them
                // here so that we can simply pop the strings off the back to get them
                // in order. We do this since we cannot otherwise "move" a string out of a
                // Vec.
                let mut reversed_read_names = if let Some(mut names_vec) = read_names {
                    names_vec.reverse();
                    Some(names_vec)
                } else {
                    None
                };

                for window in aln_boundaries.windows(2) {
                    pb.inc(1);
                    let group_start = window[0];
                    let group_end = window[1];
                    let ag = &ags[group_start..group_end];
                    let as_probs = &aprobs[group_start..group_end];
                    let read_name_opt = if let Some(ref mut names_vec) = reversed_read_names {
                        names_vec.pop()
                    } else {
                        None
                    };

                    if store.add_filtered_group(ag, as_probs, txps_mut) {
                        if let Some(ref mut nvec) = name_vec {
                            let read_name = read_name_opt.unwrap_or(EMPTY_READ_NAME.to_string());
                            nvec.push(read_name)
                                .expect("cannot push name to read name vector");
                        }
                        if ag.len() == 1 {
                            store.inc_unique_alignments();
                        }
                    }
                }
            }
            pb.finish_with_message("Finished aligning reads.");
            (store, name_vec)
        });

        // Wait for the producer to finish reading
        let total_reads = producer.join().expect("Producer thread panicked");

        let mut discard_tables: Vec<DiscardTable> = Vec::with_capacity(map_threads);
        for consumer in consumers {
            let dt = consumer.join().expect("Consumer thread panicked");
            discard_tables.push(dt);
        }

        drop(aln_group_sender);

        let (mut store, name_vec) = aln_group_consumer
            .join()
            .expect("Alignment group consumer panicked");

        info!(
            "Parsed {} total reads",
            total_reads.to_formatted_string(&Locale::en)
        );

        for dt in &discard_tables {
            store.aggregate_discard_table(dt);
        }
        (store, name_vec)
    });

    perform_inference_and_write_output(
        header,
        &mut store,
        name_vec,
        txps,
        txps_name,
        seqcol_digest,
        args,
    )
}
