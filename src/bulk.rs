use crate::alignment_parser;
use crate::em;
use crate::kde_utils;
use crate::prog_opts::Args;
use crate::util::oarfish_types::AlnInfo;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_output};
use crate::{binomial_continuous_prob, normalize_read_probs};
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};
#[allow(unused_imports)]
use minimap2_sys as mm_ffi;
use noodles_bam as bam;
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use std::io::BufRead;
use tracing::info;

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

pub fn quantify_bulk_alignments_from_bam<R: BufRead>(
    header: &noodles_sam::Header,
    filter_opts: AlignmentFilters,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: String,
) -> anyhow::Result<()> {
    // now parse the actual alignments for the reads and store the results
    // in our in-memory stor
    let mut store = InMemoryAlignmentStore::new(filter_opts, header);
    alignment_parser::parse_alignments(&mut store, header, reader, txps)?;

    // print discard table information in which the user might be interested.
    info!("\ndiscard_table: \n{}\n", store.discard_table.to_table());

    // if we are using the KDE, create that here.
    let kde_opt: Option<kders::kde::KDEModel> = if args.use_kde {
        Some(kde_utils::get_kde_model(txps, &store)?)
    } else {
        None
    };

    if store.filter_opts.model_coverage {
        //obtaining the Cumulative Distribution Function (CDF) for each transcript
        binomial_continuous_prob(txps, &args.bin_width, args.threads);
        //Normalize the probabilities for the records of each read
        normalize_read_probs(&mut store, txps, &args.bin_width);
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
        eq_map: &store,
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

    let aux_txp_counts = crate::util::aux_counts::get_aux_counts(&store, txps)?;

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
    Ok(())
}

use crate::util::oarfish_types::DiscardTable;
use crossbeam::channel::bounded;
use crossbeam::channel::Receiver;
use crossbeam::channel::Sender;
use needletail::parse_fastx_file;

#[allow(clippy::too_many_arguments)]
pub fn quantify_bulk_alignments_raw_reads(
    header: &noodles_sam::Header,
    aligner: minimap2::Aligner,
    filter_opts: AlignmentFilters,
    read_path: std::path::PathBuf,
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

    type AlignmentGroupInfo = (Vec<AlnInfo>, Vec<f32>);
    let mut store = std::thread::scope(|s| {
        const READ_CHUNK_SIZE: usize = 100;

        let (read_sender, read_receiver): (
            Sender<(Vec<u8>, Vec<usize>)>,
            Receiver<(Vec<u8>, Vec<usize>)>,
        ) = bounded(args.threads * 10);
        let (aln_group_sender, aln_group_receiver): (
            Sender<AlignmentGroupInfo>,
            Receiver<AlignmentGroupInfo>,
        ) = bounded(args.threads * 10);

        // Producer thread: reads sequences and sends them to the channel
        let producer = s.spawn(move || {
            let mut reader = parse_fastx_file(read_path).expect("valid path/file");
            let mut ctr = 0_usize;
            let mut chunk_size = 0_usize;
            let mut reads_vec: Vec<u8> = Vec::new();
            let mut read_boundaries: Vec<usize> = Vec::new();
            read_boundaries.push(0);

            while let Some(result) = reader.next() {
                let record = result.expect("Error reading record");

                chunk_size += 1;
                ctr += 1;

                // put this read on the current chunk
                reads_vec.extend_from_slice(&record.seq());
                read_boundaries.push(reads_vec.len());

                // send off the next chunks of reads to a thread
                if chunk_size >= READ_CHUNK_SIZE {
                    read_sender
                        .send((reads_vec.clone(), read_boundaries.clone()))
                        .expect("Error sending sequence");
                    // prepare for the next chunk
                    reads_vec.clear();
                    read_boundaries.clear();
                    read_boundaries.push(0);
                    chunk_size = 0;
                }
            }

            // if any reads remain, send them off
            if chunk_size > 0 {
                read_sender
                    .send((reads_vec, read_boundaries))
                    .expect("Error sending sequence");
            }
            ctr
        });

        // Consumer threads: receive sequences and perform alignment
        let consumers: Vec<_> = (0..map_threads)
            .map(|_| {
                let receiver = read_receiver.clone();
                let mut filter = filter_opts.clone();
                let mut loc_aligner = aligner.clone().with_cigar();

                let my_txp_info_view = &txp_info_view;
                let aln_group_sender = aln_group_sender.clone();
                s.spawn(move || {
                    let mut discard_table = DiscardTable::new();
                    // get the next chunk of reads
                    for (seq, boundaries) in receiver {
                        // iterate over every read
                        for window in boundaries.windows(2) {
                            let map_res_opt = loc_aligner.map(
                                &seq[window[0]..window[1]],
                                true,
                                false,
                                None,
                                None,
                            );
                            let map_res = map_res_opt.ok();
                            if let Some(mut mappings) = map_res {
                                let (ag, aprobs) = filter.filter(
                                    &mut discard_table,
                                    header,
                                    my_txp_info_view,
                                    &mut mappings,
                                );
                                aln_group_sender
                                    .send((ag, aprobs))
                                    .expect("Error sending alignment group");
                            }
                        }
                    }
                    // NOTE: because `clone()` clones the raw pointer, if it is
                    // still set when this tread goes out of scope, the underlying
                    // raw pointer will be freed and the other aligners will have
                    // references to already freed memory and this will lead to a
                    // double-free when they are dropped. So, to avoid this, here
                    // we set the idx pointer to None directly. Track: https://github.com/jguhlin/minimap2-rs/issues/71
                    loc_aligner.idx = None;
                    discard_table
                })
            })
            .collect();

        #[allow(clippy::useless_asref)]
        let txps_mut = txps.as_mut();
        let filter_opts_store = filter_opts.clone();
        let aln_group_consumer = s.spawn(move || {
            let mut store = InMemoryAlignmentStore::new(filter_opts_store, header);
            for (ng, (ag, as_probs)) in aln_group_receiver.iter().enumerate() {
                if ng > 0 && (ng % 100000 == 1) {
                    info!("processed {} mapped reads", ng);
                }
                if ag.len() == 1 {
                    store.inc_unique_alignments();
                }
                store.add_filtered_group(&ag, &as_probs, txps_mut);
            }
            info!("DONE!");
            store
        });

        // Wait for the producer to finish reading
        let total_reads = producer.join().expect("Producer thread panicked");
        info!("Read Producer finished; parsed {} reads", total_reads);

        let mut discard_tables: Vec<DiscardTable> = Vec::with_capacity(map_threads);
        for consumer in consumers {
            let dt = consumer.join().expect("Consumer thread panicked");
            discard_tables.push(dt);
        }

        drop(aln_group_sender);

        let mut store = aln_group_consumer
            .join()
            .expect("Alignment group consumer panicked");

        for dt in &discard_tables {
            store.aggregate_discard_table(dt);
        }
        store
    });
    // print discard table information in which the user might be interested.
    info!("\ndiscard_table: \n{}\n", store.discard_table.to_table());

    // if we are using the KDE, create that here.
    let kde_opt: Option<kders::kde::KDEModel> = if args.use_kde {
        Some(kde_utils::get_kde_model(txps, &store)?)
    } else {
        None
    };

    if store.filter_opts.model_coverage {
        //obtaining the Cumulative Distribution Function (CDF) for each transcript
        binomial_continuous_prob(txps, &args.bin_width, args.threads);
        //Normalize the probabilities for the records of each read
        normalize_read_probs(&mut store, txps, &args.bin_width);
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
        eq_map: &store,
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

    let aux_txp_counts = crate::util::aux_counts::get_aux_counts(&store, txps)?;

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
    Ok(())
}
