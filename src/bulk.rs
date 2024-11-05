use crate::alignment_parser;
use crate::em;
use crate::kde_utils;
use crate::prog_opts::Args;
use crate::util::oarfish_types::AlnInfo;
use crate::util::oarfish_types::DiscardTable;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_output, write_out_prob};
use crate::{binomial_continuous_prob, normalize_read_probs};
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};
use crossbeam::channel::bounded;
use crossbeam::channel::Receiver;
use crossbeam::channel::Sender;
#[allow(unused_imports)]
use minimap2_sys as mm_ffi;
use minimap2_temp as minimap2;
use needletail::parse_fastx_file;
use noodles_bam as bam;
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use std::io::BufRead;
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
        "short_quant": &args.short_quant,
        "num_bootstraps": &args.num_bootstraps,
        "seqcol_digest": seqcol_digest
    })
}

fn perform_inference_and_write_output(
    header: &noodles_sam::header::Header,
    store: &mut InMemoryAlignmentStore,
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

    if args.aln_prob {
        write_out_prob(&args.output, &emi, &txps_name)?;
    }
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
    alignment_parser::parse_alignments(&mut store, header, reader, txps, args.quiet)?;
    perform_inference_and_write_output(header, &mut store, txps, txps_name, seqcol_digest, args)
}

#[derive(Clone)]
struct ReadChunkWithNames {
    read_seq: Vec<u8>,
    read_names: Vec<u8>,
    seq_sep: Vec<usize>,
    name_sep: Vec<usize>,
}

impl ReadChunkWithNames {
    pub fn new() -> Self {
        Self {
            read_seq: Vec::new(),
            read_names: Vec::new(),
            seq_sep: vec![0usize],
            name_sep: vec![0usize],
        }
    }

    #[inline(always)]
    pub fn add_id_and_read(&mut self, id: &[u8], read: &[u8]) {
        self.read_names.extend_from_slice(id);
        self.read_names.push(b'\0');
        self.read_seq.extend_from_slice(read);
        self.name_sep.push(self.read_names.len());
        self.seq_sep.push(self.read_seq.len());
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.read_names.clear();
        self.read_seq.clear();
        self.name_sep.clear();
        self.name_sep.push(0);
        self.seq_sep.clear();
        self.seq_sep.push(0);
    }

    pub fn iter(&self) -> ReadChunkIter {
        ReadChunkIter {
            chunk: self,
            pos: 0,
        }
    }
}

struct ReadChunkIter<'a> {
    chunk: &'a ReadChunkWithNames,
    pos: usize,
}

impl<'a> Iterator for ReadChunkIter<'a> {
    type Item = (&'a [u8], &'a [u8]);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos < self.chunk.seq_sep.len() - 1 {
            let i = self.pos;
            let name: &[u8] =
                &self.chunk.read_names[self.chunk.name_sep[i]..self.chunk.name_sep[i + 1]];
            let seq: &[u8] = &self.chunk.read_seq[self.chunk.seq_sep[i]..self.chunk.seq_sep[i + 1]];
            self.pos += 1;
            Some((name, seq))
        } else {
            None
        }
    }

    #[inline(always)]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let rem = (self.chunk.seq_sep.len() - 1) - self.pos;
        (rem, Some(rem))
    }
}

impl<'a> ExactSizeIterator for ReadChunkIter<'a> {}

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
    let mut store = std::thread::scope(|s| {
        const READ_CHUNK_SIZE: usize = 200;
        const ALN_GROUP_CHUNK_LIMIT: usize = 100;

        let (read_sender, read_receiver): (Sender<ReadGroup>, Receiver<ReadGroup>) =
            bounded(args.threads * 10);
        let (aln_group_sender, aln_group_receiver): (
            Sender<AlignmentGroupInfo>,
            Receiver<AlignmentGroupInfo>,
        ) = bounded(args.threads * 100);

        // Producer thread: reads sequences and sends them to the channel
        let producer = s.spawn(move || {
            let mut ctr = 0_usize;
            let mut chunk_size = 0_usize;
            let mut read_chunk = ReadChunkWithNames::new();

            for read_path in read_paths {
                let mut reader =
                    parse_fastx_file(read_path).expect("valid path/file to read sequences");

                while let Some(result) = reader.next() {
                    let record = result.expect("Error reading record");

                    chunk_size += 1;
                    ctr += 1;

                    // put this read on the current chunk
                    read_chunk.add_id_and_read(record.id(), &record.seq());

                    // send off the next chunks of reads to a thread
                    if chunk_size >= READ_CHUNK_SIZE {
                        read_sender
                            .send(read_chunk.clone())
                            .expect("Error sending sequence");
                        // prepare for the next chunk
                        read_chunk.clear();
                        chunk_size = 0;
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

                    let mut chunk_size = 0_usize;
                    let mut aln_group_alns: Vec<AlnInfo> = Vec::new();
                    let mut aln_group_probs: Vec<f32> = Vec::new();
                    let mut aln_group_boundaries: Vec<usize> = Vec::new();
                    let mut aln_group_read_names = args.aln_prob.then(|| Vec::new());
                    aln_group_boundaries.push(0);

                    // get the next chunk of reads
                    for read_chunk in receiver {
                        // iterate over every read
                        for (name, seq) in read_chunk.iter() {
                            // map the next read, with cigar string
                            let map_res_opt =
                                loc_aligner.map_with_name(name, seq, true, false, None, None);
                            if let Ok(mut mappings) = map_res_opt {
                                if args.aln_prob {
                                    if let Some(first_part) = String::from_utf8_lossy(name).split_whitespace().next() {
                                        let first_part_str = first_part.to_string();
                                        if let Some(last_mapping) = mappings.last_mut() {
                                            last_mapping.query_name = Some(first_part_str);
                                        }
                                    }
                                }
                                
                                let (ag, aprobs, read_name) = filter.filter(
                                    &mut discard_table,
                                    header,
                                    my_txp_info_view,
                                    &mut mappings,
                                );

                                if !ag.is_empty() {
                                    aln_group_alns.extend_from_slice(&ag);
                                    aln_group_probs.extend_from_slice(&aprobs);
                                    aln_group_boundaries.push(aln_group_alns.len());
                                    //eprintln!("read_top: {:?}", read_name);
                                    if let Some(ref mut names_vec) = aln_group_read_names {
                                        if let Some(names) = read_name {
                                            names_vec.push(names);
                                        }
                                    }
                                    chunk_size += 1;
                                }
                                if chunk_size >= ALN_GROUP_CHUNK_LIMIT {
                                    aln_group_sender
                                        .send((
                                            aln_group_alns.clone(),
                                            aln_group_probs.clone(),
                                            aln_group_boundaries.clone(),
                                            aln_group_read_names.clone(),
                                        ))
                                        .expect("Error sending alignment group");
                                    aln_group_alns.clear();
                                    aln_group_probs.clear();
                                    aln_group_boundaries.clear();
                                    aln_group_boundaries.push(0);
                                    if let Some(ref mut names_vec) = aln_group_read_names {
                                        names_vec.clear();
                                    }
                                    chunk_size = 0;
                                }
                            } else {
                                warn!(
                                    "Error encountered mapping read : {}",
                                    map_res_opt.unwrap_err()
                                );
                            }
                        }
                    }
                    if chunk_size > 0 {
                        aln_group_sender
                            .send((aln_group_alns, aln_group_probs, aln_group_boundaries, aln_group_read_names))
                            .expect("Error sending alignment group");
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

            for (ags, aprobs, aln_boundaries, read_name) in aln_group_receiver {
                for (index_receiver, window) in aln_boundaries.windows(2).enumerate() {
                    pb.inc(1);
                    let group_start = window[0];
                    let group_end = window[1];
                    let ag = &ags[group_start..group_end];
                    let as_probs = &aprobs[group_start..group_end];
                    let read_name_slice = read_name.as_ref().map(|names| &names[index_receiver]);

                    //eprintln!("as_probs: {:?}", as_probs);
                    //eprintln!("read_name_slice: {:?}", read_name_slice);

                    if ag.len() == 1 {
                        store.inc_unique_alignments();
                    }
                    store.add_filtered_group(ag, as_probs, read_name_slice.cloned(), txps_mut);
                }
            }
            pb.finish_with_message("Finished aligning reads.");
            store
        });

        // Wait for the producer to finish reading
        let total_reads = producer.join().expect("Producer thread panicked");

        let mut discard_tables: Vec<DiscardTable> = Vec::with_capacity(map_threads);
        for consumer in consumers {
            let dt = consumer.join().expect("Consumer thread panicked");
            discard_tables.push(dt);
        }

        drop(aln_group_sender);

        let mut store = aln_group_consumer
            .join()
            .expect("Alignment group consumer panicked");

        info!(
            "Parsed {} total reads",
            total_reads.to_formatted_string(&Locale::en)
        );

        for dt in &discard_tables {
            store.aggregate_discard_table(dt);
        }
        store
    });

    perform_inference_and_write_output(header, &mut store, txps, txps_name, seqcol_digest, args)
}
