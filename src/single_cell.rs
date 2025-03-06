use crate::alignment_parser;
use crate::em;
use crate::prog_opts::Args;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::write_function;
use crossbeam::queue::ArrayQueue;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use path_tools::WithAdditionalExtension;
use serde_json::json;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::{error, info};

struct QuantOutputInfo {
    barcode_file: std::io::BufWriter<File>,
    row_ids: Vec<u32>,
    col_ids: Vec<u32>,
    vals: Vec<f32>,
    row_index: usize,
}

/// Produce a [serde_json::Value] that encodes the relevant arguments and
/// parameters of the run that we wish to record to file. Ultimately, this
/// will be written to the corresponding `meta_info.json` file for this run.
fn get_single_cell_json_info(
    args: &Args,
    seqcol_digest: &seqcol_rs::DigestResult,
) -> serde_json::Value {
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
        "digest": seqcol_digest.to_json()
    })
}

pub fn quantify_single_cell_from_collated_bam<R: BufRead>(
    header: &noodles_sam::Header,
    filter_opts: &AlignmentFilters,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    args: &Args,
    seqcol_digest: seqcol_rs::DigestResult,
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

    let nthreads = args.threads;
    std::thread::scope(|s| {
        let bc_path = args.output.with_additional_extension(".barcodes.txt");
        let bc_file = File::create(bc_path)?;
        let bc_writer = Arc::new(Mutex::new(QuantOutputInfo {
            barcode_file: std::io::BufWriter::new(bc_file),
            row_ids: Vec::new(),
            col_ids: Vec::new(),
            vals: Vec::new(),
            row_index: 0usize,
        }));

        // the element consists of the vector of records corresponding
        // to this cell, a (read-only) copy of TranscriptInfo which will
        // be copied and modified by the thread, and the barcode
        // (represented as a Vec<u8>) for the cell.
        type QueueElement<'a> = (Vec<RecordBuf>, &'a [TranscriptInfo], Vec<u8>);

        let q: Arc<ArrayQueue<QueueElement>> = Arc::new(ArrayQueue::new(4 * nthreads));
        let done_parsing = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let mut thread_handles: Vec<std::thread::ScopedJoinHandle<'_, anyhow::Result<usize>>> =
            Vec::with_capacity(nthreads);

        for _worker_id in 0..nthreads {
            let in_q = q.clone();
            let done_parsing = done_parsing.clone();
            let num_txps = txps.len();
            let bc_out = bc_writer.clone();
            let bin_width = args.bin_width;
            let filter_opts = filter_opts.clone();

            let handle = s.spawn(move || {
                let mut col_ids = Vec::with_capacity(num_txps);
                let mut row_ids = Vec::with_capacity(num_txps);
                let mut vals = Vec::with_capacity(num_txps);
                let mut num_cells = 0_usize;
                let mut records_for_read = Vec::<RecordBuf>::with_capacity(16);

                // while the queue might still be being filled
                while !done_parsing.load(std::sync::atomic::Ordering::SeqCst) {
                    // get the next cell
                    while let Some(elem) = in_q.pop() {
                        let mut recs = elem.0;
                        // new copy of txp info for this barcode
                        let mut txps = Vec::with_capacity(num_txps);
                        txps.extend_from_slice(elem.1);
                        // the barcode of this cell
                        let barcode = elem.2;
                        // where we will store the relevant alignment records
                        let mut store = InMemoryAlignmentStore::new(filter_opts.clone(), header);

                        // sort by read name and then parse the records for this cell
                        alignment_parser::sort_and_parse_barcode_records(
                            &mut recs,
                            &mut store,
                            &mut txps,
                            &mut records_for_read,
                        )?;

                        if store.filter_opts.model_coverage {
                            //obtaining the Cumulative Distribution Function (CDF) for each transcript
                            crate::binomial_continuous_prob(&mut txps, &bin_width, 1);
                            //Normalize the probabilities for the records of each read
                            crate::normalize_read_probs(&mut store, &txps, &bin_width);
                        }

                        // wrap up all of the relevant information we need for estimation
                        // in an EMInfo struct and then call the EM algorithm.
                        let emi = EMInfo {
                            eq_map: &store,
                            txp_info: &txps,
                            max_iter: args.max_em_iter,
                            convergence_thresh: args.convergence_thresh,
                            init_abundances: None,
                            kde_model: None,
                        };
                        // run the EM for this cell
                        let counts = em::em(&emi, 1);
                        // clear out the vectors where we will store
                        // the count information for this cell
                        col_ids.clear();
                        vals.clear();
                        for (col_idx, v) in counts.iter().enumerate() {
                            if *v > 0.0 {
                                col_ids.push(col_idx as u32);
                                vals.push((*v) as f32);
                            }
                        }
                        // fill the row ids for this cell; fist
                        // we size the vector to the correct length
                        // and fill it with 0s and below we
                        // fill with the appropriate number (i.e. the
                        // cell/barcode ID).
                        row_ids.resize(col_ids.len(), 0_u32);
                        num_cells += 1;

                        let row_index: usize;
                        {
                            // grab a lock and fill out the count info for
                            // this cell.
                            let writer_deref = bc_out.lock();
                            let writer = &mut *writer_deref.unwrap();
                            writeln!(&mut writer.barcode_file, "{}", unsafe {
                                std::str::from_utf8_unchecked(&barcode)
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

        // get the data for the next cell
        let mut peekable_bam_iter = reader.record_bufs(header).peekable();
        const CB_TAG: [u8; 2] = [b'C', b'B'];
        let mut num_cells = 0_usize;
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

            let records_for_barcode =
                alignment_parser::parse_alignments_for_barcode(&mut peekable_bam_iter, &barcode)?;

            num_cells += 1;
            if num_cells > 1 && num_cells % 100 == 0 {
                info!("Processed {} cells.", num_cells);
            }

            let mut astore = (records_for_barcode, &(*txps), barcode);

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
        let info = get_single_cell_json_info(args, &seqcol_digest);
        write_function::write_single_cell_output(&args.output, info, header, &trimat)?;
        Ok(())
    })
}
