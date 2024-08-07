use crate::alignment_parser;
use crate::em;
use crate::prog_opts::Args;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::write_function;
use crossbeam::queue::ArrayQueue;
use noodles_bam as bam;
use path_tools::WithAdditionalExtension;
use serde_json::json;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::{error, subscriber, Level};

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

pub fn quantify_single_cell_from_collated_bam<R: BufRead>(
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
        let bc_writer = Arc::new(Mutex::new(QuantOutputInfo {
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
        type QueueElement<'a> = (InMemoryAlignmentStore<'a>, Vec<u8>);

        let q: Arc<ArrayQueue<Arc<QueueElement>>> = Arc::new(ArrayQueue::new(4 * nthreads));
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
                            max_iter: args.max_em_iter,
                            convergence_thresh: args.convergence_thresh,
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

            let mut aln_store = InMemoryAlignmentStore::new(filter_opts.clone(), header);
            let _num_rec_processed = alignment_parser::parse_alignments_for_barcode(
                &mut aln_store,
                txps,
                &mut peekable_bam_iter,
                &barcode,
                &mut records_for_read,
            )?;

            let mut astore = std::sync::Arc::new((aln_store, barcode));

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
