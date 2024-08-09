use crate::alignment_parser;
use crate::em;
use crate::kde_utils;
use crate::prog_opts::Args;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, TranscriptInfo,
};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{write_infrep_file, write_output};
use crate::{binomial_continuous_prob, normalize_read_probs};
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};
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
