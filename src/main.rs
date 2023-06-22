use clap::Parser;

use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{self, BufReader, BufWriter, Write},
};

use num_format::{Locale, ToFormattedString};
use tracing::info;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

use bio_types::annot::loc::Loc;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::Strand;
use noodles_bam as bam;
use noodles_gtf as gtf;
use noodles_gtf::record::Strand as NoodlesStrand;

use bio_types::annot::contig::Contig;
use coitrees::{COITree, IntervalNode};
mod util;
use crate::util::oarfish_types::{AlignmentFilters, InMemoryAlignmentStore, TranscriptInfo};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser, required = true)]
    alignments: String,
    /// Location where output quantification file should be written
    #[clap(short, long, value_parser, required = true)]
    output: String,
    // Maximum allowable distance of the right-most end of an alignment from the 3' transcript end
    #[clap(short, long, value_parser, default_value_t = u32::MAX as i64)]
    three_prime_clip: i64,
    // Maximum allowable distance of the left-most end of an alignment from the 5' transcript end
    #[clap(short, long, value_parser, default_value_t = u32::MAX)]
    five_prime_clip: u32,
    // Fraction of the best possible alignment score that a secondary alignment must have for
    // consideration
    #[clap(short, long, value_parser, default_value_t = 0.95)]
    score_threshold: f32,
    // Fraction of a query that must be mapped within an alignemnt to consider the alignemnt
    // valid
    #[clap(short, long, value_parser, default_value_t = 0.5)]
    min_aligned_fraction: f32,
    // Minimum number of nucleotides in the aligned portion of a read
    #[clap(short = 'l', long, value_parser, default_value_t = 50)]
    min_aligned_len: u32,
    // Allow both forward-strand and reverse-complement alignments
    #[clap(short = 'n', long, value_parser)]
    allow_negative_strand: bool,
    // Apply the coverage model
    #[clap(long, value_parser)]
    model_coverage: bool,
}

/// Holds the info relevant for running the EM algorithm
struct EMInfo<'eqm, 'tinfo> {
    eq_map: &'eqm InMemoryAlignmentStore,
    txp_info: &'tinfo mut Vec<TranscriptInfo>,
    max_iter: u32,
}

/*
struct FullLengthProbs {
    length_bins: Vec<usize>,
    probs: Vec<f64>,
}

impl FullLengthProbs {
    fn from_data(
        eq_map: &InMemoryAlignmentStore,
        tinfo: &[TranscriptInfo],
        len_bins: &[usize],
    ) -> Self {
        let mut num = vec![0; len_bins.len()];
        let mut denom = vec![0; len_bins.len()];

        for (alns, probs) in eq_map.iter() {
            let a = alns.first().unwrap();
            if alns.len() == 1 {
                let target_id = a.reference_sequence_id().unwrap();
                let tlen = tinfo[target_id].len.get();
                let lindex = match len_bins.binary_search(&tlen) {
                    Ok(i) => i,
                    Err(i) => i,
                };
                let is_fl = (tlen - a.alignment_span()) < 20;
                if is_fl {
                    num[lindex] += 1;
                }
                denom[lindex] += 1;
            }
        }

        let probs: Vec<f64> = num
            .iter()
            .zip(denom.iter())
            .map(|(n, d)| {
                if *d > 0 {
                    (*n as f64) / (*d as f64)
                } else {
                    1e-5
                }
            })
            .collect();

        FullLengthProbs {
            length_bins: len_bins.to_vec(),
            probs,
        }
    }

    fn get_prob_for(&self, len: usize) -> f64 {
        let lindex = match self.length_bins.binary_search(&len) {
            Ok(i) => i,
            Err(i) => i,
        };
        self.probs[lindex]
    }
}
*/

#[inline]
fn m_step(
    eq_map: &InMemoryAlignmentStore,
    tinfo: &mut [TranscriptInfo],
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) {
    for (alns, probs) in eq_map.iter() {
        let mut denom = 0.0_f64;
        for (a, p) in alns.iter().zip(probs.iter()) {
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = tinfo[target_id].coverage_prob;
            denom += prev_count[target_id] * prob * cov_prob;
        }

        if denom > 1e-8 {
            for (_i, (a, p)) in alns.iter().zip(probs.iter()).enumerate() {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = tinfo[target_id].coverage_prob;

                let inc = (prev_count[target_id] * prob * cov_prob) / denom;
                curr_counts[target_id] += inc;

                let start = a.start;
                let stop = a.end;
                tinfo[target_id].add_interval(start, stop, inc);
            }
        }
    }
}

fn em(em_info: &mut EMInfo) -> Vec<f64> {
    let eq_map = em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &mut Vec<TranscriptInfo> = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    // init
    let avg = total_weight / (tinfo.len() as f64);
    let mut prev_counts = vec![avg; tinfo.len()];
    let mut curr_counts = vec![0.0f64; tinfo.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;
    let mut _fl_prob = 0.5f64;

    let _length_bins = vec![
        200,
        500,
        1000,
        2000,
        5000,
        10000,
        15000,
        20000,
        25000,
        usize::MAX,
    ];
    //let len_probs = FullLengthProbs::from_data(eq_map, tinfo, &length_bins);

    while niter < max_iter {
        m_step(eq_map, tinfo, &mut prev_counts, &mut curr_counts);

        //std::mem::swap(&)
        for i in 0..curr_counts.len() {
            if prev_counts[i] > 1e-8 {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
            if fops.model_coverage {
                tinfo[i].compute_coverage_prob();
                tinfo[i].clear_coverage_dist();
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if (rel_diff < 1e-3) && (niter > 10) {
            break;
        }
        niter += 1;
        if niter % 10 == 0 {
            info!(
                "iteration {}; rel diff {}",
                niter.to_formatted_string(&Locale::en),
                rel_diff
            );
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| {
        if *x < 1e-8 {
            *x = 0.0
        }
    });
    m_step(eq_map, tinfo, &mut prev_counts, &mut curr_counts);
    curr_counts
}

fn main() -> io::Result<()> {
    tracing_subscriber::registry()
        .with(fmt::layer().with_writer(io::stderr))
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let args = Args::parse();

    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(bam::Reader::new)?;

    let filter_opts = AlignmentFilters::builder()
        .five_prime_clip(args.five_prime_clip)
        .three_prime_clip(args.three_prime_clip)
        .score_threshold(args.score_threshold)
        .min_aligned_fraction(args.min_aligned_fraction)
        .min_aligned_len(args.min_aligned_len)
        .allow_rc(args.allow_negative_strand)
        .model_coverage(args.model_coverage)
        .build();

    let header = reader.read_header()?;

    for (prog, _pmap) in header.programs().iter() {
        info!("program: {}", prog);
    }

    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(header.reference_sequences().len());

    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (_rseq, rmap) in header.reference_sequences().iter() {
        // println!("ref: {}, rmap : {:?}", rseq, rmap.length());
        txps.push(TranscriptInfo::with_len(rmap.length()));
    }

    //let mut rmap = HashMap<usize, ::new();
    //
    let mut prev_read = String::new();
    let mut _num_mapped = 0_u64;
    let mut records_for_read = vec![];
    let mut store = InMemoryAlignmentStore::new(filter_opts);

    for result in reader.records(&header) {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        let record_copy = record.clone();
        if let Some(rname) = record.read_name() {
            let rstring: String =
                <noodles_sam::record::read_name::ReadName as AsRef<str>>::as_ref(rname).to_owned();
            // if this is an alignment for the same read, then
            // push it onto our temporary vector.
            if prev_read == rstring {
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            } else {
                if !prev_read.is_empty() {
                    //println!("the previous read had {} mappings", records_for_read.len());
                    store.add_group(&mut txps, &mut records_for_read);
                    records_for_read.clear();
                    _num_mapped += 1;
                }
                prev_read = rstring;
                if let Some(_ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                }
            }
        }
    }
    if !records_for_read.is_empty() {
        store.add_group(&mut txps, &mut records_for_read);
        records_for_read.clear();
        _num_mapped += 1;
    }

    info!("discard_table: {:?}", store.discard_table);

    if store.filter_opts.model_coverage {
        info!("computing coverages");
        for t in txps.iter_mut() {
            t.compute_coverage_prob();
        }
    }

    info!("done");

    info!(
        "Total number of alignment records : {}",
        store.total_len().to_formatted_string(&Locale::en)
    );
    info!(
        "number of aligned reads : {}",
        store.num_aligned_reads().to_formatted_string(&Locale::en)
    );

    let mut emi = EMInfo {
        eq_map: &store,
        txp_info: &mut txps,
        max_iter: 1000,
    };

    let counts = em(&mut emi);

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(args.output)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "tname\tcoverage\tlen\tnum_reads").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (i, (_rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            _rseq,
            txps[i].coverage_prob,
            rmap.length(),
            counts[i]
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}
