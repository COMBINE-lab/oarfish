use clap::Parser;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader},
    num::NonZeroUsize,
};

use tracing::{info, warn};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
use num_format::{Locale, ToFormattedString};

use bio_types::annot::loc::Loc;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::Strand;
use noodles_bam as bam;
use noodles_gtf as gtf;
use noodles_gtf::record::Strand as NoodlesStrand;
use noodles_sam as sam;
use sam::record::data::field::tag;

use bio_types::annot::contig::Contig;
use coitrees::{COITree, IntervalNode};
use nested_intervals::IntervalSet;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser)]
    alignments: String,
}

#[derive(Debug, PartialEq)]
struct TranscriptInfo {
    len: NonZeroUsize,
    ranges: Vec<std::ops::Range<u32>>,
    coverage: f32,
}

impl TranscriptInfo {
    fn new() -> Self {
        Self {
            len: NonZeroUsize::new(0).unwrap(),
            ranges: Vec::new(),
            coverage: 0.0,
        }
    }
    fn with_len(len: NonZeroUsize) -> Self {
        Self {
            len,
            ranges: Vec::new(),
            coverage: 0.0,
        }
    }
}

#[derive(Debug)]
struct InMemoryAlignmentStore {
    alignments: Vec<sam::alignment::record::Record>,
    probabilities: Vec<f32>,
    // holds the boundaries between records for different reads
    boundaries: Vec<usize>,
}

struct InMemoryAlignmentStoreIter<'a> {
    store: &'a InMemoryAlignmentStore,
    idx: usize,
}

impl<'a> Iterator for InMemoryAlignmentStoreIter<'a> {
    type Item = (&'a [sam::alignment::record::Record], &'a [f32]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx + 1 >= self.store.boundaries.len() {
            None
        } else {
            let start = self.store.boundaries[self.idx];
            let end = self.store.boundaries[self.idx + 1];
            self.idx += 1;
            Some((
                &self.store.alignments[start..end],
                &self.store.probabilities[start..end],
            ))
        }
    }
}

impl InMemoryAlignmentStore {
    fn new() -> Self {
        InMemoryAlignmentStore {
            alignments: vec![],
            probabilities: vec![],
            boundaries: vec![0],
        }
    }

    fn iter(&self) -> InMemoryAlignmentStoreIter {
        InMemoryAlignmentStoreIter {
            store: &self,
            idx: 0,
        }
    }

    fn normalize_group_score(&mut self, ag: &mut Vec<sam::alignment::record::Record>) {
        const THRESH: f32 = 0.95;
        if ag.len() == 1 {
            self.probabilities.push(1.0)
        } else {
            let mut max_score = 0_i32;
            let mut scores = Vec::<i32>::with_capacity(ag.len());

            // get the maximum score and fill in the score array
            for a in ag.iter() {
                let score_value = a
                    .data()
                    .get(&tag::ALIGNMENT_SCORE)
                    .expect("could not get value");
                let score = score_value.as_int().unwrap() as i32;
                scores.push(score);
                if score > max_score {
                    max_score = score;
                }
            }

            for score in scores.iter_mut() {
                let fscore = *score as f32;
                let mscore = max_score as f32;
                if fscore > THRESH * mscore {
                    let f = (fscore - mscore) / 10.0_f32;
                    self.probabilities.push(f.exp());
                } else {
                    *score = -1;
                }
            }

            let mut index = 0;
            ag.retain(|_| {
                index += 1;
                scores[index - 1] > -1
            })
        }
    }

    fn add_group(&mut self, ag: &mut Vec<sam::alignment::record::Record>) {
        self.normalize_group_score(ag);
        self.alignments.extend_from_slice(&ag);
        self.boundaries.push(self.alignments.len());
    }

    fn total_len(&self) -> usize {
        self.alignments.len()
    }

    fn num_aligned_reads(&self) -> usize {
        if self.boundaries.len() > 0 {
            self.boundaries.len() - 1
        } else {
            0
        }
    }

    /*
    fn normalize_scores(&mut self) {
        self.probabilities = vec![0.0_f32; self.alignments.len()];
        for w in self.boundaries.windows(2) {
            let s: usize = w[0];
            let e: usize = w[1];
            if e - s > 1 {
                let mut max_score = 0_i32;
                let mut scores = Vec::<i32>::with_capacity(e - s);
                for a in &self.alignments[s..e] {
                    let score_value = a
                        .data()
                        .get(&tag::ALIGNMENT_SCORE)
                        .expect("could not get value");
                    let score = score_value.as_int().unwrap() as i32;
                    scores.push(score);
                    if score > max_score {
                        max_score = score;
                    }
                }
                for (i, score) in scores.iter().enumerate() {
                    let f = ((*score as f32) - (max_score as f32)) / 10.0_f32;
                    self.probabilities[s + i] = f.exp();
                }
            } else {
                self.probabilities[s] = 1.0
            }
        }
    }
    */
}

/// Holds the info relevant for running the EM algorithm
struct EMInfo<'eqm, 'tinfo> {
    eq_map: &'eqm InMemoryAlignmentStore,
    txp_info: &'tinfo Vec<TranscriptInfo>,
    max_iter: u32,
}

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

#[inline]
fn m_step(
    eq_map: &InMemoryAlignmentStore,
    tinfo: &[TranscriptInfo],
    prob_full: &FullLengthProbs,
    prev_count: &mut [f64],
    curr_counts: &mut [f64],
) {
    let mut len_probs: Vec<f64> = Vec::new();
    for (alns, probs) in eq_map.iter() {
        let mut denom = 0.0_f64;
        for (a, p) in alns.iter().zip(probs.iter()) {
            let target_id = a.reference_sequence_id().unwrap();
            let prob = *p as f64;
            /*
            let cov_prob = if tinfo[target_id].coverage < 0.15 { 1e-5 } else { 1.0 };//powf(2.0) as f64;

            let tlen = tinfo[target_id].len.get();
            let is_fl = (tlen - a.alignment_span()) < 20;
            let fl_prob = if is_fl {
                prob_full.get_prob_for(tlen)
            } else {
                1.0 - prob_full.get_prob_for(tlen)
            };
            len_probs.push(fl_prob);
            */
            let cov_prob = 1.0;
            let fl_prob = 1.0;
            len_probs.push(fl_prob);
            denom += prev_count[target_id] * prob * cov_prob * fl_prob;
        }

        if denom > 1e-8 {
            for (i, (a, p)) in alns.iter().zip(probs.iter()).enumerate() {
                let target_id = a.reference_sequence_id().unwrap();
                let prob = *p as f64;
                //let cov_prob = if tinfo[target_id].coverage < 0.15 { 1e-5 } else { 1.0 };//powf(2.0) as f64;
                let cov_prob = 1.0;
                let len_prob = len_probs[i];
                let inc = (prev_count[target_id] * prob * cov_prob * len_prob) / denom;
                curr_counts[target_id] += inc;
            }
        }
        len_probs.clear();
    }

    //eprintln!("nunique : {}, ntotal : {}", nunique, eq_map.num_aligned_reads());
    //prob_full_len / (eq_map.num_aligned_reads() as f64)
}

fn em(em_info: &EMInfo) -> Vec<f64> {
    let eq_map = em_info.eq_map;
    let tinfo = em_info.txp_info;
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    // init
    let avg = total_weight / (tinfo.len() as f64);
    let mut prev_counts = vec![avg; tinfo.len()];
    let mut curr_counts = vec![0.0f64; tinfo.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;
    let mut fl_prob = 0.5f64;

    let length_bins = vec![
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
    let len_probs = FullLengthProbs::from_data(eq_map, tinfo, &length_bins);

    while niter < max_iter {
        m_step(
            eq_map,
            tinfo,
            &len_probs,
            &mut prev_counts,
            &mut curr_counts,
        );

        //std::mem::swap(&)
        for i in 0..curr_counts.len() {
            if prev_counts[i] > 1e-8 {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if (rel_diff < 1e-3) && (niter > 50) {
            break;
        }
        niter += 1;
        if niter % 10 == 0 {
            info!("iteration {}; rel diff {}", 
                niter.to_formatted_string(&Locale::en), rel_diff);
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| {
        if *x < 1e-8 {
            *x = 0.0
        }
    });
    m_step(
        eq_map,
        tinfo,
        &len_probs,
        &mut prev_counts,
        &mut curr_counts,
    );
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
    let mut num_mapped = 0_u64;
    let mut records_for_read = vec![];
    let mut store = InMemoryAlignmentStore::new();

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
                if let (Some(ref_id), false) = (
                    record.reference_sequence_id(),
                    record.flags().is_supplementary(),
                ) {
                    records_for_read.push(record_copy);
                    txps[ref_id].ranges.push(
                        (record.alignment_start().unwrap().get() as u32)
                            ..(record.alignment_end().unwrap().get() as u32),
                    );
                }
            } else {
                if !prev_read.is_empty() {
                    //println!("the previous read had {} mappings", records_for_read.len());
                    store.add_group(&mut records_for_read);
                    records_for_read.clear();
                    num_mapped += 1;
                }
                prev_read = rstring;
                if let (Some(ref_id), false) = (
                    record.reference_sequence_id(),
                    record.flags().is_supplementary(),
                ) {
                    records_for_read.push(record_copy);
                    txps[ref_id].ranges.push(
                        (record.alignment_start().unwrap().get() as u32)
                            ..(record.alignment_end().unwrap().get() as u32),
                    );
                }
            }
        }
    }
    if !records_for_read.is_empty() {
        store.add_group(&mut records_for_read);
        records_for_read.clear();
        num_mapped += 1;
    }

    info!("computing coverages");
    for t in txps.iter_mut() {
        let interval_set = IntervalSet::new(&t.ranges).expect("couldn't build interval set");
        let mut interval_set = interval_set.merge_connected();
        let len = t.len.get() as u32;
        let covered = interval_set.covered_units();
        t.coverage = (covered as f32) / (len as f32);
    }
    info!("done");

    //eprintln!("Number of mapped reads : {}", num_mapped);
    //eprintln!("normalizing alignment scores");
    //store.normalize_scores();
    info!("Total number of alignment records : {}", store.total_len().to_formatted_string(&Locale::en));
    info!("number of aligned reads : {}", store.num_aligned_reads().to_formatted_string(&Locale::en));

    let emi = EMInfo {
        eq_map: &store,
        txp_info: &txps,
        max_iter: 1000,
    };

    let counts = em(&emi);

    println!("tname\tcoverage\tlen\tnum_reads");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (i, (_rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        println!(
            "{}\t{}\t{}\t{}",
            _rseq,
            txps[i].coverage,
            rmap.length(),
            counts[i]
        );
    }

    Ok(())
}

//
// ignore anything below this line for now
//

#[allow(unused)]
fn main_old() -> io::Result<()> {
    let args = Args::parse();

    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(gtf::Reader::new)?;
    let mut evec = Vec::new();
    let mut tvec = Vec::new();
    let mut tmap = HashMap::new();

    for result in reader.records() {
        let record = result?;
        match record.ty() {
            "exon" => {
                let s: isize = (usize::from(record.start()) as isize) - 1;
                let e: isize = usize::from(record.end()) as isize;
                let l: usize = (e - s).try_into().unwrap();
                let mut t = String::new();
                for e in record.attributes().iter() {
                    if e.key() == "transcript_id" {
                        t = e.value().to_owned();
                    }
                }

                let ni = tmap.len();
                let tid = *tmap.entry(t.clone()).or_insert(ni);

                // if this is what we just inserted
                if ni == tid {
                    tvec.push(Spliced::new(0, 1, 1, Strand::Forward));
                }

                let strand = match record.strand().unwrap() {
                    NoodlesStrand::Forward => Strand::Forward,
                    NoodlesStrand::Reverse => Strand::Reverse,
                };
                let c = Contig::new(tid, s, l, strand);
                evec.push(c);
            }
            "transcript" => {
                let mut t = String::new();
                for e in record.attributes().iter() {
                    if e.key() == "transcript_id" {
                        t = e.value().to_owned();
                    }
                }
                let ni = tmap.len();
                let tid = *tmap.entry(t.clone()).or_insert(ni);

                // if this is what we just inserted
                if ni == tid {
                    tvec.push(Spliced::new(0, 1, 1, Strand::Forward));
                }
            }
            _ => {}
        }
    }

    let mut txp_to_exon = HashMap::new();

    let mut l = 0;
    let mut max_len = 0;
    for (i, e) in evec.iter().enumerate() {
        let mut v = txp_to_exon.entry(e.refid()).or_insert(vec![]);
        v.push(i);
        l = v.len();
        if l > max_len {
            max_len = l;
        }
    }

    let mut txp_features: HashMap<usize, _> = HashMap::new();

    for (k, v) in txp_to_exon.iter_mut() {
        let strand = evec[v[0]].strand();
        v.sort_unstable_by_key(|x| evec[*x as usize].start());
        let s = evec[v[0]].start();
        let starts: Vec<usize> = v
            .iter()
            .map(|e| (evec[*e as usize].start() - s) as usize)
            .collect();
        let lens: Vec<usize> = v.iter().map(|e| evec[*e as usize].length()).collect();
        println!("lens = {:?}, starts = {:?}", lens, starts);
        txp_features.insert(
            **k,
            Spliced::with_lengths_starts(k, s, &lens, &starts, strand).unwrap(),
        );
    }

    let interval_vec: Vec<IntervalNode<usize, usize>> = evec
        .iter()
        .enumerate()
        .map(|(i, e)| {
            IntervalNode::new(
                e.start() as i32,
                (e.start() + e.length() as isize) as i32,
                i,
            )
        })
        .collect();
    let ct = COITree::new(interval_vec);

    println!("parsed {} exons", evec.len());
    println!("parsed {} transcripts", tvec.len());
    println!("max exon transcript had {} exons", max_len);
    Ok(())
}
