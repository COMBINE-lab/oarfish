use clap::Parser;

use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{self, BufReader, BufWriter, Write},
    num::NonZeroUsize,
};

use num_format::{Locale, ToFormattedString};
use tracing::{info, warn};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
use typed_builder::TypedBuilder;

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
    #[clap(short, long, value_parser, required = true)]
    alignments: String,
    /// Location where output quantification file should be written
    #[clap(short, long, value_parser, required = true)]
    output: String,
    // Maximum allowable distance of the right-most end of an alignment from the 3' transcript end
    #[clap(short, long, value_parser, default_value_t = u32::MAX)]
    three_prime_clip: u32,
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
    filter_opts: AlignmentFilters,
    alignments: Vec<sam::alignment::record::Record>,
    probabilities: Vec<f32>,
    // holds the boundaries between records for different reads
    boundaries: Vec<usize>,
    pub discard_table: DiscardTable,
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
    fn new(fo: AlignmentFilters) -> Self {
        InMemoryAlignmentStore {
            filter_opts: fo,
            alignments: vec![],
            probabilities: vec![],
            boundaries: vec![0],
            discard_table: DiscardTable::new(),
        }
    }

    fn iter(&self) -> InMemoryAlignmentStoreIter {
        InMemoryAlignmentStoreIter {
            store: self,
            idx: 0,
        }
    }

    fn add_group(
        &mut self,
        txps: &Vec<TranscriptInfo>,
        ag: &mut Vec<sam::alignment::record::Record>,
    ) {
        let probs = self.filter_opts.filter(&mut self.discard_table, txps, ag);
        if !ag.is_empty() {
            self.alignments.extend_from_slice(ag);
            self.probabilities.extend_from_slice(&probs);
            self.boundaries.push(self.alignments.len());
        }
    }

    fn total_len(&self) -> usize {
        self.alignments.len()
    }

    fn num_aligned_reads(&self) -> usize {
        if !self.boundaries.is_empty() {
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
    m_step(
        eq_map,
        tinfo,
        &len_probs,
        &mut prev_counts,
        &mut curr_counts,
    );
    curr_counts
}

#[derive(TypedBuilder, Debug)]
struct AlignmentFilters {
    five_prime_clip: u32,
    three_prime_clip: u32,
    score_threshold: f32,
    min_aligned_fraction: f32,
    min_aligned_len: u32,
}

#[derive(Debug)]
struct DiscardTable {
    discard_5p: u32,
    discard_3p: u32,
    discard_score: u32,
    discard_aln_frac: u32,
    discard_aln_len: u32,
    discard_ori: u32,
    discard_supp: u32,
}

impl DiscardTable {
    fn new() -> Self {
        DiscardTable {
            discard_5p: 0,
            discard_3p: 0,
            discard_score: 0,
            discard_aln_frac: 0,
            discard_aln_len: 0,
            discard_ori: 0,
            discard_supp: 0,
        }
    }
}

impl AlignmentFilters {
    fn filter(
        &mut self,
        discard_table: &mut DiscardTable,
        txps: &Vec<TranscriptInfo>,
        ag: &mut Vec<sam::alignment::record::Record>,
    ) -> Vec<f32> {
        ag.retain(|x| {
            if !x.flags().is_unmapped() {
                let tid = x.reference_sequence_id().unwrap();
                let aln_span = x.alignment_span();

                // the read is not aligned to the - strand
                let filt_ori = !x.flags().is_reverse_complemented();
                if !filt_ori {
                    discard_table.discard_ori += 1;
                    return false;
                }

                let filt_supp = !x.flags().is_supplementary();
                if !filt_supp {
                    discard_table.discard_supp += 1;
                    return false;
                }

                // enough absolute sequence (# of bases) is aligned
                let filt_aln_len = (aln_span as u32) > self.min_aligned_len;
                if !filt_aln_len {
                    discard_table.discard_aln_len += 1;
                    return false;
                }

                // not too far from the 5' end
                let filt_5p = (x.alignment_start().unwrap().get() as u32) < self.five_prime_clip;
                if !filt_5p {
                    discard_table.discard_5p += 1;
                    return false;
                }

                // not too far from the 3' end
                let filt_3p = (x.alignment_end().unwrap().get() as u32)
                    > (txps[tid].len.get() as u32 - self.three_prime_clip);
                if !filt_3p {
                    discard_table.discard_3p += 1;
                    return false;
                }

                // enough of the read is aligned
                let filt_aln_frac =
                    ((aln_span as f32) / (x.sequence().len() as f32)) >= self.min_aligned_fraction;
                if !filt_aln_frac {
                    discard_table.discard_aln_frac += 1;
                    return false;
                }

                true
            } else {
                false
            }
        });

        if ag.len() == 1 {
            vec![1.0_f32]
        } else {
            let mut max_score = i32::MIN;
            let mut scores = Vec::<i32>::with_capacity(ag.len());
            let mut probabilities = Vec::<f32>::with_capacity(ag.len());

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

            let mscore = max_score as f32;
            let thresh_score = mscore - (mscore.abs() * (1.0 - self.score_threshold));

            for score in scores.iter_mut() {
                let fscore = *score as f32;
                if fscore >= thresh_score {
                    let f = (fscore - mscore) / 10.0_f32;
                    probabilities.push(f.exp());
                } else {
                    *score = i32::MIN;
                    discard_table.discard_score += 1;
                }
            }

            let mut index = 0;
            ag.retain(|_| {
                index += 1;
                scores[index - 1] > i32::MIN
            });

            /*if ag.is_empty() {
                warn!("No valid scores for read!");
                warn!("max_score = {}. scores = {:?}", max_score, scores);
            }*/
            probabilities
        }
    }
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
    let mut num_mapped = 0_u64;
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
                if let Some(ref_id) = record.reference_sequence_id() {
                    records_for_read.push(record_copy);
                    txps[ref_id].ranges.push(
                        (record.alignment_start().unwrap().get() as u32)
                            ..(record.alignment_end().unwrap().get() as u32),
                    );
                }
            } else {
                if !prev_read.is_empty() {
                    //println!("the previous read had {} mappings", records_for_read.len());
                    store.add_group(&txps, &mut records_for_read);
                    records_for_read.clear();
                    num_mapped += 1;
                }
                prev_read = rstring;
                if let Some(ref_id) = record.reference_sequence_id() {
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
        store.add_group(&txps, &mut records_for_read);
        records_for_read.clear();
        num_mapped += 1;
    }

    info!("discard_table: {:?}", store.discard_table);

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
    info!(
        "Total number of alignment records : {}",
        store.total_len().to_formatted_string(&Locale::en)
    );
    info!(
        "number of aligned reads : {}",
        store.num_aligned_reads().to_formatted_string(&Locale::en)
    );

    let emi = EMInfo {
        eq_map: &store,
        txp_info: &txps,
        max_iter: 1000,
    };

    let counts = em(&emi);

    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .open(args.output)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    write!(writer, "tname\tcoverage\tlen\tnum_reads\n").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
    for (i, (_rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        write!(
            writer,
            "{}\t{}\t{}\t{}\n",
            _rseq,
            txps[i].coverage,
            rmap.length(),
            counts[i]
        )
        .expect("Couldn't write to output file.");
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
        v.sort_unstable_by_key(|x| evec[*x].start());
        let s = evec[v[0]].start();
        let starts: Vec<usize> = v.iter().map(|e| (evec[*e].start() - s) as usize).collect();
        let lens: Vec<usize> = v.iter().map(|e| evec[*e].length()).collect();
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
