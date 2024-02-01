use std::fmt;
use std::num::NonZeroUsize;
use serde::Deserialize;

use std::iter::FromIterator;
use tabled::builder::Builder;
use tabled::settings::Style;

use typed_builder::TypedBuilder;

use bio_types::strand::Strand;
use noodles_sam as sam;
use sam::record::data::field::tag;
use tracing::{info, warn};

#[derive(Clone, Debug, PartialEq)]
pub struct AlnInfo {
    pub ref_id: u32,
    pub start: u32,
    pub end: u32,
    pub prob: f64,
    pub strand: Strand,
}

impl AlnInfo {
    pub fn alignment_span(&self) -> u32 {
        self.end - self.start
    }
}

impl From<&sam::alignment::record::Record> for AlnInfo {
    fn from(aln: &sam::alignment::record::Record) -> Self {
        Self {
            ref_id: aln.reference_sequence_id().unwrap() as u32,
            start: aln.alignment_start().unwrap().get() as u32,
            end: aln.alignment_end().unwrap().get() as u32,
            prob: 0.0_f64,
            strand: if aln.flags().is_reverse_complemented() {
                Strand::Reverse
            } else {
                Strand::Forward
            },
        }
    }
}

/// Holds the info relevant for reading the short read quants
#[derive(Debug, Deserialize, Clone)]
pub struct Record {
    pub Name: String,
    pub Length: i32,
    pub EffectiveLength: f64,
    pub TPM: f64,
    pub NumReads: f64,
}

/// Holds the info relevant for running the EM algorithm
pub struct EMInfo<'eqm, 'tinfo> {
    // the read alignment infomation we'll need to
    // perform the EM
    pub eq_map: &'eqm InMemoryAlignmentStore,
    // relevant information about each target transcript
    pub txp_info: &'tinfo mut Vec<TranscriptInfo>,
    // maximum number of iterations the EM will run
    // before returning an estimate.
    pub max_iter: u32,
    // the EM will terminate early if *all* parameters
    // have converged within this threshold of relative
    // change between two subsequent iterations.
    pub convergence_thresh: f64,
}



#[derive(Debug, PartialEq)]
pub struct TranscriptInfo {
    pub len: NonZeroUsize,
    pub total_weight: f64,
    coverage_bins: Vec<f64>,
    pub ranges: Vec<std::ops::Range<u32>>,
    pub coverage_prob: Vec<f32>,
    pub lenf: f64,
}

impl TranscriptInfo {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            len: NonZeroUsize::new(0).unwrap(),
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; 10],
            ranges: Vec::new(),
            coverage_prob: Vec::new(),
            lenf: 0_f64,
        }
    }

    pub fn with_len(len: NonZeroUsize) -> Self {
        Self {
            len,
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; 10],
            ranges: Vec::new(),
            coverage_prob: Vec::new(),
            lenf: len.get() as f64,
        }
    }

    pub fn add_interval(&mut self, start: u32, stop: u32, weight: f64) {
        const NUM_INTERVALS: usize = 10_usize;
        const NUM_INTERVALS_F: f64 = 10.0_f64;
        // find the starting bin
        let start_bin = ((((start as f64) / self.lenf) * NUM_INTERVALS_F).floor() as usize)
            .min(NUM_INTERVALS - 1);
        let end_bin = ((((stop as f64) / self.lenf) * NUM_INTERVALS_F).floor() as usize)
            .min(NUM_INTERVALS - 1);
        for bin in &mut self.coverage_bins[start_bin..=end_bin] {
            *bin += weight;
        }
        self.total_weight += weight;
    }


    pub fn clear_coverage_dist(&mut self) {
        self.coverage_bins.fill(0.0_f64);
        self.total_weight = 0.0_f64;
    }
}

#[derive(Debug)]
pub struct InMemoryAlignmentStore {
    pub filter_opts: AlignmentFilters,
    pub alignments: Vec<AlnInfo>,
    pub as_probabilities: Vec<f32>,
    pub coverage_probabilities: Vec<f64>,
    // holds the boundaries between records for different reads
    boundaries: Vec<usize>,
    pub discard_table: DiscardTable,
}

pub struct InMemoryAlignmentStoreIter<'a> {
    pub store: &'a InMemoryAlignmentStore,
    pub idx: usize,
}

impl<'a> Iterator for InMemoryAlignmentStoreIter<'a> {
    type Item = (&'a [AlnInfo], &'a [f32], &'a [f64]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx + 1 >= self.store.boundaries.len() {
            None
        } else {
            let start = self.store.boundaries[self.idx];
            let end = self.store.boundaries[self.idx + 1];
            self.idx += 1;
            Some((
                &self.store.alignments[start..end],
                &self.store.as_probabilities[start..end],
                &self.store.coverage_probabilities[start..end],
            ))
        }
    }
}

impl InMemoryAlignmentStore {
    pub fn new(fo: AlignmentFilters) -> Self {
        InMemoryAlignmentStore {
            filter_opts: fo,
            alignments: vec![],
            as_probabilities: vec![],
            coverage_probabilities: vec![],
            boundaries: vec![0],
            discard_table: DiscardTable::new(),
        }
    }

    pub fn iter(&self) -> InMemoryAlignmentStoreIter {
        InMemoryAlignmentStoreIter {
            store: self,
            idx: 0,
        }
    }

    pub fn add_group(
        &mut self,
        txps: &mut [TranscriptInfo],
        ag: &mut Vec<sam::alignment::record::Record>,
    ) {
        let (alns, as_probs) = self.filter_opts.filter(&mut self.discard_table, txps, ag);
        if !alns.is_empty() {
            self.alignments.extend_from_slice(&alns);
            self.as_probabilities.extend_from_slice(&as_probs);
            self.coverage_probabilities.extend(vec![0.0_f64; self.alignments.len()]);
            self.boundaries.push(self.alignments.len());
            for a in alns {
                txps[a.ref_id as usize].ranges.push(
                    a.start..a.end,
                    );
            }
        }
    }

    pub fn total_len(&self) -> usize {
        self.alignments.len()
    }

    pub fn num_aligned_reads(&self) -> usize {
        if !self.boundaries.is_empty() {
            self.boundaries.len() - 1
        } else {
            0
        }
    }
}
#[derive(TypedBuilder, Debug)]
pub struct AlignmentFilters {
    five_prime_clip: u32,
    three_prime_clip: i64,
    score_threshold: f32,
    min_aligned_fraction: f32,
    min_aligned_len: u32,
    allow_rc: bool,
    pub model_coverage: bool,
}

#[derive(Debug)]
pub struct DiscardTable {
    discard_5p: u32,
    discard_3p: u32,
    discard_score: u32,
    discard_aln_frac: u32,
    discard_aln_len: u32,
    discard_ori: u32,
    discard_supp: u32,
    valid_best_aln: u32,
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
            valid_best_aln: 0,
        }
    }
}

impl DiscardTable {
    pub fn to_table(&self) -> tabled::tables::Table {
        let d5 = format!("{}", self.discard_5p);
        let d3 = format!("{}", self.discard_3p);
        let dscore = format!("{}", self.discard_score);
        let dfrac = format!("{}", self.discard_aln_frac);
        let dlen = format!("{}", self.discard_aln_len);
        let dori = format!("{}", self.discard_ori);
        let dsupp = format!("{}", self.discard_supp);
        let vread = format!("{}", self.valid_best_aln);

        let data = vec![
            ["reason", "count"],
            ["too far from 5' end", &d5],
            ["too far from 3' end", &d3],
            ["score too low", &dscore],
            ["aligned fraction too low", &dfrac],
            ["aligned length too short", &dlen],
            ["inconsistent orientation", &dori],
            ["supplementary alignment", &dsupp],
            ["reads with valid best alignment", &vread],
        ];
        let mut binding = Builder::from_iter(data).build();
        let table = binding.with(Style::rounded());
        table.clone()
    }
}

impl fmt::Display for DiscardTable {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "discarded because of distance from 5' {}\n",
            self.discard_5p
        )
        .expect("couldn't format discard table.");
        write!(
            f,
            "discarded because of distance from 3' {}\n",
            self.discard_3p
        )
        .expect("couldn't format discard table.");
        write!(
            f,
            "discarded because of score fraction {}\n",
            self.discard_score
        )
        .expect("couldn't format discard table.");
        write!(
            f,
            "discarded because of aligned fraction {}\n",
            self.discard_aln_frac
        )
        .expect("couldn't format discard table.");
        write!(
            f,
            "discarded because of aligned length {}\n",
            self.discard_aln_len
        )
        .expect("couldn't format discard table.");
        write!(
            f,
            "discarded because of aligned orientation {}\n",
            self.discard_ori
        )
        .expect("couldn't format discard table.");
        write!(
            f,
            "discarded because alignment is supplemental {}\n",
            self.discard_supp
        )
    }
}

impl AlignmentFilters {
    fn filter(
        &mut self,
        discard_table: &mut DiscardTable,
        txps: &mut [TranscriptInfo],
        ag: &mut Vec<sam::alignment::record::Record>,
    ) -> (Vec<AlnInfo>, Vec<f32>) {
        // track the best score of any alignment we've seen
        // so far for this read (this will designate the
        // "primary" alignment for the read).
        let mut best_retained_score = i32::MIN;
        // the fraction of the best read that is contained in
        // its alignment.
        let mut aln_frac_at_best_retained = 0_f32;
        // the number of nucleotides contained in the alignment
        // for the best read.
        let mut aln_len_at_best_retained = 0_u32;

        // because of the SAM format, we can't rely on the sequence
        // itself to be present in each alignment record (only the primary)
        // so, here we explicitly look for a non-zero sequence.
        let seq_len = ag
            .iter()
            .find_map(|x| {
                let l = x.sequence().len();
                if l > 0 {
                    Some(l as u32)
                } else {
                    None
                }
            })
            .unwrap_or(0_u32);

        ag.retain(|x| {
            if !x.flags().is_unmapped() {
                let tid = x.reference_sequence_id().unwrap();
                let aln_span = x.alignment_span() as u32;
                let score = x
                    .data()
                    .get(&tag::ALIGNMENT_SCORE)
                    .expect("could not get value")
                    .as_int()
                    .unwrap_or(i32::MIN as i64) as i32;

                // the alignment is to the - strand
                let is_rc = x.flags().is_reverse_complemented();
                if is_rc && !self.allow_rc {
                    discard_table.discard_ori += 1;
                    return false;
                }

                // the alignment is supplementary
                // *NOTE*: this removes "supplementary" alignments, *not*
                // "secondary" alignments.
                let is_supp = x.flags().is_supplementary();
                if is_supp {
                    discard_table.discard_supp += 1;
                    return false;
                }

                // enough absolute sequence (# of bases) is aligned
                let filt_aln_len = (aln_span as u32) < self.min_aligned_len;
                if filt_aln_len {
                    discard_table.discard_aln_len += 1;
                    return false;
                }

                // not too far from the 3' end
                let filt_3p = (x.alignment_end().unwrap().get() as i64)
                    <= (txps[tid].len.get() as i64 - self.three_prime_clip);
                if filt_3p {
                    discard_table.discard_3p += 1;
                    return false;
                }

                // not too far from the 5' end
                let filt_5p = (x.alignment_start().unwrap().get() as u32) >= self.five_prime_clip;
                if filt_5p {
                    discard_table.discard_5p += 1;
                    return false;
                }

                // at this point, we've committed to retaining this
                // alignment. if it has the best score, then record this
                // as the best retained alignment score seen so far
                if score > best_retained_score {
                    best_retained_score = score;
                    // and record the alignment's length
                    // and covered fraction
                    aln_len_at_best_retained = aln_span;
                    aln_frac_at_best_retained = if seq_len > 0 {
                        (aln_span as f32) / (seq_len as f32)
                    } else {
                        0_f32
                    };
                }
                true
            } else {
                false
            }
        });

        if ag.is_empty() || aln_len_at_best_retained == 0 || best_retained_score <= 0 {
            // There were no valid alignments
            return (vec![], vec![]);
        }
        if aln_frac_at_best_retained < self.min_aligned_fraction {
            // The best retained alignment did not have sufficient
            // coverage to be kept
            discard_table.discard_aln_frac += 1;
            return (vec![], vec![]);
        }

        // if we got here, then we have a valid "best" alignment
        discard_table.valid_best_aln += 1;

        let mut probabilities = Vec::<f32>::with_capacity(ag.len());
        let mscore = best_retained_score as f32;
        let inv_max_score = 1.0 / mscore;

        // get a vector of all of the scores
        let mut scores: Vec<i32> = ag
            .iter_mut()
            .map(|a| {
                a.data()
                    .get(&tag::ALIGNMENT_SCORE)
                    .expect("could not get value")
                    .as_int()
                    .unwrap_or(0) as i32
            })
            .collect();

        for (i, score) in scores.iter_mut().enumerate() {
            let fscore = *score as f32;
            let score_ok = (fscore * inv_max_score) >= self.score_threshold; //>= thresh_score;
            if score_ok {
                let f = 10_f32 * ((fscore - mscore) / mscore);
                probabilities.push(f.exp());

                let tid = ag[i].reference_sequence_id().unwrap();
                txps[tid].add_interval(
                    ag[i].alignment_start().unwrap().get() as u32,
                    ag[i].alignment_end().unwrap().get() as u32,
                    1.0_f64,
                );
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

        (ag.iter().map(|x| x.into()).collect(), probabilities)
    }
}


struct FullLengthProbs {
    length_bins: Vec<usize>,
    probs: Vec<Vec<f64>>,
    new_probs: Vec<Vec<f64>>,
}

impl FullLengthProbs {
    const NUM_BINS: usize = 10;

    fn get_bin(span: u32, len: f64) -> usize {
        let len_frac = (len - span as f64) / len;
        (len_frac * FullLengthProbs::NUM_BINS as f64)
            .round()
            .min((FullLengthProbs::NUM_BINS - 1) as f64) as usize
    }

    fn from_data(
        eq_map: &InMemoryAlignmentStore,
        tinfo: &[TranscriptInfo],
        len_bins: &[usize],
    ) -> Self {
        let mut len_probs: Vec<Vec<f64>> =
            vec![vec![0.0f64; FullLengthProbs::NUM_BINS]; len_bins.len()];
        let new_probs: Vec<Vec<f64>> =
            vec![vec![0.0f64; FullLengthProbs::NUM_BINS]; len_bins.len()];

        for (alns, probs, coverage_probs) in eq_map.iter() {
            let inc = 1f64 / probs.len() as f64;
            for a in alns {
                let target_id = a.ref_id as usize;
                let tlenf = tinfo[target_id].lenf;
                let tlen = tlenf as usize;
                let lindex = match len_bins.binary_search(&tlen) {
                    Ok(i) => i,
                    Err(i) => i,
                };
                let bin_num = FullLengthProbs::get_bin(a.alignment_span(), tlenf);
                len_probs[lindex][bin_num] += inc;
            }
        }

        for lb in len_probs.iter_mut() {
            let tot: f64 = (*lb).iter().sum();
            if tot > 0.0 {
                lb.iter_mut().for_each(|v| *v /= tot);
            }
        }

        FullLengthProbs {
            length_bins: len_bins.to_vec(),
            probs: len_probs,
            new_probs,
        }
    }

    fn get_prob_for(&self, a: &AlnInfo, len: usize) -> f64 {
        let lindex = match self.length_bins.binary_search(&len) {
            Ok(i) => i,
            Err(i) => i,
        };
        let lenf = len as f64;
        let bindex = FullLengthProbs::get_bin(a.alignment_span(), lenf);
        self.probs[lindex][bindex]
    }

    fn update_probs(&mut self, a: &AlnInfo, len: usize, inc: f64) {
        let lindex = match self.length_bins.binary_search(&len) {
            Ok(i) => i,
            Err(i) => i,
        };
        let lenf = len as f64;
        let bindex = FullLengthProbs::get_bin(a.alignment_span(), lenf);
        self.new_probs[lindex][bindex] += inc;
    }

    fn swap_probs(&mut self) {
        std::mem::swap(&mut self.probs, &mut self.new_probs);
        for lb in self.probs.iter_mut() {
            let tot: f64 = (*lb).iter().sum();
            if tot > 0.0 {
                lb.iter_mut().for_each(|v| *v /= tot);
            }
        }
        self.new_probs.fill(vec![0.0f64; FullLengthProbs::NUM_BINS]);
    }
}
