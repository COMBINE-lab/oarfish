use std::num::NonZeroUsize;
use std::fmt;

use typed_builder::TypedBuilder;

use bio_types::strand::Strand;
use noodles_sam as sam;
use sam::record::data::field::tag;

#[derive(Clone, Debug, PartialEq)]
pub struct AlnInfo {
    pub ref_id: u32,
    pub start: u32,
    pub end: u32,
    pub prob: f64,
    pub strand: Strand,
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

#[derive(Debug, PartialEq)]
pub struct TranscriptInfo {
    len: NonZeroUsize,
    total_weight: f64,
    coverage_bins: Vec<f64>,
    pub coverage_prob: f64,
    lenf: f64,
}

impl TranscriptInfo {
    #[allow(dead_code)]
    fn new() -> Self {
        Self {
            len: NonZeroUsize::new(0).unwrap(),
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; 10],
            //ranges: Vec::new(),
            coverage_prob: 1.0,
            lenf: 0_f64,
        }
    }

    pub fn with_len(len: NonZeroUsize) -> Self {
        Self {
            len,
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; 10],
            //ranges: Vec::new(),
            coverage_prob: 1.0,
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

    pub fn compute_coverage_prob(&mut self) {
        // first normalize
        const EVEN_COV: f64 = 0.1_f64;
        let sum: f64 = self.coverage_bins.iter().sum();
        let deviation: f64 = if sum > 0.0 {
            self.coverage_bins
                .iter()
                .copied()
                .fold(0.0_f64, |a, e| a + ((e / sum) - EVEN_COV).abs())
        } else {
            2.0_f64
        };
        let coverage_prob = (-deviation * 4.0).exp();
        self.coverage_prob = coverage_prob;
    }

    pub fn clear_coverage_dist(&mut self) {
        self.coverage_bins.fill(0.0_f64);
        self.total_weight = 0.0_f64;
    }
}

#[derive(Debug)]
pub struct InMemoryAlignmentStore {
    pub filter_opts: AlignmentFilters,
    alignments: Vec<AlnInfo>,
    probabilities: Vec<f32>,
    // holds the boundaries between records for different reads
    boundaries: Vec<usize>,
    pub discard_table: DiscardTable,
}

pub struct InMemoryAlignmentStoreIter<'a> {
    store: &'a InMemoryAlignmentStore,
    idx: usize,
}

impl<'a> Iterator for InMemoryAlignmentStoreIter<'a> {
    type Item = (&'a [AlnInfo], &'a [f32]);

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
    pub fn new(fo: AlignmentFilters) -> Self {
        InMemoryAlignmentStore {
            filter_opts: fo,
            alignments: vec![],
            probabilities: vec![],
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
        let (alns, probs) = self.filter_opts.filter(&mut self.discard_table, txps, ag);
        if !ag.is_empty() {
            self.alignments.extend_from_slice(&alns);
            self.probabilities.extend_from_slice(&probs);
            self.boundaries.push(self.alignments.len());
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

impl fmt::Display for DiscardTable {
    
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "discarded because of distance from 5' {}\n", self.discard_5p).expect("couldn't format discard table.");
        write!(f, "discarded because of distance from 3' {}\n", self.discard_3p).expect("couldn't format discard table.");
        write!(f, "discarded because of score fraction {}\n", self.discard_score).expect("couldn't format discard table.");
        write!(f, "discarded because of aligned fraction {}\n", self.discard_aln_frac).expect("couldn't format discard table.");
        write!(f, "discarded because of aligned length {}\n", self.discard_aln_len).expect("couldn't format discard table.");
        write!(f, "discarded because of aligned orientation {}\n", self.discard_ori).expect("couldn't format discard table.");
        write!(f, "discarded because alignment is supplemental {}\n", self.discard_supp)
    }
}

impl AlignmentFilters {
    fn filter(
        &mut self,
        discard_table: &mut DiscardTable,
        txps: &mut [TranscriptInfo],
        ag: &mut Vec<sam::alignment::record::Record>,
    ) -> (Vec<AlnInfo>, Vec<f32>) {
        ag.retain(|x| {
            if !x.flags().is_unmapped() {
                let tid = x.reference_sequence_id().unwrap();
                let aln_span = x.alignment_span();

                // the read is not aligned to the - strand
                let is_rc = x.flags().is_reverse_complemented();
                if is_rc && !self.allow_rc {
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
                let filt_3p = (x.alignment_end().unwrap().get() as i64)
                    > (txps[tid].len.get() as i64 - self.three_prime_clip);
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
            let mut av = vec![];
            if let Some(a) = ag.first() {
                let tid = a.reference_sequence_id().unwrap();
                txps[tid].add_interval(
                    a.alignment_start().unwrap().get() as u32,
                    a.alignment_end().unwrap().get() as u32,
                    1.0_f64,
                );
                av.push(a.into());
            }
            (av, vec![1.0_f32])
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

            for (i, score) in scores.iter_mut().enumerate() {
                let fscore = *score as f32;
                if fscore >= thresh_score {
                    let f = (fscore - mscore) / 10.0_f32;
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

            /*if ag.is_empty() {
                warn!("No valid scores for read!");
                warn!("max_score = {}. scores = {:?}", max_score, scores);
            }*/
            (ag.iter().map(|x| x.into()).collect(), probabilities)
        }
    }
}
