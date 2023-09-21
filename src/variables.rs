use std::num::NonZeroUsize;
use noodles_sam as sam;
use sam::record::data::field::tag;
use serde::Deserialize;

#[derive(Debug, PartialEq)]
pub struct TranscriptInfo {
    pub len: NonZeroUsize,
    pub ranges: Vec<std::ops::Range<u32>>,
    pub coverage: f64,
    pub coverage_prob: Vec<f32>,
}

impl TranscriptInfo {
    pub fn new() -> Self {
        Self {
            len: NonZeroUsize::new(0).unwrap(),
            ranges: Vec::new(),
            coverage: 0.0,
            coverage_prob: Vec::new(),
        }
    }
    pub fn with_len(len: NonZeroUsize) -> Self {
        Self {
            len,
            ranges: Vec::new(),
            coverage: 0.0,
            coverage_prob: Vec::new(),
        }
    }
}

#[derive(Debug)]
pub struct InMemoryAlignmentStore {
    pub alignments: Vec<sam::alignment::record::Record>,
    pub probabilities: Vec<f32>,
    // holds the boundaries between records for different reads
    pub boundaries: Vec<usize>,
}

pub struct InMemoryAlignmentStoreIter<'a> {
    pub store: &'a InMemoryAlignmentStore,
    pub idx: usize,
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
    pub fn new() -> Self {
        InMemoryAlignmentStore {
            alignments: vec![],
            probabilities: vec![],
            boundaries: vec![0],
        }
    }

    pub fn iter(&self) -> InMemoryAlignmentStoreIter {
        InMemoryAlignmentStoreIter {
            store: &self,
            idx: 0,
        }
    }

    pub fn add_group(&mut self, ag: &Vec<sam::alignment::record::Record>) {
        self.alignments.extend_from_slice(&ag);
        self.boundaries.push(self.alignments.len());
    }

    pub fn total_len(&self) -> usize {
        self.alignments.len()
    }

    pub fn num_aligned_reads(&self) -> usize {
        if self.boundaries.len() > 0 {
            self.boundaries.len() - 1
        } else {
            0
        }
    }

    pub fn normalize_scores(&mut self) {
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
}

/// Holds the info relevant for running the EM algorithm
pub struct EMInfo<'eqm, 'tinfo> {
    pub eq_map: &'eqm InMemoryAlignmentStore,
    pub txp_info: &'tinfo Vec<TranscriptInfo>,
    pub max_iter: u32,
}

#[derive(Debug, Deserialize, Clone)]
pub struct Record {
    pub Name: String,
    pub Length: i32,
    pub EffectiveLength: f64,
    pub TPM: f64,
    pub NumReads: f64,
}