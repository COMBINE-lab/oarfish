use serde::Deserialize;
use std::fmt;
use std::num::NonZeroUsize;

use kders::kde::KDEModel;
use std::iter::FromIterator;
use tabled::builder::Builder;
use tabled::settings::Style;

use seqcol_rs::DigestResult;
use serde::Serialize;
use serde_json::json;
use typed_builder::TypedBuilder;

use bio_types::strand::Strand;
use bstr::{B, ByteSlice};
use noodles_sam as sam;
use sam::{Header, alignment::record::data::field::tag::Tag as AlnTag};

#[allow(unused_imports)]
use tracing::{error, info, warn};

use crate::prog_opts::{ProjProbSource, ReadAssignmentProbOut};
use crate::util::constants::EMPTY_READ_NAME;

pub(crate) struct NamedDigestVec(Vec<(String, DigestResult)>);

impl From<Vec<(String, DigestResult)>> for NamedDigestVec {
    fn from(o: Vec<(String, DigestResult)>) -> Self {
        Self(o)
    }
}

impl NamedDigestVec {
    pub fn new() -> Self {
        NamedDigestVec(Vec::new())
    }

    pub fn push(&mut self, x: (String, DigestResult)) {
        self.0.push(x);
    }

    pub fn to_json(&self) -> serde_json::Value {
        let mut op = serde_json::Map::new();
        op.extend(self.0.iter().map(|(k, v)| (k.clone(), v.to_json())));
        json!(op)
    }
}

// how we can get our raw input
pub(crate) enum InputSourceType {
    Ubam,
    Fastx,
    Unknown,
}

// need both FASTX and UBAM to be able to act as an
// input source
pub(crate) trait ReadSource {
    fn add_to_read_group(&self, rg: &mut ReadChunkWithNames);
}

impl ReadSource for needletail::parser::SequenceRecord<'_> {
    #[inline(always)]
    fn add_to_read_group(&self, rg: &mut ReadChunkWithNames) {
        let read_str = B(self.id())
            .fields_with(|ch| ch.is_ascii_whitespace())
            .next();
        let read_name = if let Some(first_part) = read_str {
            first_part
        } else {
            EMPTY_READ_NAME.as_bytes()
        };

        // put this read on the current chunk
        rg.add_id_and_read(read_name, &self.seq());
    }
}

impl ReadSource for noodles_sam::alignment::RecordBuf {
    #[inline(always)]
    fn add_to_read_group(&self, rg: &mut ReadChunkWithNames) {
        let read_name = if let Some(name) = self.name() {
            name.fields_with(|ch| ch.is_ascii_whitespace())
                .next()
                .unwrap_or(EMPTY_READ_NAME.as_bytes())
        } else {
            EMPTY_READ_NAME.as_bytes()
        };

        // put this read on the current chunk
        rg.add_id_and_read(read_name, self.sequence().as_ref());
    }
}

#[derive(Clone)]
pub(crate) struct ReadChunkWithNames {
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

    pub fn iter(&self) -> ReadChunkIter<'_> {
        ReadChunkIter {
            chunk: self,
            pos: 0,
        }
    }
}

pub struct ReadChunkIter<'a> {
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
        // there's the - 1 because the separator always has
        // one more entry than the actual number of sequences
        // because it starts with a 0.
        let rem = (self.chunk.seq_sep.len() - 1) - self.pos;
        (rem, Some(rem))
    }
}

impl ExactSizeIterator for ReadChunkIter<'_> {}

pub trait AlnRecordLike {
    fn opt_sequence_len(&self) -> Option<usize>;
    fn is_reverse_complemented(&self) -> bool;
    fn is_unmapped(&self) -> bool;
    fn ref_id(&self, header: &Header) -> anyhow::Result<usize>;
    fn aln_span(&self) -> Option<usize>;
    fn aln_score(&self) -> Option<i64>;
    fn aln_start(&self) -> u32;
    fn aln_end(&self) -> u32;
    fn is_supp(&self) -> bool;
    #[allow(dead_code)]
    fn name(&self) -> Option<String>;
}

#[derive(Debug, Clone, PartialEq)]
pub enum CigarOp {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SequenceMatch,
    SequenceMismatch,
}

impl From<u8> for CigarOp {
    fn from(e: u8) -> Self {
        match e {
            0 => CigarOp::Match,
            1 => CigarOp::Insertion,
            2 => CigarOp::Deletion,
            3 => CigarOp::Skip,
            4 => CigarOp::SoftClip,
            5 => CigarOp::HardClip,
            6 => CigarOp::Pad,
            7 => CigarOp::SequenceMatch,
            8 => CigarOp::SequenceMismatch,
            x => {
                error!("invalid cigar code {}", x);
                CigarOp::Skip
            }
        }
    }
}

/// from noodles: https://docs.rs/noodles-sam/latest/src/noodles_sam/alignment/record/cigar/op/kind.rs.html
impl CigarOp {
    #[allow(dead_code)]
    pub fn consumes_read(&self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Insertion
                | Self::SoftClip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }

    #[allow(dead_code)]
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Deletion
                | Self::Skip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }
}

pub trait NoodlesAlignmentLike {}
impl NoodlesAlignmentLike for noodles_sam::alignment::record_buf::RecordBuf {}

/// implement the AlnRecordLike trait for the underlying noodles Record type
impl<T: NoodlesAlignmentLike + noodles_sam::alignment::Record> AlnRecordLike for T {
    fn opt_sequence_len(&self) -> Option<usize> {
        Some(self.sequence().len())
    }

    fn is_reverse_complemented(&self) -> bool {
        self.flags().expect("valid flags").is_reverse_complemented()
    }

    fn is_unmapped(&self) -> bool {
        self.flags()
            .expect("alignment record should have flags")
            .is_unmapped()
    }

    fn ref_id(&self, header: &Header) -> anyhow::Result<usize> {
        self.reference_sequence_id(header)
            .unwrap()
            .map_err(anyhow::Error::from)
    }

    fn aln_span(&self) -> Option<usize> {
        // NOTE: The transpose here is because noodles_sam
        // switched from returning Result<Option<usize>> to
        // Option<Result<usize>> --- think if there is a
        // better way to handle the new return type directly
        self.alignment_span().transpose().expect("valid span")
    }

    fn aln_score(&self) -> Option<i64> {
        self.data()
            .get(&AlnTag::ALIGNMENT_SCORE)
            .unwrap()
            .expect("could not get value")
            .as_int()
    }

    fn aln_start(&self) -> u32 {
        self.alignment_start()
            .unwrap()
            .expect("valid aln start")
            .get() as u32
    }

    fn aln_end(&self) -> u32 {
        self.alignment_end()
            .unwrap()
            .expect("valid aln start")
            .get() as u32
    }

    fn is_supp(&self) -> bool {
        self.flags()
            .expect("alignment record should have flags")
            .is_supplementary()
    }

    fn name(&self) -> Option<String> {
        self.name().map(|n| n.to_string())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AlnInfo {
    pub ref_id: u32,
    pub start: u32,
    pub end: u32,
    pub prob: f64,
    pub strand: Strand,
}

impl AlnInfo {
    #[allow(dead_code)]
    pub fn alignment_span(&self) -> u32 {
        self.end - self.start
    }
}

impl AlnInfo {
    fn from_aln_rec_like<T: AlnRecordLike>(aln: &T, aln_header: &Header) -> Self {
        Self {
            ref_id: aln.ref_id(aln_header).expect("valid ref_id") as u32,
            start: aln.aln_start(),
            end: aln.aln_end(),
            prob: 0.0_f64,
            strand: if aln.is_reverse_complemented() {
                Strand::Reverse
            } else {
                Strand::Forward
            },
        }
    }
}
/*
impl<T: sam::alignment::record::Record> From<&T> for AlnInfo {
    fn from(aln: &T) -> Self {
        Self {
            ref_id: aln.reference_sequence_id().unwrap() as u32,
            start: aln.alignment_start().unwrap().expect("valid aln start").get() as u32,
            end: aln.alignment_end().unwrap().expect("valid aln end").get() as u32,
            prob: 0.0_f64,
            strand: if aln.flags().expect("valid flags").is_reverse_complemented() {
                Strand::Reverse
            } else {
                Strand::Forward
            },
        }
    }
}
*/

/// Holds the info relevant for reading the short read quants
#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "PascalCase")]
pub struct ShortReadRecord {
    pub name: String,
    #[allow(dead_code)]
    pub length: i32,
    #[allow(dead_code)]
    pub effective_length: f64,
    #[allow(dead_code)]
    #[serde(rename = "TPM")]
    pub tpm: f64,
    pub num_reads: f64,
}

#[allow(dead_code)]
impl ShortReadRecord {
    pub fn empty(name: &str) -> Self {
        Self {
            name: name.to_owned(),
            length: 0,
            effective_length: 0.0,
            tpm: 0.0,
            num_reads: 0.0,
        }
    }
}

/// Holds the info relevant for running the EM algorithm
pub struct EMInfo<'eqm, 'tinfo, 'h> {
    // the read alignment infomation we'll need to
    // perform the EM
    pub eq_map: &'eqm InMemoryAlignmentStore<'h>,
    // relevant information about each target transcript
    pub txp_info: &'tinfo [TranscriptInfo],
    // maximum number of iterations the EM will run
    // before returning an estimate.
    pub max_iter: u32,
    // the EM will terminate early if *all* parameters
    // have converged within this threshold of relative
    // change between two subsequent iterations.
    pub convergence_thresh: f64,
    // An optional vector of abundances from which
    // to initalize the EM, otherwise, a default
    // uniform initalization is used.
    pub init_abundances: Option<Vec<f64>>,
    /// holds the KDE model if we will be using one
    /// and [None] otherwise
    pub kde_model: Option<KDEModel>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TranscriptInfo {
    pub len: NonZeroUsize,
    pub total_weight: f64,
    pub coverage_bins: Vec<f64>,
    pub coverage_prob: Vec<f64>,
    pub lenf: f64,
}

impl TranscriptInfo {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            len: NonZeroUsize::new(0).unwrap(),
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; 10],
            coverage_prob: Vec::new(),
            lenf: 0_f64,
        }
    }

    pub fn with_len(len: NonZeroUsize) -> Self {
        Self {
            len,
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; 10],
            coverage_prob: Vec::new(),
            lenf: len.get() as f64,
        }
    }
    pub fn with_len_and_bin_width(len: NonZeroUsize, bin_width: u32) -> Self {
        Self {
            len,
            total_weight: 0.0_f64,
            coverage_bins: vec![0.0_f64; ((len.get() as f64) / (bin_width as f64)).ceil() as usize],
            coverage_prob: Vec::new(),
            lenf: len.get() as f64,
        }
    }

    #[inline(always)]
    pub fn get_normalized_counts_and_lengths(&self) -> (Vec<f32>, Vec<f32>) {
        let num_intervals = self.coverage_bins.len();
        let num_intervals_f = num_intervals as f64;
        let tlen_f = self.lenf;
        let bin_width = (tlen_f / num_intervals_f).round() as f32;

        let cov_f32 = self.coverage_bins.iter().map(|f| *f as f32).collect();
        let mut widths_f32 = Vec::<f32>::with_capacity(num_intervals);
        for bidx in 0..num_intervals {
            let bidxf = bidx as f32;
            let bin_start = bidxf * bin_width;
            let bin_end = ((bidxf + 1.0) * bin_width).min(self.lenf as f32);
            widths_f32.push(bin_end - bin_start);
            if bin_end < bin_start {
                warn!("tlen_f: {:?}", tlen_f);
                warn!("num_intervals_f: {:?}", num_intervals_f);
                warn!("bin_width: {:?}", bin_width);
                warn!("bidxf: {:?}", bidxf);
            }
            assert!(bin_end > bin_start);
        }
        (cov_f32, widths_f32)
    }

    #[inline(always)]
    pub fn add_interval(&mut self, start: u32, stop: u32, weight: f64) {
        const ONE_PLUS_EPSILON: f64 = 1.0_f64 + f64::EPSILON;
        let num_intervals = self.coverage_bins.len();
        let num_intervals_f = num_intervals as f64;
        let tlen_f = self.lenf;
        let bin_width = (tlen_f / num_intervals_f).round();
        let start = start.min(stop);
        let stop = start.max(stop);
        let start_bin = (((start as f64) / tlen_f) * num_intervals_f).floor() as usize;
        let end_bin = (((stop as f64) / tlen_f) * num_intervals_f).floor() as usize;

        let get_overlap = |s1: u32, e1: u32, s2: u32, e2: u32| -> u32 {
            if s1 <= e2 {
                e1.min(e2) - s1.max(s2)
            } else {
                0_u32
            }
        };

        for (bidx, bin) in self.coverage_bins[start_bin..end_bin]
            .iter_mut()
            .enumerate()
        {
            let bidxf = (start_bin + bidx) as f64;
            let curr_bin_start = (bidxf * bin_width) as u32;
            let curr_bin_end = ((bidxf + 1.0) * bin_width).min(tlen_f) as u32;

            let olap = get_overlap(start, stop, curr_bin_start, curr_bin_end);
            let olfrac = (olap as f64) / ((curr_bin_end - curr_bin_start) as f64);
            *bin += olfrac;
            if olfrac > ONE_PLUS_EPSILON {
                error!("first_bin = {start_bin}, last_bin = {end_bin}");
                error!(
                    "bin = {}, olfrac = {}, olap = {}, curr_bin_start = {}, curr_bin_end = {}, start = {start}, stop = {stop}",
                    *bin, olfrac, olap, curr_bin_start, curr_bin_end
                );
                panic!(
                    "coverage computation error; please report this error at https://github.com/COMBINE-lab/oarfish."
                )
            }
        }
        self.total_weight += weight;
    }

    #[allow(dead_code)]
    pub fn clear_coverage_dist(&mut self) {
        self.coverage_bins.fill(0.0_f64);
        self.total_weight = 0.0_f64;
    }
}

#[derive(Debug)]
pub struct InMemoryAlignmentStore<'h> {
    pub filter_opts: AlignmentFilters,
    pub aln_header: &'h Header,
    pub alignments: Vec<AlnInfo>,
    pub as_probabilities: Vec<f32>,
    pub coverage_probabilities: Vec<f64>,
    // holds the boundaries between records for different reads
    boundaries: Vec<usize>,
    pub discard_table: DiscardTable,
    pub num_unique_alignments: usize,
}

impl InMemoryAlignmentStore<'_> {
    #[inline]
    pub fn len(&self) -> usize {
        self.boundaries.len().saturating_sub(1)
    }

    pub fn aggregate_discard_table(&mut self, table: &DiscardTable) {
        self.discard_table.aggregate(table);
    }
}

pub struct InMemoryAlignmentStoreSamplingWithReplacementIter<'a, 'h, 'b> {
    pub store: &'a InMemoryAlignmentStore<'h>,
    pub rand_inds: std::slice::Iter<'b, usize>,
}

impl<'a> Iterator for InMemoryAlignmentStoreSamplingWithReplacementIter<'a, '_, '_> {
    type Item = (&'a [AlnInfo], &'a [f32], &'a [f64]);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next_ind) = self.rand_inds.next() {
            let start = self.store.boundaries[*next_ind];
            let end = self.store.boundaries[*next_ind + 1];
            Some((
                &self.store.alignments[start..end],
                &self.store.as_probabilities[start..end],
                &self.store.coverage_probabilities[start..end],
            ))
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.rand_inds.size_hint()
    }
}

impl ExactSizeIterator for InMemoryAlignmentStoreSamplingWithReplacementIter<'_, '_, '_> {}

pub struct InMemoryAlignmentStoreIter<'a, 'h> {
    pub store: &'a InMemoryAlignmentStore<'h>,
    pub idx: usize,
}

impl<'a> Iterator for InMemoryAlignmentStoreIter<'a, '_> {
    type Item = (&'a [AlnInfo], &'a [f32], &'a [f64]);

    #[inline]
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

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (
            self.store.len() - self.idx,
            Some(self.store.len() - self.idx),
        )
    }
}

impl ExactSizeIterator for InMemoryAlignmentStoreIter<'_, '_> {}

impl<'h> InMemoryAlignmentStore<'h> {
    pub fn new(fo: AlignmentFilters, header: &'h Header) -> Self {
        InMemoryAlignmentStore {
            filter_opts: fo.clone(),
            aln_header: header,
            alignments: vec![],
            as_probabilities: vec![],
            coverage_probabilities: vec![],
            boundaries: vec![0],
            discard_table: DiscardTable::new(),
            num_unique_alignments: 0,
        }
    }

    pub fn iter(&self) -> InMemoryAlignmentStoreIter<'_, '_> {
        InMemoryAlignmentStoreIter {
            store: self,
            idx: 0,
        }
    }

    pub fn random_sampling_iter<'a, 'b>(
        &'a self,
        inds: &'b [usize],
    ) -> InMemoryAlignmentStoreSamplingWithReplacementIter<'a, 'h, 'b>
    where
        'b: 'a,
    {
        InMemoryAlignmentStoreSamplingWithReplacementIter {
            store: self,
            rand_inds: inds.iter(),
        }
    }

    #[inline(always)]
    pub fn add_group<T: NoodlesAlignmentLike + sam::alignment::record::Record + std::fmt::Debug>(
        &mut self,
        txps: &mut [TranscriptInfo],
        ag: &mut Vec<T>,
    ) -> bool {
        if !ag.is_empty() {
            let (alns, as_probs) =
                self.filter_opts
                    .filter(&mut self.discard_table, self.aln_header, txps, ag);
            self.add_filtered_group(&alns, &as_probs, txps)
        } else {
            false
        }
    }

    /// Filter and add one read's projected (transcriptome-space) alignments.
    ///
    /// This is the genome-mode analogue of [`add_group`](Self::add_group): it
    /// runs [`AlignmentFilters::filter_projected`] (updating the discard table)
    /// and then [`add_filtered_group`](Self::add_filtered_group). `read_len` is
    /// the query length (for the aligned-fraction filter) and `beta` the
    /// similarity→probability spread.
    #[inline(always)]
    pub fn add_projected_group(
        &mut self,
        txps: &mut [TranscriptInfo],
        recs: &[ProjectedAlnRecord],
        read_len: usize,
        beta: f32,
        prob_source: ProjProbSource,
    ) -> bool {
        if recs.is_empty() {
            return false;
        }
        let (alns, as_probs) = self.filter_opts.filter_projected(
            &mut self.discard_table,
            txps,
            recs,
            read_len,
            beta,
            prob_source,
        );
        self.add_filtered_group(&alns, &as_probs, txps)
    }

    #[inline(always)]
    pub fn add_filtered_group(
        &mut self,
        alns: &[AlnInfo],
        as_probs: &[f32],
        txps: &mut [TranscriptInfo],
    ) -> bool {
        if !alns.is_empty() {
            for a in alns.iter() {
                let tid = a.ref_id as usize;
                txps[tid].add_interval(a.start, a.end, 1.0_f64);
            }
            self.alignments.extend_from_slice(alns);
            self.as_probabilities.extend_from_slice(as_probs);
            self.coverage_probabilities
                .extend(vec![0.0_f64; alns.len()]);
            self.boundaries.push(self.alignments.len());
            true
        } else {
            false
        }
    }

    #[inline(always)]
    pub fn total_len(&self) -> usize {
        self.alignments.len()
    }

    #[inline(always)]
    pub fn num_aligned_reads(&self) -> usize {
        self.len()
    }

    #[inline(always)]
    pub fn inc_unique_alignments(&mut self) {
        self.num_unique_alignments += 1;
    }

    #[inline(always)]
    pub fn unique_alignments(&self) -> usize {
        self.num_unique_alignments
    }
}

/// The parameters controling the filters that will
/// be applied to alignments
#[derive(TypedBuilder, Clone, Debug, Serialize)]
pub struct AlignmentFilters {
    /// How far an alignment can start from the
    /// 5' end of the transcript and still be
    /// considered valid
    five_prime_clip: u32,
    /// How far an alignment can end from the
    /// 3' end of the transcript and still be
    /// considered valid
    three_prime_clip: i64,
    /// The fraction of the best recorded alignment
    /// score for this *read* that the current alignment
    /// must obtain in order to be retained. For
    /// example 0.95 means that it must be within
    /// 95% of the best alignment.
    score_threshold: f32,
    /// The minimum fraction of the read length that
    /// must be part of the alignment in order
    /// for this alignment to be retained.
    min_aligned_fraction: f32,
    /// The minimum absolute length of this read that
    /// must be aligned under this alignment in
    /// order for this alignment to be  retained
    min_aligned_len: u32,
    /// Determines which alignments we should consider
    /// as valid during quantificationwhich_strand.
    which_strand: bio_types::strand::Strand,
    // True if we are enabling our coverage model and
    // false otherwise.
    pub model_coverage: bool,
    // The growth rate (or `k`) parameter of the logistic
    // function. This only matters if `model_coverage` is true.
    pub logistic_growth_rate: f64,
    // True if we are enabling to output the alignment probability and
    // false otherwise.
    pub write_assignment_probs: bool,
    pub write_assignment_probs_type: Option<ReadAssignmentProbOut>,
    /// Denominator `D` in the score→probability conversion
    /// `exp((score - best_score) / D)` used to weight a read's alignments in the
    /// EM (transcriptome mode and the `score`/`combined` projected sources).
    /// Larger flattens, smaller sharpens. Exposed via `--score-prob-denom`.
    #[builder(default = 5.0)]
    pub score_prob_denom: f32,
}

/// This structure records information about
/// the number of alignments (and reads) discarded
/// due to the application of `AlignmentFilters`.
#[derive(Debug, Serialize)]
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
    pub fn new() -> Self {
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

    pub fn aggregate(&mut self, other: &Self) {
        self.discard_5p += other.discard_5p;
        self.discard_3p += other.discard_3p;
        self.discard_score += other.discard_score;
        self.discard_aln_frac += other.discard_aln_frac;
        self.discard_aln_len += other.discard_aln_len;
        self.discard_ori += other.discard_ori;
        self.discard_supp += other.discard_supp;
        self.valid_best_aln += other.valid_best_aln;
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
        writeln!(
            f,
            "discarded because of distance from 5' {}",
            self.discard_5p
        )
        .expect("couldn't format discard table.");
        writeln!(
            f,
            "discarded because of distance from 3' {}",
            self.discard_3p
        )
        .expect("couldn't format discard table.");
        writeln!(
            f,
            "discarded because of score fraction {}",
            self.discard_score
        )
        .expect("couldn't format discard table.");
        writeln!(
            f,
            "discarded because of aligned fraction {}",
            self.discard_aln_frac
        )
        .expect("couldn't format discard table.");
        writeln!(
            f,
            "discarded because of aligned length {}",
            self.discard_aln_len
        )
        .expect("couldn't format discard table.");
        writeln!(
            f,
            "discarded because of aligned orientation {}",
            self.discard_ori
        )
        .expect("couldn't format discard table.");
        writeln!(
            f,
            "discarded because alignment is supplemental {}",
            self.discard_supp
        )
    }
}

impl AlignmentFilters {
    /// Applies the filters defined by this AlignmentFilters struct
    /// to the alignments provided in `ag`, a vector of alignments representing
    /// a group of contiguous alignments for the same target.
    ///
    /// The details about what alignments have been filtered and why will
    /// be added to the provided `discard_table`.  
    ///
    /// This function returns a vector of the `AlnInfo` structs for alignments
    /// that pass the filter, the associated probabilities for each and, if the
    /// user requested per-read alignment probabilities, the read name.
    pub fn filter<T: AlnRecordLike + std::fmt::Debug>(
        &mut self,
        discard_table: &mut DiscardTable,
        aln_header: &Header,
        txps: &[TranscriptInfo],
        ag: &mut Vec<T>,
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
            .find_map(|x| x.opt_sequence_len().map(|y| y as u32))
            .unwrap_or(0_u32);

        // apply the filter criteria to determine what alignments to retain
        ag.retain(|x| {
            // we ony want to retain mapped reads
            if !x.is_unmapped() {
                let tid = x.ref_id(aln_header).expect("valid ref id");

                // get the alignment span
                let aln_span = x.aln_span().unwrap() as u32;

                // get the alignment score, as computed by the aligner
                let score = x.aln_score().unwrap_or(i32::MIN as i64) as i32;

                // the alignment is to the - strand
                let is_rc = x.is_reverse_complemented();

                // if we are keeping only forward strand alignments
                // filter this alignment if it is rc
                match (is_rc, self.which_strand) {
                    (_, bio_types::strand::Strand::Unknown) => { /*do nothing*/ }
                    // is rc and we want rc
                    (true, bio_types::strand::Strand::Reverse) => { /*do nothing*/ }
                    // is fw and we want fw
                    (false, bio_types::strand::Strand::Forward) => { /*do nothing*/ }
                    // is fw and we want rc
                    (false, bio_types::strand::Strand::Reverse) => {
                        discard_table.discard_ori += 1;
                        return false;
                    }
                    // is rc and we want fw
                    (true, bio_types::strand::Strand::Forward) => {
                        discard_table.discard_ori += 1;
                        return false;
                    }
                }

                // the alignment is supplementary
                // *NOTE*: this removes "supplementary" alignments, *not*
                // "secondary" alignments.
                let is_supp = x.is_supp();
                if is_supp {
                    discard_table.discard_supp += 1;
                    return false;
                }

                // enough absolute sequence (# of bases) is aligned
                let filt_aln_len = aln_span < self.min_aligned_len;
                if filt_aln_len {
                    discard_table.discard_aln_len += 1;
                    return false;
                }

                // not too far from the 3' end
                let filt_3p =
                    (x.aln_end() as i64) <= (txps[tid].len.get() as i64 - self.three_prime_clip);
                if filt_3p {
                    discard_table.discard_3p += 1;
                    return false;
                }

                // not too far from the 5' end
                let filt_5p = x.aln_start() >= self.five_prime_clip;
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
            .map(|a| a.aln_score().unwrap_or(0) as i32)
            .collect();

        let _min_allowed_score = self.score_threshold * mscore;

        let score_prob_denom = self.score_prob_denom;
        for score in scores.iter_mut() {
            let fscore = *score as f32;
            let score_ok = (fscore * inv_max_score) >= self.score_threshold; //>= thresh_score;
            if score_ok {
                let f = (fscore - mscore) / score_prob_denom;
                probabilities.push(f.exp());
            } else {
                *score = i32::MIN;
                discard_table.discard_score += 1;
            }
        }

        let mut score_it = scores.iter();
        ag.retain(|_| *score_it.next().unwrap() > i32::MIN);
        assert_eq!(ag.len(), probabilities.len());

        (
            ag.iter()
                .map(|x| AlnInfo::from_aln_rec_like(x, aln_header))
                .collect(),
            probabilities,
        )
    }
}

/// A single transcriptome-space alignment produced by projecting a genomic
/// alignment through bramble. This is the neutral hand-off type between the
/// projection bridge (`util::projection`) and the alignment filters, so that
/// [`AlignmentFilters`] need not depend on bramble's `ProjectedAlignment`.
///
/// Coordinates follow the same convention as the BAM path: `start`/`end` are
/// 1-based (matching noodles `alignment_start`/`alignment_end`), and
/// `aligned_len` is the number of transcript bases the alignment spans (the
/// analogue of `aln_span`).
#[derive(Clone, Debug)]
pub struct ProjectedAlnRecord {
    /// Transcript id == reference index in the (transcriptome) header.
    pub ref_id: u32,
    /// 1-based start on the transcript (forward-strand coordinates).
    pub start: u32,
    /// 1-based inclusive end on the transcript (forward-strand coordinates).
    pub end: u32,
    /// Transcript bases spanned by the alignment (`end - start + 1`).
    pub aligned_len: u32,
    /// Query (read) bases in the aligned portion (M/I/=/X), for aligned-fraction.
    pub query_aligned_len: u32,
    /// True if the read maps to the reverse strand of the transcript.
    pub is_reverse: bool,
    /// bramble similarity score (higher is better; best in a group anchors the
    /// probability), used in place of the integer alignment score.
    pub similarity: f64,
    /// alignment score of the *source genomic alignment* this was
    /// projected from (shared by all transcripts at the same genomic locus).
    /// Used by the `score`/`combined` probability sources to discriminate
    /// paralogous loci, which similarity cannot. 0 if unavailable.
    pub aln_score: i32,
}

impl AlignmentFilters {
    /// Filter a single read's projected (transcriptome-space) alignments and
    /// compute per-alignment probabilities.
    ///
    /// This mirrors [`AlignmentFilters::filter`] but operates on projected
    /// alignments (which have no integer alignment score). Probabilities are
    /// derived from bramble's similarity score relative to the read's best
    /// projected hit: `prob = exp((similarity - best_similarity) * beta)`, so
    /// the best hit gets weight 1.0 and weaker hits decay (the analogue of the
    /// BAM path's `exp((score - max_score) / 5.0)`). `read_len` is the query
    /// length used for the aligned-fraction filter; `beta` is the spread knob
    /// (`--projected-prob-beta`). Discard reasons are recorded in
    /// `discard_table` exactly as the BAM path does.
    pub fn filter_projected(
        &self,
        discard_table: &mut DiscardTable,
        txps: &[TranscriptInfo],
        recs: &[ProjectedAlnRecord],
        read_len: usize,
        beta: f32,
        prob_source: ProjProbSource,
    ) -> (Vec<AlnInfo>, Vec<f32>) {
        let mut best_sim = f64::MIN;
        let mut best_score = i32::MIN;
        let mut aln_frac_at_best = 0_f32;
        let mut kept: Vec<&ProjectedAlnRecord> = Vec::with_capacity(recs.len());

        for r in recs.iter() {
            // strand orientation
            match (r.is_reverse, self.which_strand) {
                (_, Strand::Unknown) => { /* keep both */ }
                (true, Strand::Reverse) => { /* keep */ }
                (false, Strand::Forward) => { /* keep */ }
                (false, Strand::Reverse) | (true, Strand::Forward) => {
                    discard_table.discard_ori += 1;
                    continue;
                }
            }

            // enough absolute sequence (# of transcript bases) is aligned
            if r.aligned_len < self.min_aligned_len {
                discard_table.discard_aln_len += 1;
                continue;
            }

            let tid = r.ref_id as usize;

            // not too far from the 3' end
            if (r.end as i64) <= (txps[tid].len.get() as i64 - self.three_prime_clip) {
                discard_table.discard_3p += 1;
                continue;
            }

            // not too far from the 5' end
            if r.start >= self.five_prime_clip {
                discard_table.discard_5p += 1;
                continue;
            }

            // committed to retaining; track the best (highest-similarity) hit
            if r.similarity > best_sim {
                best_sim = r.similarity;
                aln_frac_at_best = if read_len > 0 {
                    (r.query_aligned_len as f32) / (read_len as f32)
                } else {
                    0_f32
                };
            }
            if r.aln_score > best_score {
                best_score = r.aln_score;
            }
            kept.push(r);
        }

        if kept.is_empty() || best_sim <= 0.0 {
            return (vec![], vec![]);
        }
        if aln_frac_at_best < self.min_aligned_fraction {
            discard_table.discard_aln_frac += 1;
            return (vec![], vec![]);
        }

        discard_table.valid_best_aln += 1;

        let msim = best_sim;
        let inv_msim = 1.0 / msim;
        let mut out_alns = Vec::<AlnInfo>::with_capacity(kept.len());
        let mut probabilities = Vec::<f32>::with_capacity(kept.len());

        for r in kept {
            let sim_ok = ((r.similarity * inv_msim) as f32) >= self.score_threshold;
            if !sim_ok {
                discard_table.discard_score += 1;
                continue;
            }
            // Defensively clamp the projected interval to the transcript bounds.
            // With a correct projection `1 <= start <= end <= len` always holds;
            // this guard just prevents a stray out-of-range coordinate from
            // panicking the coverage binning during a long run.
            let tlen = txps[r.ref_id as usize].len.get() as u32;
            let start = r.start.clamp(1, tlen);
            let end = r.end.clamp(start, tlen);

            // log-weight of this alignment relative to the read's best hit.
            // `score` discriminates paralogous genomic loci (similarity saturates
            // near 1.0 for any well-covered transcript and cannot); `similarity`
            // discriminates isoforms sharing a locus (same genomic score).
            let score_prob_denom = self.score_prob_denom;
            let f = match prob_source {
                ProjProbSource::Similarity => ((r.similarity - msim) as f32) * beta,
                ProjProbSource::Score => ((r.aln_score - best_score) as f32) / score_prob_denom,
                ProjProbSource::Combined => {
                    ((r.aln_score - best_score) as f32) / score_prob_denom
                        + beta * ((r.similarity - msim) as f32)
                }
            };
            probabilities.push(f.exp());
            out_alns.push(AlnInfo {
                ref_id: r.ref_id,
                start,
                end,
                prob: 0.0_f64,
                strand: if r.is_reverse {
                    Strand::Reverse
                } else {
                    Strand::Forward
                },
            });
        }

        (out_alns, probabilities)
    }
}

#[cfg(test)]
mod tests {
    use crate::util::oarfish_types::AlnInfo;
    use bio_types::strand::Strand;

    #[test]
    fn aln_span_is_correct() {
        let ainf = AlnInfo {
            ref_id: 0,
            start: 0,
            end: 100,
            prob: 0.5,
            strand: Strand::Forward,
        };
        assert_eq!(ainf.alignment_span(), 100);
    }
}
