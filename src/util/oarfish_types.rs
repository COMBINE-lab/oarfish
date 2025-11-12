use serde::Deserialize;
use std::fmt;
use std::num::NonZeroUsize;

use statrs::distribution::{Continuous, Normal};
use statrs::statistics::Distribution;

use kders::kde::KDEModel;
use std::iter::FromIterator;
use tabled::builder::Builder;
use tabled::settings::Style;

use seqcol_rs::DigestResult;
use serde::Serialize;
use serde_json::json;
use serde_with::serde_as;
use typed_builder::TypedBuilder;

use bio_types::strand::Strand;
use bstr::{B, ByteSlice};
//use minimap2_temp as minimap2;
use minimap2;
use noodles_sam as sam;
use sam::{Header, alignment::record::data::field::tag::Tag as AlnTag};

#[allow(unused_imports)]
use tracing::{error, info, warn};

use crate::prog_opts::{FragmentEndModel, ReadAssignmentOut, ReadAssignmentProbOut};
use crate::util::constants::EMPTY_READ_NAME;

use crate::util::probs::LogSpace;

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

    pub fn iter(&self) -> ReadChunkIter {
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

impl AlnRecordLike for minimap2::Mapping {
    fn opt_sequence_len(&self) -> Option<usize> {
        self.query_len.map(|x| x.get() as usize)
    }

    fn is_reverse_complemented(&self) -> bool {
        self.strand == minimap2::Strand::Reverse
    }

    fn is_unmapped(&self) -> bool {
        self.target_name.is_none()
    }

    fn ref_id(&self, _header: &Header) -> anyhow::Result<usize> {
        if let Some(ref tgt_name) = self.target_name {
            let tgt_str = tgt_name.as_bytes();
            if let Some(id) = _header.reference_sequences().get_index_of(tgt_str) {
                return Ok(id);
            }
            anyhow::bail!("Could not get ref_id of target {}", tgt_name);
        }
        anyhow::bail!("Could not get ref_id of mapping without target name")
    }

    fn aln_span(&self) -> Option<usize> {
        if let Some(ref aln) = self.alignment {
            return if let Some(ref cigar) = aln.cigar {
                let mut span = 0_usize;
                for (len, op) in cigar.iter() {
                    let co: CigarOp = (*op).into();
                    if co.consumes_reference() {
                        span += *len as usize;
                    }
                }
                Some(span)
            } else {
                error!("Had an alignment but no CIGAR!");
                None
            };
        }
        None
    }

    fn aln_score(&self) -> Option<i64> {
        if let Some(ref aln) = self.alignment {
            aln.alignment_score.map(|x| x as i64)
        } else {
            None
        }
    }

    fn aln_start(&self) -> u32 {
        self.target_start as u32
    }

    fn aln_end(&self) -> u32 {
        self.target_end as u32
    }

    fn is_supp(&self) -> bool {
        self.is_supplementary
    }

    fn name(&self) -> Option<String> {
        self.query_name.as_ref().map(|q| q.to_string())
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
    pub eq_map: &'eqm mut InMemoryAlignmentStore<'h>,
    // relevant information about each target transcript
    pub txp_info: &'tinfo mut [TranscriptInfo],
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

impl<'eqm, 'tinfo, 'h> EMInfo<'eqm, 'tinfo, 'h> {
    pub fn num_txps(&self) -> usize {
        self.txp_info.len()
    }
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
    pub fn get_normalized_counts_and_lengths(&self, verbose: bool) -> (Vec<f32>, Vec<f32>) {
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
        if verbose {
            eprintln!(
                "num_intervals: {}, bin_width: {}, cov: {:#?}",
                num_intervals_f, bin_width, cov_f32
            );
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
            *bin += weight * olfrac;
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

impl InMemoryAlignmentStore<'_> {
    #[inline]
    pub fn get_alignments_for_read(
        &self,
        read_rank: usize,
    ) -> anyhow::Result<(&'_ [AlnInfo], &'_ [f32], &'_ [f64])> {
        if read_rank + 1 < self.boundaries.len() {
            let start = self.boundaries[read_rank];
            let end = self.boundaries[read_rank + 1];
            Ok((
                &self.alignments[start..end],
                &self.as_probabilities[start..end],
                &self.coverage_probabilities[start..end],
            ))
        } else {
            anyhow::bail!(
                "asked for read of rank {read_rank}, but there is only information for {} reads.",
                self.boundaries.len()
            )
        }
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

    pub fn iter(&self) -> InMemoryAlignmentStoreIter {
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

    #[inline(always)]
    pub fn add_filtered_group(
        &mut self,
        alns: &[AlnInfo],
        as_probs: &[f32],
        _txps: &mut [TranscriptInfo],
    ) -> bool {
        if !alns.is_empty() {
            /*
            for a in alns.iter() {
                let tid = a.ref_id as usize;
                txps[tid].add_interval(a.start, a.end, 1.0_f64);
            }
            */
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

serde_with::serde_conv!(
    DistrAsPair,
    statrs::distribution::Normal,
    |dist: &statrs::distribution::Normal| (
        dist.mean().expect("valid mean"),
        dist.std_dev().expect("valid std_dev")
    ),
    |value: (f64, f64)| -> Result<_, std::convert::Infallible> {
        Ok(statrs::distribution::Normal::new(value.0, value.1).expect("can construct"))
    }
);

/// The model that will be applied to reads starting far from the sequenced end of the
/// underlying transcript
#[serde_as]
#[derive(Clone, Debug, Serialize)]
pub struct FragmentEndFalloffDist {
    #[serde_as(as = "DistrAsPair")]
    pub dist: statrs::distribution::Normal,
    pub thresh: f64,
    pub model: FragmentEndModel,
    pub log_min_prob: f64,
}

// inspired by https://docs.rs/bio-types/latest/src/bio_types/strand.rs.html#81
fn compatible_ori<T: AlnRecordLike>(a: &T, s: &bio_types::strand::Strand) -> bool {
    match (a.is_reverse_complemented(), *s) {
        (true, Strand::Forward) => false,
        (false, Strand::Forward) => true,
        (true, Strand::Reverse) => true,
        (false, Strand::Reverse) => false,
        (_, Strand::Unknown) => true,
    }
}

impl FragmentEndFalloffDist {
    pub fn new(
        mu: f64,
        std_dev: f64,
        thresh: f64,
        model: FragmentEndModel,
        log_min_prob: f64,
    ) -> Self {
        FragmentEndFalloffDist {
            dist: Normal::new(mu, std_dev)
                .expect("should be able to construct a normal distribution"),
            thresh,
            model,
            log_min_prob,
        }
    }

    pub fn eval_alignment_ln<T: AlnRecordLike>(&self, a: &T, tlenf: f64) -> LogSpace {
        let thresh = self.thresh;
        let extra_dist = match self.model {
            FragmentEndModel::ThreePrimeStart(ori) => {
                //           end      clip   txp_end
                // ========== * ====== ( =======]
                //               dist
                if compatible_ori::<T>(a, &ori) {
                    ((tlenf - thresh) - (a.aln_end() as f64)).max(0.)
                } else {
                    tlenf
                }
            }
            FragmentEndModel::FivePrimeStart(ori) => {
                if compatible_ori::<T>(a, &ori) {
                    (a.aln_start() as f64 - thresh).max(0.)
                } else {
                    tlenf
                }
            }
            FragmentEndModel::EitherStart(ori) => {
                // TODO: It would be nice if we can make an informed choice
                // here based on alignment orientation, but right now we
                // will just chose the alignment part closer to the end
                if compatible_ori::<T>(a, &ori) {
                    let three_prime_dist = ((tlenf - thresh) - (a.aln_end() as f64)).max(0.);
                    let five_prime_dist = (a.aln_start() as f64 - thresh).max(0.);
                    three_prime_dist.min(five_prime_dist)
                } else {
                    tlenf
                }
            }
            FragmentEndModel::None => 0.,
        };
        let p = self.dist.ln_pdf(extra_dist) - self.dist.ln_pdf(0.);
        let minp = self.log_min_prob;
        LogSpace::new_from_ln(p.max(minp))
    }

    /*
    pub fn eval_norm(&self, distance: i64) -> f64 {
        let cdist: f64 = distance.max(0) as f64;
        self.dist.pdf(cdist) / self.dist.pdf(0.)
    }

    pub fn eval_norm_thresh(&self, distance: i64, thresh: i64) -> f64 {
        let cdist: f64 = (distance - thresh).max(0) as f64;
        self.dist.pdf(cdist) / self.dist.pdf(0.)
    }
    */
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
    // the denominator in the heat kernel function used to downweight
    // suboptimal alignments (i.e. the `t` in the equation
    // exp^((A - Amax)/t)
    pub alignment_score_denom: f32,
    // To evaluate soft probabilities for fragments starting far
    // from the end.
    pub falloff_dist: Option<FragmentEndFalloffDist>,
    // True if we are enabling to output the alignment probability and
    // false otherwise.
    pub write_assignment_probs: bool,
    pub write_assignment_probs_type: Option<ReadAssignmentProbOut>,
    // True if we are enabling to output the alignment assignments
    // false otherwise.
    pub write_assignments: bool,
    pub write_assignments_type: Option<ReadAssignmentOut>,
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
                let filt_3p = (!is_rc)
                    && ((x.aln_end() as i64)
                        <= (txps[tid].len.get() as i64 - self.three_prime_clip));
                if filt_3p {
                    discard_table.discard_3p += 1;
                    return false;
                }

                // not too far from the 5' end
                let filt_5p = (is_rc) && (x.aln_start() >= self.five_prime_clip);
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
        let mut m_dist_score = crate::util::probs::MIN_LOG_P;
        let mut tot_dist_score = crate::util::probs::MIN_LOG_P;

        // get a vector of all of the scores
        let mut scores: Vec<(i32, LogSpace)> = ag
            .iter_mut()
            .map(|a| {
                // the alignment is to the - strand
                let p = if let Some(falloff_dist) = &self.falloff_dist {
                    let tid = a.ref_id(aln_header).expect("valid ref id");
                    let tlenf = txps[tid].lenf;
                    let p = falloff_dist.eval_alignment_ln::<T>(a, tlenf);
                    m_dist_score = m_dist_score.max(p);
                    p
                } else {
                    LogSpace::new_from_linear(1.0_f64)
                };
                tot_dist_score += p;
                (a.aln_score().unwrap_or(0) as i32, p)
            })
            .collect();

        scores.iter_mut().for_each(|x| x.1 /= tot_dist_score);

        let _min_allowed_score = self.score_threshold * mscore;

        for score in scores.iter_mut() {
            let score_prob_denom = self.alignment_score_denom;
            let fscore = score.0 as f32;
            let score_ok = (fscore * inv_max_score) >= self.score_threshold; //>= thresh_score;
            if score_ok {
                //let f = ((fscore - mscore) / (mscore - min_allowed_score)) * SCORE_PROB_DENOM;
                let f = (fscore - mscore) / score_prob_denom;
                probabilities.push(f.exp() * score.1.get_linear() as f32);
            } else {
                score.0 = i32::MIN;
                discard_table.discard_score += 1;
            }
        }

        let mut score_it = scores.iter();
        ag.retain(|_| score_it.next().unwrap().0 > i32::MIN);
        assert_eq!(ag.len(), probabilities.len());

        (
            ag.iter()
                .map(|x| AlnInfo::from_aln_rec_like(x, aln_header))
                .collect(),
            probabilities,
        )
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
