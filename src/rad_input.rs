//! Bulk quantification from a RAD file of `bulk_long_read` records (one record per read with
//! all of its alignments; see the `long-read-rad` crate), as `spearfish map --rad` writes.
//! Any mapper that writes the schema can be quantified this way.
//!
//! The file's tags stand in for the BAM header: the target names and lengths, and the mapper's
//! suggested score scale (used for `--score-prob-denom` when it is not given).

use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Mutex;

use crossbeam::channel::{Receiver, Sender, bounded};
use long_read_rad::{Aln, Header as RadHeader};
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;
use num_format::{Locale, ToFormattedString};
use tracing::{info, warn};

use crate::prog_opts::Args;
use crate::util::oarfish_types::{
    AlignmentFilters, AlnInfo, AlnRecordLike, DiscardTable, InMemoryAlignmentStore, NamedDigestVec,
    TranscriptInfo,
};

/// RAD scores are real-valued and oarfish's filter takes integers: scores are scaled by this
/// factor, as for the in-process spearfish backend.
pub const RAD_SCORE_SCALE: f32 = crate::util::mapper::SPEARFISH_SCORE_SCALE;

/// Whether `path` is a RAD file rather than a BAM file.
pub fn is_rad(path: &Path) -> anyhow::Result<bool> {
    long_read_rad::is_rad(path)
}

/// Read the file's header and fill in the options it implies: the score denominator, if the
/// user gave none and the mapper suggests one.
pub fn read_header(path: &Path, args: &mut Args) -> anyhow::Result<RadHeader> {
    let h = long_read_rad::read_header(path)?;
    let m = &h.meta;
    info!(
        "RAD file {}: {} targets; {} of {} reads mapped by {} {} (preset {:?}, chunk codec {})",
        path.display(),
        h.ref_names.len().to_formatted_string(&Locale::en),
        h.counts.mapped.to_formatted_string(&Locale::en),
        h.counts.reads.to_formatted_string(&Locale::en),
        m.mapper,
        m.mapper_version,
        m.preset,
        h.codec.as_str()
    );
    if args.score_prob_denom.is_none() && m.score_denom > 0.0 {
        let d = m.score_denom * RAD_SCORE_SCALE;
        info!(
            "using the mapper's score scale: --score-prob-denom {d} (score differences of {} \
             in the mapper's units)",
            m.score_denom
        );
        args.score_prob_denom = Some(d);
    }
    if args.model_truncation && !m.ends_extended {
        warn!(
            "--model-truncation: the mapper did not extend alignment spans to the read ends \
             (for spearfish, map with --extend-ends true); positions are seed spans"
        );
    }
    if args.write_assignment_probs.is_some() {
        anyhow::bail!("RAD input carries no read names; --write-assignment-probs needs a BAM");
    }
    if args.single_cell {
        anyhow::bail!("RAD input is bulk; single-cell mode needs a collated BAM");
    }
    Ok(h)
}

/// A SAM header over the RAD file's targets, with the mapper as the program.
pub fn sam_header(h: &RadHeader) -> anyhow::Result<noodles_sam::Header> {
    let mut b = noodles_sam::Header::builder();
    for (name, &len) in h.ref_names.iter().zip(&h.ref_lengths) {
        b = b.add_reference_sequence(
            name.clone(),
            HeaderMap::<header_val::map::ReferenceSequence>::new(NonZeroUsize::try_from(
                len.max(1) as usize,
            )?),
        );
    }
    let prog = if h.meta.mapper.is_empty() {
        "unknown"
    } else {
        h.meta.mapper.as_str()
    };
    b = b.add_program(prog, HeaderMap::<header_val::map::Program>::default());
    Ok(b.build())
}

/// An alignment from a RAD record, as the filter sees it.
#[derive(Debug, Clone)]
struct RadAln {
    a: Aln,
    read_len: u32,
}

impl AlnRecordLike for RadAln {
    fn opt_sequence_len(&self) -> Option<usize> {
        (self.read_len > 0).then_some(self.read_len as usize)
    }
    fn is_reverse_complemented(&self) -> bool {
        !self.a.fw
    }
    fn is_unmapped(&self) -> bool {
        false
    }
    fn ref_id(&self, _header: &noodles_sam::Header) -> anyhow::Result<usize> {
        Ok(self.a.tid as usize)
    }
    fn aln_span(&self) -> Option<usize> {
        Some(self.a.ref_end.saturating_sub(self.a.ref_start) as usize)
    }
    fn aln_score(&self) -> Option<i64> {
        Some((self.a.score * RAD_SCORE_SCALE).round() as i64)
    }
    fn aln_start(&self) -> u32 {
        self.a.ref_start
    }
    fn aln_end(&self) -> u32 {
        self.a.ref_end
    }
    fn is_supp(&self) -> bool {
        false
    }
    fn name(&self) -> Option<String> {
        None
    }
}

/// Filtered alignments of a batch of reads: alignments, their probabilities, and the read
/// boundaries.
type Batch = (Vec<AlnInfo>, Vec<f32>, Vec<usize>);

/// Parse, filter and store the RAD file's alignments, then quantify.
#[allow(clippy::too_many_arguments)]
pub fn quantify(
    path: &Path,
    header: &noodles_sam::Header,
    rad: &RadHeader,
    filter_opts: AlignmentFilters,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: NamedDigestVec,
) -> anyhow::Result<()> {
    for (i, (t, &len)) in txps.iter().zip(&rad.ref_lengths).enumerate() {
        if t.len.get() != len.max(1) as usize {
            anyhow::bail!(
                "target {i}: length {} in the header but {len} in ref_lengths",
                t.len
            );
        }
    }
    let start = std::time::Instant::now();
    let workers = NonZeroUsize::new(args.threads.saturating_sub(2).max(1)).unwrap();
    let (_, mut reader) = long_read_rad::open(path, workers)?;
    let view: Vec<TranscriptInfo> = txps.to_vec();
    let (tx, rx): (Sender<Batch>, Receiver<Batch>) = bounded(args.threads * 16);
    let discards = Mutex::new(Vec::new());
    let mut store = std::thread::scope(|s| -> anyhow::Result<InMemoryAlignmentStore> {
        let store_opts = filter_opts.clone();
        let txps_mut = &mut *txps;
        let consumer = s.spawn(move || {
            let mut store = InMemoryAlignmentStore::new(store_opts, header);
            for (alns, probs, bounds) in rx {
                for w in bounds.windows(2) {
                    let (ag, ap) = (&alns[w[0]..w[1]], &probs[w[0]..w[1]]);
                    if store.add_filtered_group(ag, ap, txps_mut) && ag.len() == 1 {
                        store.inc_unique_alignments();
                    }
                }
            }
            store
        });
        let view = &view;
        let discards = &discards;
        reader.process_parallel(workers, |mc| {
            let mut filter = filter_opts.clone();
            let mut dt = DiscardTable::new();
            let mut batch: Batch = (Vec::new(), Vec::new(), vec![0]);
            let mut group: Vec<RadAln> = Vec::new();
            for chunk in mc.iter() {
                for rec in &chunk.reads {
                    group.clear();
                    group.extend(rec.alns().map(|a| RadAln {
                        a,
                        read_len: rec.read_len,
                    }));
                    let (ag, ap) = filter.filter(&mut dt, header, view, &mut group);
                    if !ag.is_empty() {
                        batch.0.extend_from_slice(&ag);
                        batch.1.extend_from_slice(&ap);
                        batch.2.push(batch.0.len());
                    }
                }
            }
            if batch.2.len() > 1 {
                tx.send(batch).expect("RAD alignment consumer stopped");
            }
            discards.lock().unwrap().push(dt);
        })?;
        drop(tx);
        Ok(consumer.join().expect("RAD alignment consumer panicked"))
    })?;
    for dt in discards.into_inner().unwrap() {
        store.aggregate_discard_table(&dt);
    }
    let elapsed = start.elapsed();
    info!(
        "{} of {} mapped reads in the RAD file kept alignments after filtering ({})",
        store.len().to_formatted_string(&Locale::en),
        rad.counts.mapped.to_formatted_string(&Locale::en),
        humantime::format_duration(elapsed)
    );
    crate::bulk::perform_inference_and_write_output(
        header,
        &mut store,
        None,
        txps,
        txps_name,
        seqcol_digest,
        elapsed,
        args,
    )
}
