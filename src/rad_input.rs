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
    let index_info = long_read_rad::read_index_info(path)?;
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
    let mut settings = crate::prog_opts::SpearfishSettings::from_rad(m);
    let genes = index_genes(&h, index_info.as_ref());
    match &index_info {
        Some(i) => info!(
            "spearfish index: {} transcripts, {} genes (gene labels: {})",
            i.num_transcripts.to_formatted_string(&Locale::en),
            i.num_genes.to_formatted_string(&Locale::en),
            genes.source
        ),
        None => info!(
            "the RAD file does not record the index's gene labels (written before spearfish \
             0.3.1); gene labels: unknown"
        ),
    }
    // Also in the run's warnings at the end, and in meta_info.json.
    if let Some(w) = crate::util::run_warnings::check_index_genes(&genes) {
        warn!("{}", w.message);
    }
    settings.index_genes = Some(genes);
    info!("{}", settings.span_extension().describe());
    if m.quant_scores {
        info!("alignments are weighed by the mapper's quantification scores (qscore)");
    }
    args.spearfish_settings = Some(settings);
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

/// The index's gene labels as the RAD file records them: the file tags `index_genes_from`,
/// `index_num_genes`, `index_num_transcripts` and `index_gene_groups`, which spearfish 0.3.1
/// and later write (read as `long_read_rad::IndexInfo`).
///
/// `source` is `index_genes_from`, except that `given` (genes passed in through spearfish's
/// Python API) is recorded as `provided`, and an empty value (an index without a manifest) as
/// `unknown`; the counts and `has_gene_groups` are recorded in either case. A file without the
/// tags (spearfish 0.3.0) gives `source: unknown`, no gene count and no `has_gene_groups`, with
/// the transcript count of the header, and so no warning.
pub fn index_genes(
    h: &RadHeader,
    info: Option<&long_read_rad::IndexInfo>,
) -> crate::prog_opts::IndexGenes {
    let Some(info) = info else {
        return crate::prog_opts::IndexGenes::unknown(Some(h.ref_names.len()));
    };
    let source = match info.genes_from.as_str() {
        "" => "unknown",
        "given" => "provided",
        s => s,
    };
    crate::prog_opts::IndexGenes {
        source: source.to_string(),
        genes: Some(info.num_genes as usize),
        transcripts: Some(info.num_transcripts as usize),
        has_gene_groups: Some(info.gene_groups),
    }
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
    /// The quantification score (`qscore`), which is the score in files without it.
    fn aln_score(&self) -> Option<i64> {
        Some((self.a.quant_score * RAD_SCORE_SCALE).round() as i64)
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
        seqcol_digest.into(),
        elapsed,
        args,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use long_read_rad::{ChunkCodec, Counts, IndexInfo, Meta, RadWriter};

    /// A RAD file over three targets, with the index description `info` when given (as
    /// spearfish 0.3.1 writes it) and without (as spearfish 0.3.0 did).
    fn write_rad(path: &Path, info: Option<&IndexInfo>) {
        let refs = [("tA", 2000u32), ("tB", 1500), ("tC", 3000)];
        let meta = Meta {
            mapper: "spearfish".into(),
            mapper_version: "0.3.1".into(),
            ..Default::default()
        };
        let w = match info {
            Some(i) => RadWriter::create_with_index(path, &refs, &meta, i, ChunkCodec::None),
            None => RadWriter::create(path, &refs, &meta, ChunkCodec::None),
        }
        .unwrap();
        w.finish(Counts::default()).unwrap();
    }

    fn settings_of(path: &Path) -> crate::prog_opts::SpearfishSettings {
        let p = path.to_str().unwrap();
        let mut args = Args::try_parse_from(["oarfish", "-a", p, "-o", "o"]).unwrap();
        read_header(path, &mut args).unwrap();
        args.spearfish_settings.unwrap()
    }

    /// The index's gene labels are read from the RAD file tags: a gene-less index is recorded
    /// and warned about, a grouped one is recorded and not; `given` reads as `provided` and an
    /// empty source as `unknown`; a file without the tags is `unknown`, without a warning.
    #[test]
    fn rad_index_genes_from_file_tags() {
        use crate::util::run_warnings::check_index_genes;
        let dir = tempfile::tempdir().unwrap();
        let genes_of = |name: &str, info: Option<IndexInfo>| {
            let path = dir.path().join(name);
            write_rad(&path, info.as_ref());
            settings_of(&path).index_genes.unwrap()
        };

        let none = genes_of("none.rad", Some(IndexInfo::new("none", 3, 3, false)));
        assert_eq!(none.source, "none");
        assert_eq!((none.genes, none.transcripts), (Some(3), Some(3)));
        assert_eq!(none.has_gene_groups, Some(false));
        let w = check_index_genes(&none).unwrap();
        assert_eq!(w.code, "index_without_gene_groups");
        assert!(w.message.contains("carry no gene IDs") && w.message.contains("--t2g"));

        let fasta = genes_of("fasta.rad", Some(IndexInfo::new("fasta", 2, 3, true)));
        assert_eq!(fasta.source, "fasta");
        assert_eq!((fasta.genes, fasta.transcripts), (Some(2), Some(3)));
        assert_eq!(fasta.has_gene_groups, Some(true));
        assert!(check_index_genes(&fasta).is_none());

        let given = genes_of("given.rad", Some(IndexInfo::new("given", 3, 3, false)));
        assert_eq!(given.source, "provided");
        let w = check_index_genes(&given).unwrap();
        assert!(
            w.message
                .contains("gene labels given when the index was built")
        );

        let empty = genes_of("empty.rad", Some(IndexInfo::new("", 2, 3, true)));
        assert_eq!(empty.source, "unknown");
        assert_eq!((empty.genes, empty.has_gene_groups), (Some(2), Some(true)));

        let old = genes_of("old.rad", None);
        assert_eq!(old, crate::prog_opts::IndexGenes::unknown(Some(3)));
        assert!(check_index_genes(&old).is_none());
    }
}
