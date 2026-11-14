use crate::prog_opts::{Args, SequencingTech};
use crate::util::digest_utils;
use crate::util::file_utils::{get_ref_source, is_fasta};
use crate::util::oarfish_types::NamedDigestVec;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;

use num_format::{Locale, ToFormattedString};

use std::fs::File;
use std::num::NonZeroUsize;

use tracing::{info, warn};

pub(crate) type HeaderReaderAlignerDigest = (
    noodles_sam::header::Header,
    Option<bam::io::Reader<bgzf::io::MultithreadedReader<File>>>,
    Option<crate::util::mapper::Mapper>,
    NamedDigestVec,
);

// ===========================================================================
// rammap backend
// ===========================================================================
//
// Builds a pure-Rust `rammap` aligner. The FASTA (`--annotated`/`--novel`,
// `--genome`) construction path and prebuilt-index loading (rammap RMMI or
// legacy `.mmi`) are both supported.

/// Per-base match score of the rammap transcriptome preset used for `seq_tech` (rammap follows
/// minimap2: 1 for `map-hifi`, 2 for `map-ont` and `map-pb`).
pub fn rammap_match_score(seq_tech: &Option<SequencingTech>) -> Option<i32> {
    match rammap_txome_preset(seq_tech).ok()? {
        rammap::api::Preset::MapHifi => Some(1),
        rammap::api::Preset::MapOnt | rammap::api::Preset::MapPb => Some(2),
        _ => None,
    }
}

/// Map an oarfish read-mode sequencing tech onto a rammap transcriptome preset.
fn rammap_txome_preset(seq_tech: &Option<SequencingTech>) -> anyhow::Result<rammap::api::Preset> {
    match seq_tech {
        Some(SequencingTech::OntCDNA) | Some(SequencingTech::OntDRNA) => {
            Ok(rammap::api::Preset::MapOnt)
        }
        Some(SequencingTech::PacBio) => Ok(rammap::api::Preset::MapPb),
        Some(SequencingTech::PacBioHifi) => Ok(rammap::api::Preset::MapHifi),
        None => anyhow::bail!("sequencing tech must be provided in read mode, but it was not!"),
    }
}

/// Hint rammap's (rayon) index-build parallelism to use `n` threads, unless the
/// user has already pinned `RAYON_NUM_THREADS`. rammap has no library-level
/// thread setter, so the global rayon pool is the only lever; it is initialized
/// lazily on first use (the index build), so setting the env var here takes
/// effect as long as we do it before constructing the aligner.
fn rammap_set_index_threads(n: usize) {
    if std::env::var_os("RAYON_NUM_THREADS").is_none() {
        // SAFETY: called from the single-threaded setup path before any rayon
        // work (and thus any other thread that might read the environment) runs.
        unsafe { std::env::set_var("RAYON_NUM_THREADS", n.max(1).to_string()) };
    }
}

/// Build a SAM header from a built rammap index, emitting one `@SQ` record per
/// reference sequence in `target_id` order (so a mapping's `target_id` equals
/// its header reference id).
fn rammap_make_header(
    aligner: &rammap::api::Aligner,
) -> anyhow::Result<noodles_sam::header::Header> {
    let mut header = noodles_sam::header::Header::builder();
    for ts in aligner.index().seqs.iter() {
        header = header.add_reference_sequence(
            ts.name.clone(),
            HeaderMap::<header_val::map::ReferenceSequence>::new(NonZeroUsize::try_from(ts.len)?),
        );
    }
    header = header.add_program("rammap", HeaderMap::<header_val::map::Program>::default());
    Ok(header.build())
}

/// Load a pre-built index for the rammap backend. rammap's `from_index` accepts
/// both rammap's native RMMI format and legacy `.mmi` files (it detects the
/// magic and parses the binary layout), so an index built by an earlier oarfish
/// release (`--index-out`) can be reused here. Indexing-related `--seq-tech`
/// parameters are ignored: the index's own k/w are used. The oarfish digest
/// footer (if present) is read back; otherwise a name/length-based digest is
/// recorded.
fn get_aligner_from_index(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    let idx_file = args.index.clone().expect("index file should exist");
    let idx_str = idx_file
        .to_str()
        .expect("could not convert index path to &str");
    if is_spearfish_index(&idx_file) {
        anyhow::bail!(
            "{} is a spearfish index; pass --mapper spearfish to use it",
            idx_file.display()
        );
    }

    warn!(
        "You are using an existing index with the rammap backend (rammap reads both its own RMMI \
         format and legacy .mmi files). Indexing-related parameters implied by --seq-tech are \
         ignored; the index's own k/w are used."
    );

    rammap_set_index_threads(args.threads.max(1));
    let preset = rammap_txome_preset(&args.seq_tech)?;
    let mut aligner = rammap::api::Aligner::from_index(idx_str, preset)?;
    aligner.options_mut().filtering.best_n = args.best_n as i32;
    aligner.options_mut().filtering.seed = 11;

    let n_seq = aligner.index().seqs.len();
    info!(
        "index contains {} sequences",
        n_seq.to_formatted_string(&Locale::en)
    );

    let header = rammap_make_header(&aligner)?;

    // Prefer the oarfish digest footer appended at index-build time. If absent,
    // recompute the full reference signature from the index sequences (rammap
    // retains them); only if the index is sequence-stripped do we fall
    // back to a name/length-only digest.
    let digests = match digest_utils::read_digest_footer(idx_str) {
        Ok(d) => d,
        _ => {
            if aligner.index().has_sequences() {
                warn!(
                    "the provided index has no oarfish digest footer; recomputing the \
                     reference signature from the index sequences."
                );
                let d = digest_utils::digest_from_rammap_index(&aligner)?;
                vec![("prexisting_index".to_string(), d)].into()
            } else {
                warn!(
                    "the provided index has no oarfish digest footer and no sequences; \
                     recording a sequence-name/length digest instead of the full \
                     reference signature."
                );
                let d = digest_utils::digest_from_header(&header)?;
                vec![("prexisting_index".to_string(), d)].into()
            }
        }
    };

    Ok((
        header,
        None,
        Some(crate::util::mapper::Mapper::Rammap(std::sync::Arc::new(
            aligner,
        ))),
        digests,
    ))
}

// ===========================================================================
// spearfish backend
// ===========================================================================

/// Map an oarfish read-mode sequencing tech onto a spearfish preset.
fn spearfish_preset(seq_tech: &Option<SequencingTech>) -> anyhow::Result<spearfish::Preset> {
    match seq_tech {
        Some(SequencingTech::OntCDNA) => Ok(spearfish::Preset::Ont),
        Some(SequencingTech::OntDRNA) => Ok(spearfish::Preset::Drna),
        Some(SequencingTech::PacBioHifi) => Ok(spearfish::Preset::Hifi),
        Some(SequencingTech::PacBio) => {
            warn!("spearfish has no PacBio CLR preset; using its ONT preset.");
            Ok(spearfish::Preset::Ont)
        }
        None => anyhow::bail!("sequencing tech must be provided in read mode, but it was not!"),
    }
}

/// Whether `path` is a spearfish index (`.sfi`): its first bytes are spearfish's magic.
pub(crate) fn is_spearfish_index(path: &std::path::Path) -> bool {
    use std::io::Read;
    let mut magic = [0u8; 8];
    File::open(path)
        .and_then(|mut f| f.read_exact(&mut magic))
        .is_ok_and(|_| &magic == b"SPFISH01")
}

/// Build (from `--annotated`/`--novel`) or open (`--index`, a spearfish `.sfi`) a spearfish
/// index and mapper. The header lists transcripts in index order (= FASTA order), with
/// spearfish's transcript names (FASTA IDs cut at the first `|`). With `--index-out`, the
/// built index is written there (with the reference digests in oarfish's footer) and mapped
/// from the file.
fn get_spearfish_mapper(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    info!("oarfish is operating in read-based mode (spearfish backend)");
    let mut digests = NamedDigestVec::new();
    let index = if let Some(idx) = &args.index {
        let index = spearfish::Index::open(idx)?;
        // Prefer the digests oarfish wrote into the footer when it built the index; an index
        // from `spearfish index` has none, so fall back to a name/length digest.
        let idx_str = idx
            .to_str()
            .expect("could not convert --index path to &str");
        match digest_utils::read_digest_footer(idx_str) {
            Ok(d) => digests = d,
            Err(_) => {
                warn!(
                    "the provided spearfish index has no oarfish digest footer; recording a \
                     sequence-name/length digest instead of the full reference signature."
                );
                digests.push((
                    "prexisting_index".to_string(),
                    digest_utils::digest_from_header(&spearfish_header(&index)?)?,
                ));
            }
        }
        index
    } else {
        rammap_set_index_threads(args.threads.max(1));
        // Reference digests hash the FASTAs on background threads, overlapping with
        // reading and indexing (as on the rammap path); joined once the index is built.
        let digest_handles: Vec<_> = [
            ("annotated_transcripts_digest", &args.annotated),
            ("novel_transcripts_digest", &args.novel),
        ]
        .into_iter()
        .filter_map(|(label, f)| {
            f.as_ref()
                .map(|f| (label, digest_utils::get_digest_from_fasta(f)))
        })
        .collect();
        let t_read = std::time::Instant::now();
        let mut txs = Vec::new();
        for f in [&args.annotated, &args.novel].into_iter().flatten() {
            let mut reader = needletail::parse_fastx_file(f)?;
            while let Some(rec) = reader.next() {
                let rec = rec?;
                let (name, gene) = spearfish::name_and_gene(&String::from_utf8_lossy(rec.id()));
                txs.push((name, gene, rec.seq().into_owned()));
            }
        }
        info!(
            "read {} transcripts in {:.1}s",
            txs.len(),
            t_read.elapsed().as_secs_f64()
        );
        let t = std::time::Instant::now();
        let idx_opts = spearfish::IndexOptions::default();
        let (index, rep) = if let Some(out) = &args.index_out {
            let input = spearfish::provenance::InputInfo {
                path: [&args.annotated, &args.novel]
                    .into_iter()
                    .flatten()
                    .map(|f| spearfish::provenance::absolute(f))
                    .collect::<Vec<_>>()
                    .join(","),
                genes_from: if txs.iter().any(|(name, gene, _)| name != gene) {
                    "fasta"
                } else {
                    "none"
                }
                .to_string(),
                ..Default::default()
            };
            let (_, rep, _) = spearfish::Index::build_to_file_with(
                txs,
                &idx_opts,
                out,
                input,
                std::env::args().collect(),
            )?;
            (None, rep)
        } else {
            let (index, rep) = spearfish::Index::build(txs, &idx_opts)?;
            (Some(index), rep)
        };
        info!(
            "built spearfish index over {} transcripts ({} genes) in {:.1}s ({})",
            rep.transcripts.to_formatted_string(&Locale::en),
            rep.genes.to_formatted_string(&Locale::en),
            t.elapsed().as_secs_f64(),
            rep.seconds
                .iter()
                .map(|(k, v)| format!("{k} {v:.1}s"))
                .collect::<Vec<_>>()
                .join(", ")
        );
        for (label, h) in digest_handles {
            digests.push((label.to_string(), h.join().expect("valid digest")?));
        }
        match index {
            Some(index) => index,
            None => {
                let out = args
                    .index_out
                    .as_ref()
                    .expect("index written to --index-out");
                let out_str = out
                    .to_str()
                    .expect("could not convert --index-out path to &str");
                digest_utils::append_digest_footer(out_str, &digests)?;
                info!("wrote spearfish index to {}", out.display());
                spearfish::Index::open(out)?
            }
        }
    };
    let header = spearfish_header(&index)?;
    let opts = spearfish_map_options(args)?;
    let settings = spearfish_settings(args, &opts)?;
    info!("spearfish: {}", settings.span_extension().describe());
    args.spearfish_settings = Some(settings);
    let mapper = spearfish::Mapper::new(std::sync::Arc::new(index), &opts)?;
    Ok((
        header,
        None,
        Some(crate::util::mapper::Mapper::Spearfish(std::sync::Arc::new(
            mapper,
        ))),
        digests,
    ))
}

/// The spearfish mapping options: the preset `--seq-tech` names, with oarfish's changes and
/// the spearfish options given.
fn spearfish_map_options(args: &Args) -> anyhow::Result<spearfish::MapOptions> {
    let preset = spearfish_preset(&args.seq_tech)?;
    let mut opts = spearfish::MapOptions::preset(preset);
    // The truncation model needs where reads start and end on each transcript, so spans are
    // extended to the read ends even where the preset leaves them at the outermost seeds.
    // This also has the end fit measure every candidate by alignment (drna's chain measurement
    // only saves alignments that spans left unextended would otherwise run for it alone),
    // unless `--end-fit-from chain` asks for the chain.
    if args.model_truncation {
        opts.extend_ends(preset);
    }
    // The end fit (the quantification score): the preset's settings unless given.
    if let Some(tb) = &mut opts.params.tiebreak {
        if let Some(w) = args.end_fit {
            tb.end_fit = w;
        }
        if let Some(t) = args.end_fit_tol {
            tb.end_fit_tol = t;
        }
        if let Some(from) = args.end_fit_from {
            tb.end_fit_chain = from == crate::prog_opts::EndFitFrom::Chain;
        }
    }
    // The orientation evidence for antisense pairs: the preset's unless given (or unless
    // `ORIENT_EVIDENCE_DEFAULT` names a default of oarfish's own).
    if let (Some(ev), Some(tb)) = (
        args.orient_evidence
            .or(crate::prog_opts::ORIENT_EVIDENCE_DEFAULT),
        &mut opts.params.tiebreak,
    ) {
        use crate::prog_opts::{ORIENT_TAIL_AFTER_OFF, OrientEvidence};
        match ev {
            OrientEvidence::Quant | OrientEvidence::Mapping => {
                tb.orient_quant = ev == OrientEvidence::Quant;
            }
            // As spearfish's `--no-orient-same-exact --orient-tail-after=80`.
            OrientEvidence::Off => {
                tb.orient_same_exact = false;
                tb.orient_tail_after = ORIENT_TAIL_AFTER_OFF;
            }
        }
    }
    // Only hits within the score filter of the best usually reach the EM, so only those are
    // extended. (Reporting only those lost reads whose best hit fails the span filters: oarfish
    // then rescales the threshold to the best hit that passes.)
    if let Some(e) = &mut opts.params.ends {
        e.within = 1.0 - args.score_threshold.try_as_f32().unwrap_or(0.95);
    }
    // Reads known to be oriented need one orientation mapped (as the drna preset maps).
    if let Some(strand) = args.map_strand {
        use bio_types::strand::Strand;
        opts.params.strands = match strand {
            Strand::Forward => spearfish::Strands::Forward,
            Strand::Reverse => spearfish::Strands::Reverse,
            Strand::Unknown => spearfish::Strands::Both,
        };
    }
    Ok(opts)
}

/// The spearfish settings recorded in `meta_info.json`: the preset, span extension and why,
/// the orientations mapped and the quantification score's settings, as resolved from the
/// preset and oarfish's options.
fn spearfish_settings(
    args: &Args,
    opts: &spearfish::MapOptions,
) -> anyhow::Result<crate::prog_opts::SpearfishSettings> {
    use crate::prog_opts::{EndFitFrom, SpanExtension, SpearfishSettings};
    let preset = spearfish_preset(&args.seq_tech)?;
    let span = SpanExtension::in_process(
        opts.preset.name().to_string(),
        spearfish::MapOptions::preset(preset).params.ends.is_some(),
        args.model_truncation,
    );
    debug_assert_eq!(span.extended, opts.params.ends.is_some());
    let tb = opts.params.tiebreak.as_ref();
    Ok(SpearfishSettings {
        source: "in_process",
        mapper: None,
        mapper_version: None,
        preset: span.preset,
        extend_ends: opts.params.ends.is_some(),
        extend_ends_reason: span.reason,
        strands: opts.params.strands.name().to_string(),
        quant_scores: opts.quant_scores(),
        end_fit: Some(tb.map_or(0.0, |tb| tb.end_fit)),
        end_fit_tol: tb.map(|tb| tb.end_fit_tol),
        end_fit_from: tb.map(|tb| {
            if tb.end_fit_chain {
                EndFitFrom::Chain
            } else {
                EndFitFrom::Alignment
            }
        }),
        orient_quant: tb.map(|tb| tb.orient_quant),
        orient_evidence: tb.map(|tb| {
            crate::prog_opts::OrientEvidence::of(
                tb.orient_same_exact,
                tb.orient_tail_after,
                tb.orient_quant,
            )
        }),
    })
}

fn spearfish_header(index: &spearfish::Index) -> anyhow::Result<noodles_sam::header::Header> {
    let mut header = noodles_sam::header::Header::builder();
    for t in 0..index.n_transcripts() as u32 {
        header = header.add_reference_sequence(
            index.tx_name(t).to_string(),
            HeaderMap::<header_val::map::ReferenceSequence>::new(NonZeroUsize::try_from(
                index.tx_len(t).max(1) as usize,
            )?),
        );
    }
    header = header.add_program(
        "spearfish",
        HeaderMap::<header_val::map::Program>::default(),
    );
    Ok(header.build())
}

pub(crate) fn get_aligner_from_args(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    if args.mapper == crate::prog_opts::MapperBackend::Spearfish {
        return get_spearfish_mapper(args);
    }
    info!("oarfish is operating in read-based mode (rammap backend)");
    if args.index.is_some() {
        return get_aligner_from_index(args);
    }

    assert!(
        args.annotated
            .as_ref()
            .is_none_or(|f| is_fasta(f).expect("couldn't read input file to `--annotated`."))
    );
    assert!(
        args.novel
            .as_ref()
            .is_none_or(|f| is_fasta(f).expect("couldn't read input file to `--novel`."))
    );

    // compute the reference digests on background threads.
    let ref_digest_handle = args
        .annotated
        .clone()
        .map(|refs| digest_utils::get_digest_from_fasta(&refs));
    let novel_digest_handle = args
        .novel
        .clone()
        .map(|refs| digest_utils::get_digest_from_fasta(&refs));

    let thread_sub = if ref_digest_handle
        .as_ref()
        .xor(novel_digest_handle.as_ref())
        .is_some()
    {
        1
    } else {
        2
    };
    let idx_threads = args.threads.saturating_sub(thread_sub).max(1);
    rammap_set_index_threads(idx_threads);

    // combine annotated + novel references (via a fifo) when both are present.
    let input_source = get_ref_source(args.annotated.clone(), args.novel.clone())?;
    let input_path = input_source
        .file_path()
        .to_str()
        .expect("could not convert reference path to &str")
        .to_owned();

    let preset = rammap_txome_preset(&args.seq_tech)?;
    let mut aligner = rammap::api::Aligner::from_fasta(&input_path, preset)?;
    // up to best_n secondary mappings; fixed mapping seed.
    aligner.options_mut().filtering.best_n = args.best_n as i32;
    aligner.options_mut().filtering.seed = 11;

    let n_seq = aligner.index().seqs.len();
    info!(
        "index contains {} sequences",
        n_seq.to_formatted_string(&Locale::en)
    );

    input_source.join_if_needed()?;

    let header = rammap_make_header(&aligner)?;

    let mut digests = NamedDigestVec::new();
    if let Some(handle) = ref_digest_handle {
        let digest = handle.join().expect("valid digest")?;
        digests.push(("annotated_transcripts_digest".to_string(), digest));
    }
    if let Some(handle) = novel_digest_handle {
        let digest = handle.join().expect("valid digest")?;
        digests.push(("novel_transcripts_digest".to_string(), digest));
    }

    // if the user requested a persistent index alongside quantification, write it
    // (with its reference-signature footer) now that the in-memory index is built.
    if let Some(idx_out) = args.index_out.clone() {
        let idx_str = idx_out
            .to_str()
            .expect("could not convert --index-out path to &str");
        write_rammap_index_with_footer(&aligner, idx_str, &digests)?;
    }

    Ok((
        header,
        None,
        Some(crate::util::mapper::Mapper::Rammap(std::sync::Arc::new(
            aligner,
        ))),
        digests,
    ))
}

/// Build a *spliced* rammap aligner over a genome (FASTA) for genome-read mode,
/// returning the reference (chromosome) names in `target_id` order. The splice
/// preset is chosen from the sequencing technology (`splice:hq` for PacBio HiFi,
/// `splice` otherwise).
///
/// Annotated splice junctions (when `junc_bed` is provided) are loaded via
/// rammap's `load_junctions_bed` (the `--junc-bed` equivalent).
pub(crate) fn get_genome_aligner_from_args(
    args: &Args,
    junc_bed: Option<&std::path::Path>,
) -> anyhow::Result<(crate::util::mapper::RammapMapper, Vec<String>)> {
    let genome = args
        .genome
        .clone()
        .expect("--genome is required in genome read mode");
    info!(
        "oarfish is operating in genome read mode (rammap backend); building a spliced index over {}",
        genome.display()
    );

    rammap_set_index_threads(args.threads.max(1));

    let preset = match args.seq_tech {
        Some(SequencingTech::PacBioHifi) => rammap::api::Preset::SpliceHq,
        Some(_) => rammap::api::Preset::Splice,
        None => {
            anyhow::bail!("--seq-tech is required in genome read mode, but it was not provided!")
        }
    };

    let genome_str = genome
        .to_str()
        .expect("could not convert genome path to &str");
    // `--genome` may be a FASTA (build the index now) or a prebuilt rammap index
    // (load it). A loaded index that retains its reference sequences also enables
    // soft-clip rescue without a separate FASTA (see `rescue_db_from_genome_index`).
    let mut aligner = if is_fasta(&genome).unwrap_or(false) {
        rammap::api::Aligner::from_fasta(genome_str, preset)?
    } else {
        info!(
            "loading prebuilt rammap genome index from {}",
            genome.display()
        );
        rammap::api::Aligner::from_index(genome_str, preset)?
    };
    aligner.options_mut().filtering.best_n = args.best_n as i32;
    aligner.options_mut().filtering.seed = 11;

    // Load annotated splice junctions (the `--junc-bed` equivalent) if provided.
    if let Some(bed) = junc_bed {
        let bs = bed.to_str().expect("junction BED path must be valid UTF-8");
        aligner
            .load_junctions_bed(bs)
            .map_err(|e| anyhow::anyhow!("failed to load splice junctions {}: {}", bs, e))?;
        info!("loaded splice junctions into the genome index from {}", bs);
    }

    let refnames: Vec<String> = aligner
        .index()
        .seqs
        .iter()
        .map(|t| t.name.clone())
        .collect();
    info!(
        "genome index contains {} reference sequences",
        refnames.len().to_formatted_string(&Locale::en)
    );
    Ok((std::sync::Arc::new(aligner), refnames))
}

/// Build a rammap index from the provided reference FASTA(s) and (when
/// `--index-out` is set) persist it to disk as a native RMMI file with the
/// oarfish reference-signature footer appended. This backs the `--only-index`
/// path; the resulting signature is computed from the FASTA via
/// [`digest_utils::get_digest_from_fasta`].
pub(crate) fn get_aligner_from_fastas(
    args: &mut Args,
) -> anyhow::Result<HeaderReaderAlignerDigest> {
    assert!(
        args.annotated
            .as_ref()
            .is_none_or(|f| is_fasta(f).expect("couldn't read input file to `--annotated`."))
    );
    assert!(
        args.novel
            .as_ref()
            .is_none_or(|f| is_fasta(f).expect("couldn't read input file to `--novel`."))
    );

    // compute the reference digests on background threads.
    let ref_digest_handle = args
        .annotated
        .clone()
        .map(|refs| digest_utils::get_digest_from_fasta(&refs));
    let novel_digest_handle = args
        .novel
        .clone()
        .map(|refs| digest_utils::get_digest_from_fasta(&refs));

    let thread_sub = if ref_digest_handle
        .as_ref()
        .xor(novel_digest_handle.as_ref())
        .is_some()
    {
        1
    } else {
        2
    };
    let idx_threads = args.threads.saturating_sub(thread_sub).max(1);
    rammap_set_index_threads(idx_threads);

    // combine annotated + novel references (via a fifo) when both are present.
    let input_source = get_ref_source(args.annotated.clone(), args.novel.clone())?;
    let input_path = input_source
        .file_path()
        .to_str()
        .expect("could not convert reference path to &str")
        .to_owned();

    let preset = rammap_txome_preset(&args.seq_tech)?;
    let mut aligner = rammap::api::Aligner::from_fasta(&input_path, preset)?;
    aligner.options_mut().filtering.best_n = args.best_n as i32;
    aligner.options_mut().filtering.seed = 11;

    let n_seq = aligner.index().seqs.len();
    info!(
        "index contains {} sequences",
        n_seq.to_formatted_string(&Locale::en)
    );

    input_source.join_if_needed()?;

    let header = rammap_make_header(&aligner)?;

    let mut digests = NamedDigestVec::new();
    if let Some(handle) = ref_digest_handle {
        let digest = handle.join().expect("valid digest")?;
        digests.push(("annotated_transcripts_digest".to_string(), digest));
    }
    if let Some(handle) = novel_digest_handle {
        let digest = handle.join().expect("valid digest")?;
        digests.push(("novel_transcripts_digest".to_string(), digest));
    }

    // persist the index (+ signature footer) if the user asked for it.
    if let Some(idx_out) = args.index_out.clone() {
        let idx_str = idx_out
            .to_str()
            .expect("could not convert --index-out path to &str");
        write_rammap_index_with_footer(&aligner, idx_str, &digests)?;
    }

    Ok((
        header,
        None,
        Some(crate::util::mapper::Mapper::Rammap(std::sync::Arc::new(
            aligner,
        ))),
        digests,
    ))
}

/// Persist a rammap aligner's index to `idx_out` as a native RMMI file and
/// append the oarfish reference-signature footer. The footer is tolerated by
/// rammap's loader (which stops at the end of the bincode structure), so the
/// index remains a single self-contained file that round-trips through
/// [`rammap::api::Aligner::from_index`] while carrying its signature.
fn write_rammap_index_with_footer(
    aligner: &rammap::api::Aligner,
    idx_out: &str,
    digests: &NamedDigestVec,
) -> anyhow::Result<()> {
    aligner
        .save_index(idx_out)
        .map_err(|e| anyhow::anyhow!("failed to write rammap index to {}: {}", idx_out, e))?;
    digest_utils::append_digest_footer(idx_out, digests)?;
    info!(
        "wrote rammap index to {} with oarfish reference-signature footer",
        idx_out
    );
    Ok(())
}

/// Build a bramble rescue reference (`FastaDb`) from the reference sequences
/// already resident in the loaded genome aligner index, so genome-mode soft-clip
/// rescue works from a prebuilt index (no separate FASTA). Returns `None` if the
/// index stores no sequences.
pub(crate) fn rescue_db_from_genome_index(
    aligner: &rammap::api::Aligner,
) -> Option<bramble_rs::fasta::FastaDb> {
    let idx = aligner.index();
    if !idx.has_sequences() {
        info!(
            "the genome index stores no reference sequences; soft-clip rescue is off \
             (pass --genome-fasta, or rebuild the index with sequences, to enable it)"
        );
        return None;
    }
    info!("sourcing soft-clip-rescue reference sequences from the loaded genome index");
    // Build bramble's own map type so it's moved into `from_seq_map` without a
    // rehash and picks up bramble's deterministic hasher (see bramble_rs SeqMap).
    let mut seqs = bramble_rs::fasta::FastaDb::new_seq_map(idx.seqs.len());
    let mut nt4: Vec<u8> = Vec::new();
    for (rid, ts) in idx.seqs.iter().enumerate() {
        nt4.clear();
        nt4.resize(ts.len, 0u8);
        idx.extract_nt4_into(rid, 0, ts.len, &mut nt4);
        // nt4 (0..=4) -> ASCII; bramble's get_slice upper-cases/N-maps identically
        // for a file-loaded FASTA, so the result is interchangeable.
        let ascii: Vec<u8> = nt4
            .iter()
            .map(|&b| match b {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            })
            .collect();
        let name = ts.name.split_whitespace().next().unwrap_or("").to_string();
        seqs.insert(name, ascii);
    }
    Some(bramble_rs::fasta::FastaDb::from_seq_map(seqs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    fn opts(extra: &[&str]) -> spearfish::MapOptions {
        let mut v = vec!["oarfish", "--mapper", "spearfish", "--reads", "r.fq"];
        v.extend_from_slice(&["--index", "t.sfi", "-o", "out"]);
        v.extend_from_slice(extra);
        spearfish_map_options(&Args::try_parse_from(v).unwrap()).unwrap()
    }

    fn end_fit(o: &spearfish::MapOptions) -> (f32, u32, bool) {
        let tb = o.params.tiebreak.as_ref().unwrap();
        (tb.end_fit, tb.end_fit_tol, tb.end_fit_chain)
    }

    /// Without the end-fit options the options are the preset's, and each option sets its
    /// setting alone.
    #[test]
    fn end_fit_options_default_to_the_preset() {
        for (tech, p) in [
            ("ont-cdna", spearfish::Preset::Ont),
            ("ont-drna", spearfish::Preset::Drna),
            ("pac-bio-hifi", spearfish::Preset::Hifi),
        ] {
            let o = opts(&["--seq-tech", tech]);
            assert_eq!(
                end_fit(&o),
                end_fit(&spearfish::MapOptions::preset(p)),
                "{tech}"
            );
        }
        let drna = |extra: &[&str]| {
            let mut v = vec!["--seq-tech", "ont-drna"];
            v.extend_from_slice(extra);
            end_fit(&opts(&v))
        };
        let (w, t, chain) = drna(&[]);
        assert!(w > 0.0 && chain);
        assert_eq!(drna(&["--end-fit-from", "alignment"]), (w, t, false));
        assert_eq!(drna(&["--end-fit", "0"]), (0.0, t, true));
        assert_eq!(
            drna(&["--end-fit", "2.5", "--end-fit-tol", "6"]),
            (2.5, 6, true)
        );
        let ont = end_fit(&opts(&[
            "--seq-tech",
            "ont-cdna",
            "--end-fit-from",
            "chain",
        ]));
        assert!(ont.2);
    }

    /// `--model-truncation` extends spans, and the end fit then measures by alignment in
    /// every preset, unless `--end-fit-from chain` is given.
    #[test]
    fn truncation_model_measures_the_end_fit_by_alignment() {
        for tech in ["ont-cdna", "ont-drna", "pac-bio-hifi"] {
            let o = opts(&["--seq-tech", tech, "--model-truncation"]);
            assert!(o.params.ends.is_some(), "{tech}");
            assert!(!end_fit(&o).2, "{tech}");
        }
        let o = opts(&[
            "--seq-tech",
            "ont-drna",
            "--model-truncation",
            "--end-fit-from",
            "chain",
        ]);
        assert!(o.params.ends.is_some() && end_fit(&o).2);
    }

    fn settings(extra: &[&str]) -> serde_json::Value {
        let mut v = vec!["oarfish", "--mapper", "spearfish", "--reads", "r.fq"];
        v.extend_from_slice(&["--index", "t.sfi", "-o", "out"]);
        v.extend_from_slice(extra);
        let args = Args::try_parse_from(v).unwrap();
        let o = spearfish_map_options(&args).unwrap();
        serde_json::to_value(spearfish_settings(&args, &o).unwrap()).unwrap()
    }

    #[test]
    fn end_fit_options_are_recorded() {
        let s = settings(&["--seq-tech", "ont-drna", "--end-fit-tol", "3"]);
        assert_eq!(s["source"], "in_process");
        assert_eq!(s["preset"], "drna");
        assert_eq!(s["end_fit_tol"], 3);
        assert_eq!(s["end_fit_from"], "chain");
        assert_eq!(s["extend_ends"], false);
        assert_eq!(s["extend_ends_reason"], "preset_off");
        assert_eq!(s["strands"], "fwd");
        assert_eq!(s["quant_scores"], true);
        assert!(s.get("mapper").is_none());
        let s = settings(&["--seq-tech", "ont-cdna", "--end-fit", "0"]);
        assert_eq!(s["end_fit"], 0.0);
        assert_eq!(s["end_fit_from"], "alignment");
        assert_eq!(s["extend_ends_reason"], "preset");
        assert_eq!(s["strands"], "both");
        // The orientation evidence still acts through the quantification score.
        assert_eq!(s["orient_quant"], true);
        assert_eq!(s["quant_scores"], true);
    }

    /// Span extension and why are recorded with the other settings (and agree with the
    /// options the mapper is made with).
    #[test]
    fn span_extension_is_recorded_with_the_settings() {
        for (tech, extra, ext, reason) in [
            ("ont-cdna", None, true, "preset"),
            ("ont-cdna", Some("--model-truncation"), true, "preset"),
            ("ont-drna", None, false, "preset_off"),
            (
                "ont-drna",
                Some("--model-truncation"),
                true,
                "model_truncation",
            ),
            ("pac-bio-hifi", None, false, "preset_off"),
            (
                "pac-bio-hifi",
                Some("--model-truncation"),
                true,
                "model_truncation",
            ),
        ] {
            let mut v = vec!["--seq-tech", tech];
            v.extend(extra);
            let s = settings(&v);
            assert_eq!(s["extend_ends"], ext, "{tech} {extra:?}");
            assert_eq!(s["extend_ends_reason"], reason, "{tech} {extra:?}");
        }
    }

    fn orient(o: &spearfish::MapOptions) -> (bool, u32, bool) {
        let tb = o.params.tiebreak.as_ref().unwrap();
        (tb.orient_same_exact, tb.orient_tail_after, tb.orient_quant)
    }

    /// Without `--orient-evidence` the orientation evidence is the preset's (through the
    /// quantification score); `off` is spearfish's `--no-orient-same-exact
    /// --orient-tail-after=80`, which with `--end-fit 0` leaves no quantification score of its
    /// own; `mapping` is spearfish's `--no-orient-quant`.
    #[test]
    fn orient_evidence_option() {
        for (tech, p) in [
            ("ont-cdna", spearfish::Preset::Ont),
            ("ont-drna", spearfish::Preset::Drna),
            ("pac-bio-hifi", spearfish::Preset::Hifi),
        ] {
            let preset = spearfish::MapOptions::preset(p);
            assert_eq!(
                orient(&opts(&["--seq-tech", tech])),
                orient(&preset),
                "{tech}"
            );
            assert_eq!(
                settings(&["--seq-tech", tech])["orient_evidence"],
                "quant",
                "{tech}"
            );
            let off = opts(&["--seq-tech", tech, "--orient-evidence", "off"]);
            assert_eq!(orient(&off), (false, 80, orient(&preset).2), "{tech}");
            let off0 = opts(&[
                "--seq-tech",
                tech,
                "--orient-evidence",
                "off",
                "--end-fit",
                "0",
            ]);
            assert!(!off0.quant_scores(), "{tech}");
            let s = settings(&[
                "--seq-tech",
                tech,
                "--orient-evidence",
                "off",
                "--end-fit",
                "0",
            ]);
            assert_eq!(s["orient_evidence"], "off");
            assert_eq!(s["quant_scores"], false);
            let m = opts(&["--seq-tech", tech, "--orient-evidence", "mapping"]);
            let (exact, after, quant) = orient(&m);
            assert!(exact && after == orient(&preset).1 && !quant, "{tech}");
            assert_eq!(
                settings(&["--seq-tech", tech, "--orient-evidence", "mapping"])["orient_evidence"],
                "mapping"
            );
        }
        // `off` is what spearfish's own presets without the evidence (the `-fast` ones) use.
        let fast = spearfish::MapOptions::preset(spearfish::Preset::OntFast);
        let (exact, after, _) = orient(&fast);
        assert!(!exact);
        assert_eq!(after, crate::prog_opts::ORIENT_TAIL_AFTER_OFF);
        assert_eq!(crate::prog_opts::ORIENT_EVIDENCE_DEFAULT, None);
        assert!(
            Args::try_parse_from([
                "oarfish",
                "-a",
                "x.rad",
                "-o",
                "o",
                "--orient-evidence",
                "off"
            ])
            .is_err()
        );
    }

    #[test]
    fn end_fit_options_are_checked() {
        let parse = |extra: &[&str]| {
            let mut v = vec!["oarfish", "--mapper", "spearfish", "--reads", "r.fq"];
            v.extend_from_slice(&["--index", "t.sfi", "-o", "out", "--seq-tech", "ont-drna"]);
            v.extend_from_slice(extra);
            Args::try_parse_from(v)
        };
        assert!(parse(&["--end-fit", "-1"]).is_err());
        assert!(parse(&["--end-fit", "inf"]).is_err());
        assert!(parse(&["--end-fit-tol", "-1"]).is_err());
        assert!(parse(&["--end-fit-from", "both"]).is_err());
        // They need reads to map.
        assert!(
            Args::try_parse_from(["oarfish", "-a", "x.rad", "-o", "o", "--end-fit", "1"]).is_err()
        );
    }
}
