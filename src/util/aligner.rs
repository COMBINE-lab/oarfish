use crate::prog_opts::{Args, SequencingTech};
use crate::util::digest_utils;
use crate::util::file_utils::{get_ref_source, is_fasta};
use crate::util::oarfish_types::NamedDigestVec;

#[cfg(not(feature = "rammap"))]
use minimap2::Aligner;
#[cfg(not(feature = "rammap"))]
use minimap2_sys::MmIdx;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;

use num_format::{Locale, ToFormattedString};

#[cfg(not(feature = "rammap"))]
use std::ffi::CStr;
use std::fs::File;
use std::num::NonZeroUsize;
#[cfg(not(feature = "rammap"))]
use std::sync::Arc;

use tracing::{info, warn};

pub(crate) type HeaderReaderAlignerDigest = (
    noodles_sam::header::Header,
    Option<bam::io::Reader<bgzf::io::MultithreadedReader<File>>>,
    Option<crate::util::mapper::Mapper>,
    NamedDigestVec,
);

/// Obtain a `HeaderReaderAlignerDigest` object from the provided arguments. This
/// function properly dispatches on the appropriate argument type, so that it will just
/// *do the right thing* if the user provides an --annotated (and/or a --nove) set of sequences
/// or if they provide an existing index via --index.  It has all of the logic to verify
/// that the provided input configuration is correct, and to take the appropriate actions
/// (i.e. generate the correct digests) if it is.  It returns a `Result` so that an
/// `anyhow::Error` will be returned describing any error encountered during the process
/// of constructing the aligner, or an `Ok(Aligner)` if successful.
#[cfg(not(feature = "rammap"))]
pub(crate) fn get_aligner_from_args(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    info!("oarfish is operating in read-based mode");
    if args.index.is_some() {
        get_aligner_from_index(args)
    } else {
        get_aligner_from_fastas(args)
    }
}

/// The user provided FASTA files (possibly gzipped) via --annotated (and/or --novel), so
/// create the appropriate digests and build an aligner based on these sequences.
#[cfg(not(feature = "rammap"))]
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

    let ref_digest_handle = args
        .annotated
        .clone()
        .map(|refs| digest_utils::get_digest_from_fasta(&refs));

    let novel_digest_handle = args
        .novel
        .clone()
        .map(|refs| digest_utils::get_digest_from_fasta(&refs));

    // we are using either 1 or 2 background threads to compute the digests
    let thread_sub = if ref_digest_handle
        .as_ref()
        .xor(novel_digest_handle.as_ref())
        .is_some()
    {
        1
    } else {
        2
    };

    // set the number of indexing threads
    let idx_threads = &args.threads.saturating_sub(thread_sub).max(1);

    // if we need to combine the annotated and novel sequences into the index,
    // spawn off a thread to do that
    let input_source = get_ref_source(args.annotated.clone(), args.novel.clone())?;

    // if the user requested to write the output index to disk, prepare for that
    let idx_out_as_str = args.index_out.clone().map(|x| {
        x.to_str()
            .expect("could not convert PathBuf to &str")
            .to_owned()
    });
    let idx_output: Option<&str> = idx_out_as_str.as_deref();

    // create the aligner
    let mut aligner = make_aligner(
        &args.seq_tech,
        *idx_threads,
        input_source.file_path(),
        idx_output,
    )?;

    info!("created aligner index opts : {:?}", aligner.idxopt);
    // get up to the best_n hits for each read
    // default value is 100.
    aligner.mapopt.best_n = args.best_n as i32;
    // set the seed to be the same as what command-line
    // minimap2 uses.
    aligner.mapopt.seed = 11;

    let n_seq = aligner.n_seq();

    info!(
        "index contains {} sequences",
        n_seq.to_formatted_string(&Locale::en)
    );

    // if we concatenated reference and novel transcripts to build the index, remove the fifo here
    input_source.join_if_needed()?;

    let header = make_header(&mut aligner)?;

    let mut digests = NamedDigestVec::new();
    // we have a reference file
    if let Some(digest_handle_inner) = ref_digest_handle {
        let digest_res = digest_handle_inner.join().expect("valid digest");
        let digest = digest_res?;
        digests.push(("annotated_transcripts_digest".to_string(), digest));
    };
    // we have a novel transcripts file
    if let Some(digest_handle_inner) = novel_digest_handle {
        let digest_res = digest_handle_inner.join().expect("valid digest");
        let digest = digest_res?;
        digests.push(("novel_transcripts_digest".to_string(), digest));
    };

    // if we created an index, append the digest
    if let Some(idx_file) = idx_output {
        digest_utils::append_digest_footer(idx_file, &digests)?;
    }

    Ok((header, None, Some(aligner), digests))
}

/// The user provided an existing index, so build the proper aligner based on this index.
#[cfg(not(feature = "rammap"))]
fn get_aligner_from_index(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    let idx_file = args.index.clone().expect("index file should exist");

    warn!(
        "You are using an existing minimap2 index. This means that the parameters provided at index construction time will be applied."
    );
    warn!(
        "Thus, any *indexing-related* parameters implied by your `--seq-tech` option will be ignored. If you have not done this intentionally, please make sure the proper parameters were used when building the index."
    );

    // set the number of indexing threads
    let idx_threads = &args.threads;

    // if the user requested to write the output index to disk, prepare for that
    let idx_out_as_str = args.index_out.clone().map_or(String::new(), |x| {
        x.to_str()
            .expect("could not convert PathBuf to &str")
            .to_owned()
    });
    let idx_output = args.index_out.as_ref().map(|_| idx_out_as_str.as_str());

    // create the aligner
    let mut aligner = make_aligner(&args.seq_tech, *idx_threads, &idx_file, idx_output)?;

    info!("created aligner index opts : {:?}", aligner.idxopt);
    // get up to the best_n hits for each read
    // default value is 100.
    aligner.mapopt.best_n = args.best_n as i32;
    // set the seed to be the same as what command-line
    // minimap2 uses.
    aligner.mapopt.seed = 11;

    let n_seq = aligner.n_seq();

    info!(
        "index contains {} sequences",
        n_seq.to_formatted_string(&Locale::en)
    );

    let header = make_header(&mut aligner)?;

    let digests = match digest_utils::read_digest_footer(
        idx_file.to_str().expect("could not convert to string"),
    ) {
        // we read a pre-computed digest from an oarfish-constructed
        // minimap2 index
        Ok(d) => d,
        _ => {
            // We have been given a minimap2 index, but without the oarfish
            // footer. Now, we can build the digest we want from the index
            // itself.
            warn!(
                "computing sequence signatures from a minimap2 index that was not built with oarfish."
            );
            warn!(
                "if you are quantifying multiple samples, it will save time to let oarfish build a minimap2 index from the transcriptome reference, so that the reference signature can be reused."
            );
            let mmi: Arc<MmIdx> = Arc::clone(aligner.idx.as_ref().unwrap());
            let d = digest_utils::digest_from_index(&mmi)?;
            vec![("prexisting_mm2_index".to_string(), d)].into()
        }
    };

    Ok((header, None, Some(aligner), digests))
}

/// Build the actual aligner based on the sequencing technology, number
/// of threads, and the provided input (and optional output) file.
#[cfg(not(feature = "rammap"))]
fn make_aligner(
    seq_tech: &Option<SequencingTech>,
    idx_threads: usize,
    input: &std::path::Path,
    idx_output: Option<&str>,
) -> anyhow::Result<Aligner<minimap2::Built>> {
    // create the aligner
    match seq_tech {
        Some(SequencingTech::OntCDNA) | Some(SequencingTech::OntDRNA) => {
            minimap2::Aligner::builder()
                .map_ont()
                .with_index_threads(idx_threads)
                .with_cigar()
                .with_cigar_clipping()
                .with_index(input.to_path_buf().clone(), idx_output)
                .map_err(anyhow::Error::msg)
        }
        Some(SequencingTech::PacBio) => minimap2::Aligner::builder()
            .map_pb()
            .with_index_threads(idx_threads)
            .with_cigar()
            .with_cigar_clipping()
            .with_index(input.to_path_buf().clone(), idx_output)
            .map_err(anyhow::Error::msg),

        Some(SequencingTech::PacBioHifi) => minimap2::Aligner::builder()
            .map_hifi()
            .with_index_threads(idx_threads)
            .with_cigar()
            .with_cigar_clipping()
            .with_index(input.to_path_buf().clone(), idx_output)
            .map_err(anyhow::Error::msg),
        None => {
            anyhow::bail!("sequencing tech must be provided in read mode, but it was not!");
        }
    }
}

/// Make a proper header object from the index (this allows us to rely on the same
/// header structure in either raw sequence or alignment-based mode).
#[cfg(not(feature = "rammap"))]
fn make_header(
    aligner: &mut Aligner<minimap2::Built>,
) -> anyhow::Result<noodles_sam::header::Header> {
    let n_seq = aligner.n_seq();
    let mut header = noodles_sam::header::Header::builder();

    // TODO: better creation of the header
    {
        for i in 0..n_seq {
            let seq = aligner.get_seq(i as usize).unwrap_or_else(|| {
                panic!(
                    "{} was not a valid reference sequence index. (n_seq = {})",
                    i, n_seq
                )
            });
            let c_str = unsafe { CStr::from_ptr(seq.name) };
            let rust_str = c_str.to_str().unwrap().to_string();
            header = header.add_reference_sequence(
                rust_str,
                HeaderMap::<header_val::map::ReferenceSequence>::new(NonZeroUsize::try_from(
                    seq.len as usize,
                )?),
            );
        }
    }

    header = header.add_program(
        "minimap2-rs",
        HeaderMap::<header_val::map::Program>::default(),
    );

    Ok(header.build())
}

/// Build a *spliced* minimap2 aligner over a genome (FASTA or `.mmi` index) for
/// genome-read mode, and return the list of reference (chromosome) names in
/// target-id order. `refnames[i]` is the chromosome whose 0-based `target_id`
/// (and thus bramble `RefId`) is `i`.
///
/// The splice preset is chosen from the sequencing technology (`splice:hq` for
/// PacBio HiFi, `splice` otherwise). When `--junctions` is provided, the BED of
/// known junctions is loaded (equivalent to minimap2 `--junc-bed`).
#[cfg(not(feature = "rammap"))]
pub(crate) fn get_genome_aligner_from_args(
    args: &Args,
    junc_bed: Option<&std::path::Path>,
) -> anyhow::Result<(Aligner<minimap2::Built>, Vec<String>)> {
    let genome = args
        .genome
        .clone()
        .expect("--genome is required in genome read mode");
    info!(
        "oarfish is operating in genome read mode; building a spliced minimap2 index over {}",
        genome.display()
    );

    let idx_threads = args.threads.max(1);
    let builder = match args.seq_tech {
        Some(SequencingTech::PacBioHifi) => Aligner::builder().splice_hq(),
        Some(_) => Aligner::builder().splice(),
        None => anyhow::bail!("--seq-tech is required in genome read mode, but it was not provided!"),
    };

    let mut aligner = builder
        .with_index_threads(idx_threads)
        .with_cigar()
        .with_cigar_clipping()
        .with_index(genome.clone(), None)
        .map_err(anyhow::Error::msg)?;

    // up to best_n secondary mappings, and the same seed as command-line minimap2.
    aligner.mapopt.best_n = args.best_n as i32;
    aligner.mapopt.seed = 11;
    // Splice-junction loading (annotated and/or user-supplied) happens in the
    // caller, after the annotation transcripts have been parsed.

    // Load annotated splice junctions (minimap2 `--junc-bed`) if provided.
    if let Some(bed) = junc_bed {
        let bs = bed.to_str().expect("junction BED path must be valid UTF-8");
        aligner
            .read_junction_lr(bs)
            .map_err(|code| anyhow::anyhow!("failed to load splice junctions {} (code {})", bs, code))?;
        info!("loaded splice junctions into the genome index from {}", bs);
    }

    let refnames = genome_aligner_ref_names(&aligner);
    info!(
        "genome index contains {} reference sequences",
        refnames.len().to_formatted_string(&Locale::en)
    );
    Ok((aligner, refnames))
}

/// Extract the reference (chromosome) names from a built genome aligner, in
/// `target_id` order (0..n_seq).
#[cfg(not(feature = "rammap"))]
fn genome_aligner_ref_names(aligner: &Aligner<minimap2::Built>) -> Vec<String> {
    let n_seq = aligner.n_seq();
    let mut names = Vec::with_capacity(n_seq as usize);
    for i in 0..n_seq {
        let seq = aligner.get_seq(i as usize).unwrap_or_else(|| {
            panic!("{} was not a valid reference sequence index (n_seq = {})", i, n_seq)
        });
        let c_str = unsafe { CStr::from_ptr(seq.name) };
        names.push(c_str.to_str().unwrap().to_string());
    }
    names
}

// ===========================================================================
// rammap backend
// ===========================================================================
//
// Mirrors the minimap2 implementations above, but builds a pure-Rust `rammap`
// aligner. Only the FASTA (`--annotated`/`--novel`, `--genome`) construction
// path is supported; loading a pre-built minimap2 `.mmi` index (`--index`) is
// not, because rammap uses its own index format.

/// Map an oarfish read-mode sequencing tech onto a rammap transcriptome preset.
#[cfg(feature = "rammap")]
fn rammap_txome_preset(seq_tech: &Option<SequencingTech>) -> anyhow::Result<rammap::api::Preset> {
    match seq_tech {
        Some(SequencingTech::OntCDNA) | Some(SequencingTech::OntDRNA) => Ok(rammap::api::Preset::MapOnt),
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
#[cfg(feature = "rammap")]
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
#[cfg(feature = "rammap")]
fn rammap_make_header(aligner: &rammap::api::Aligner) -> anyhow::Result<noodles_sam::header::Header> {
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
/// both rammap's native RMMI format and minimap2 `.mmi` files (it detects the
/// magic and parses the minimap2 binary layout), so an index built by oarfish's
/// minimap2 backend (`--index-out`) can be reused here. As with the minimap2
/// backend, indexing-related `--seq-tech` parameters are ignored: the index's
/// own k/w are used. The oarfish digest footer (if present) is read back;
/// otherwise a name/length-based digest is recorded.
#[cfg(feature = "rammap")]
fn get_aligner_from_index(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    let idx_file = args.index.clone().expect("index file should exist");
    let idx_str = idx_file
        .to_str()
        .expect("could not convert index path to &str");

    warn!(
        "You are using an existing index with the rammap backend (rammap reads both its own RMMI \
         format and minimap2 .mmi files). Indexing-related parameters implied by --seq-tech are \
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
    // retains them), matching what the minimap2 backend does via
    // `digest_from_index`; only if the index is sequence-stripped do we fall
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

    Ok((header, None, Some(std::sync::Arc::new(aligner)), digests))
}

#[cfg(feature = "rammap")]
pub(crate) fn get_aligner_from_args(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
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

    // compute the reference digests on background threads (identical to the
    // minimap2 path, so the resulting signatures are comparable).
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
    // up to best_n secondary mappings; same seed as command-line minimap2.
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

    Ok((header, None, Some(std::sync::Arc::new(aligner)), digests))
}

/// Build a *spliced* rammap aligner over a genome (FASTA) for genome-read mode,
/// returning the reference (chromosome) names in `target_id` order. The splice
/// preset is chosen from the sequencing technology (`splice:hq` for PacBio HiFi,
/// `splice` otherwise).
///
/// Annotated splice junctions (when `junc_bed` is provided) are loaded via
/// rammap's `load_junctions_bed`, the equivalent of minimap2's `--junc-bed`.
#[cfg(feature = "rammap")]
pub(crate) fn get_genome_aligner_from_args(
    args: &Args,
    junc_bed: Option<&std::path::Path>,
) -> anyhow::Result<(crate::util::mapper::Mapper, Vec<String>)> {
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
        None => anyhow::bail!("--seq-tech is required in genome read mode, but it was not provided!"),
    };

    let genome_str = genome
        .to_str()
        .expect("could not convert genome path to &str");
    let mut aligner = rammap::api::Aligner::from_fasta(genome_str, preset)?;
    aligner.options_mut().filtering.best_n = args.best_n as i32;
    aligner.options_mut().filtering.seed = 11;

    // Load annotated splice junctions (rammap's `--junc-bed` equivalent) if provided.
    if let Some(bed) = junc_bed {
        let bs = bed.to_str().expect("junction BED path must be valid UTF-8");
        aligner
            .load_junctions_bed(bs)
            .map_err(|e| anyhow::anyhow!("failed to load splice junctions {}: {}", bs, e))?;
        info!("loaded splice junctions into the genome index from {}", bs);
    }

    let refnames: Vec<String> = aligner.index().seqs.iter().map(|t| t.name.clone()).collect();
    info!(
        "genome index contains {} reference sequences",
        refnames.len().to_formatted_string(&Locale::en)
    );
    Ok((std::sync::Arc::new(aligner), refnames))
}

/// Build a rammap index from the provided reference FASTA(s) and (when
/// `--index-out` is set) persist it to disk as a native RMMI file with the
/// oarfish reference-signature footer appended. This backs the `--only-index`
/// path and mirrors the minimap2 [`get_aligner_from_fastas`] above; the
/// resulting signature is computed from the FASTA via
/// [`digest_utils::get_digest_from_fasta`], so it is byte-identical to the
/// signature a minimap2 index built from the same FASTA would carry.
#[cfg(feature = "rammap")]
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

    // compute the reference digests on background threads (identical to the
    // minimap2 path, so the resulting signatures are comparable).
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

    Ok((header, None, Some(std::sync::Arc::new(aligner)), digests))
}

/// Persist a rammap aligner's index to `idx_out` as a native RMMI file and
/// append the oarfish reference-signature footer. The footer is tolerated by
/// rammap's loader (which stops at the end of the bincode structure), so the
/// index remains a single self-contained file that round-trips through
/// [`rammap::api::Aligner::from_index`] while carrying its signature.
#[cfg(feature = "rammap")]
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
#[cfg(feature = "rammap")]
pub(crate) fn rescue_db_from_genome_index(
    aligner: &rammap::api::Aligner,
) -> Option<bramble_rs::fasta::FastaDb> {
    use std::collections::HashMap;
    let idx = aligner.index();
    if !idx.has_sequences() {
        info!(
            "the genome index stores no reference sequences; soft-clip rescue is off \
             (pass --genome-fasta, or rebuild the index with sequences, to enable it)"
        );
        return None;
    }
    info!("sourcing soft-clip-rescue reference sequences from the loaded genome index");
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::with_capacity(idx.seqs.len());
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
    Some(bramble_rs::fasta::FastaDb::from_seqs(seqs))
}

/// minimap2 backend: `.mmi` sequence extraction is not implemented (the backend
/// is being retired). A FASTA `--genome`/`--genome-fasta` still provides the
/// rescue reference; a prebuilt `.mmi` simply leaves rescue off.
#[cfg(not(feature = "rammap"))]
pub(crate) fn rescue_db_from_genome_index(
    _aligner: &Aligner<minimap2::Built>,
) -> Option<bramble_rs::fasta::FastaDb> {
    None
}
