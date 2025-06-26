use crate::prog_opts::{Args, SequencingTech};
use crate::util::digest_utils;
use crate::util::file_utils::{get_ref_source, is_fasta};
use crate::util::oarfish_types::NamedDigestVec;

use minimap2::Aligner;
use minimap2_sys::MmIdx;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;

use num_format::{Locale, ToFormattedString};

use std::ffi;
use std::fs::File;
use std::num::NonZeroUsize;
use std::sync::Arc;

use tracing::{info, warn};

pub(crate) type HeaderReaderAlignerDigest = (
    noodles_sam::header::Header,
    Option<bam::io::Reader<bgzf::MultithreadedReader<File>>>,
    Option<minimap2::Aligner<minimap2::Built>>,
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
pub(crate) fn get_aligner_from_args(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    info!("oarfish is operating in read-based mode");
    if args.index.is_some() {
        get_aligner_from_index(args)
    } else {
        assert!(
            args.annotated
                .as_ref()
                .is_none_or(|f| is_fasta(f).expect("couldn't read input file."))
        );
        assert!(
            args.novel
                .as_ref()
                .is_none_or(|f| is_fasta(f).expect("couldn't read input file."))
        );
        get_aligner_from_fastas(args)
    }
}

/// The user provided FASTA files (possibly gzipped) via --annotated (and/or --novel), so
/// create the appropriate digests and build an aligner based on these sequences.
fn get_aligner_from_fastas(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
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
    let idx_out_as_str = args.index_out.clone().map_or(String::new(), |x| {
        x.to_str()
            .expect("could not convert PathBuf to &str")
            .to_owned()
    });
    let idx_output = args.index_out.as_ref().map(|_| idx_out_as_str.as_str());

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
        digest_utils::append_digest_to_mm2_index(idx_file, &digests)?;
    }

    Ok((header, None, Some(aligner), digests))
}

/// The user provided an existing index, so build the proper aligner based on this index.
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

    let digests = match digest_utils::read_digest_from_mm2_index(
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
            let c_str = unsafe { ffi::CStr::from_ptr(seq.name) };
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
