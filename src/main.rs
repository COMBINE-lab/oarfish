use anyhow::bail;
use clap::Parser;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use util::oarfish_types::NamedDigestVec;

use core::ffi;
use minimap2_sys::MmIdx;
// Or now
// use minimap2::ffi as mm_ffi;
//use minimap2_temp as minimap2;
use num_format::{Locale, ToFormattedString};
use std::io::Read;
use std::sync::Arc;
use std::{fs::File, io};
use tracing::{error, info, warn};
use tracing_subscriber::{EnvFilter, filter::LevelFilter, fmt, prelude::*};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::header::record::value as header_val;
use noodles_sam::header::record::value::Map as HeaderMap;

mod alignment_parser;
mod bootstrap;
mod bulk;
mod em;
mod prog_opts;
mod single_cell;
mod util;

use crate::prog_opts::{Args, FilterGroup, SequencingTech};
use crate::util::digest_utils;
use crate::util::file_utils::create_fifo_if_absent;
use crate::util::normalize_probability::normalize_read_probs;
use crate::util::oarfish_types::{AlignmentFilters, TranscriptInfo};
use crate::util::read_function::is_fasta;
use crate::util::{
    binomial_probability::binomial_continuous_prob, kde_utils, logistic_probability::logistic_prob,
};

type HeaderReaderAlignerDigest = (
    noodles_sam::header::Header,
    Option<bam::io::Reader<bgzf::MultithreadedReader<File>>>,
    Option<minimap2::Aligner<minimap2::Built>>,
    NamedDigestVec,
);

#[allow(dead_code)]
enum SourceType {
    Fastx(PathBuf),
    ExistingMM2Index(PathBuf),
    ExistingOarfishIndex(PathBuf),
}

impl SourceType {
    /// determine the type of the source and return the path
    /// wrapped in the correct variant.
    pub fn from_path<P: AsRef<Path>>(p: P) -> Self {
        if is_fasta(p.as_ref()).unwrap_or(false) {
            Self::Fastx(PathBuf::from(p.as_ref()))
        } else {
            match digest_utils::read_digest_from_mm2_index(
                p.as_ref().to_str().expect("can be represented as a str"),
            ) {
                // we read a pre-computed digest from an oarfish-constructed
                // minimap2 index
                Ok(_d) => Self::ExistingOarfishIndex(PathBuf::from(p.as_ref())),
                _ => Self::ExistingMM2Index(PathBuf::from(p.as_ref())),
            }
        }
    }

    #[allow(dead_code)]
    pub fn is_raw_mm2_index(&self) -> bool {
        matches!(&self, Self::ExistingMM2Index(_))
    }

    #[allow(dead_code)]
    pub fn is_oarfish_index(&self) -> bool {
        matches!(&self, Self::ExistingOarfishIndex(_))
    }

    #[allow(dead_code)]
    pub fn is_fasta(&self) -> bool {
        matches!(&self, Self::Fastx(_))
    }
}

fn get_digest_from_fasta(
    fpath: &std::path::PathBuf,
) -> std::thread::JoinHandle<anyhow::Result<seqcol_rs::DigestResult>> {
    let fpath_clone = fpath.clone();
    std::thread::spawn(|| {
        info!("generating reference digest for {}", fpath_clone.display());
        let mut seqcol_obj = seqcol_rs::SeqCol::try_from_fasta_file(fpath_clone).unwrap();
        let digest = seqcol_obj.digest(seqcol_rs::DigestConfig {
            level: seqcol_rs::DigestLevel::Level1,
            additional_attr: vec![seqcol_rs::KnownAttr::SortedNameLengthPairs],
        });
        info!("done");
        digest
    })
}

#[derive(Debug)]
struct RefSource {
    file_path: std::path::PathBuf,
    concat_handle: Option<std::thread::JoinHandle<anyhow::Result<()>>>,
}

fn get_ref_source(
    reference: Option<PathBuf>,
    novel_transcripts: Option<PathBuf>,
) -> anyhow::Result<RefSource> {
    let concat_handle: Option<std::thread::JoinHandle<_>>;

    if reference.as_ref().or(novel_transcripts.as_ref()).is_none() {
        bail!("at least one of --reference or --novel-transcripts but be provided");
    }

    // The `ref_file` input argument is either a FASTA file with reference
    // sequences, in which case we will compute the proper digest in a separate
    // thread, OR an existing minimap2 index, in which case we won't attempt
    // to treat it as a FASTA file and we will later get the digest from
    // the index.
    let input_path = if reference.is_some() && novel_transcripts.is_some() {
        let pid = std::process::id();
        let fifo_fname = format!("combine_transcripts_{}.fifo", pid);
        create_fifo_if_absent(&fifo_fname)?;

        let fifo_fname_clone = fifo_fname.to_string().clone();
        let ref_paths = reference.clone().expect("references should exist");
        let novel_paths = novel_transcripts.clone().expect("novel txps should exist");
        // a thread that will concatenate the reference transcripts and then the novel
        // trancsripts
        concat_handle = Some(std::thread::spawn(move || {
            let fifo_path = std::path::Path::new(fifo_fname_clone.as_str());
            let mut ff = std::fs::File::options().write(true).open(fifo_path)?;
            let ref_read = std::io::BufReader::new(std::fs::File::open(ref_paths)?);
            let novel_read = std::io::BufReader::new(std::fs::File::open(novel_paths)?);
            let mut reader = ref_read.chain(novel_read);
            info!("before copy");
            match std::io::copy(&mut reader, &mut ff) {
                Err(e) => {
                    error!("Error: {:#?}, in copying input to output", e);
                }
                Ok(nb) => {
                    info!("copied {} bytes from input to output across fifo", nb)
                }
            }
            drop(reader);
            drop(ff);
            Ok(())
        }));

        PathBuf::from(&fifo_fname)
    } else {
        concat_handle = None;
        reference
            .clone()
            .or(novel_transcripts.clone())
            .expect("either reference or novel transcripts must be provided")
    };

    Ok(RefSource {
        file_path: input_path.clone(),
        concat_handle,
    })
}

fn get_aligner_from_fastas(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    let ref_digest_handle = args
        .reference
        .clone()
        .map(|refs| get_digest_from_fasta(&refs));

    let novel_digest_handle = args
        .novel_transcripts
        .clone()
        .map(|refs| get_digest_from_fasta(&refs));

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

    // if we need to combine the reference and novel sequences into the index,
    // spawn off a thread to do that
    let input_source = get_ref_source(args.reference.clone(), args.novel_transcripts.clone())?;

    // if the user requested to write the output index to disk, prepare for that
    let idx_out_as_str = args.index_out.clone().map_or(String::new(), |x| {
        x.to_str()
            .expect("could not convert PathBuf to &str")
            .to_owned()
    });
    let idx_output = args.index_out.as_ref().map(|_| idx_out_as_str.as_str());

    // create the aligner
    let mut aligner = match args.seq_tech {
        Some(SequencingTech::OntCDNA) | Some(SequencingTech::OntDRNA) => {
            minimap2::Aligner::builder()
                .map_ont()
                .with_index_threads(*idx_threads)
                .with_cigar()
                .with_index(&input_source.file_path, idx_output)
                .expect("could not construct minimap2 index")
        }
        Some(SequencingTech::PacBio) => minimap2::Aligner::builder()
            .map_pb()
            .with_index_threads(*idx_threads)
            .with_cigar()
            .with_index(&input_source.file_path, idx_output)
            .expect("could not construct minimap2 index"),
        Some(SequencingTech::PacBioHifi) => minimap2::Aligner::builder()
            .map_hifi()
            .with_index_threads(*idx_threads)
            .with_cigar()
            .with_index(&input_source.file_path, idx_output)
            .expect("could not construct minimap2 index"),
        None => {
            anyhow::bail!("sequencing tech must be provided in read mode, but it was not!");
        }
    };

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
    if let Some(concat_handle) = input_source.concat_handle {
        if let Err(e) = concat_handle.join() {
            bail!(
                "Failed to concatenate reference and novel input transcript sequences: {:#?}",
                e
            );
        } else {
            info!("joined successfully!");
        }
        std::fs::remove_file(input_source.file_path)?;
    }

    let mut header = noodles_sam::header::Header::builder();

    #[derive(Debug, PartialEq, Eq)]
    pub struct SeqMetaData {
        pub name: String,
        pub length: u32,
        pub is_alt: bool,
    }

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

    let header = header.build();

    let mut digests = NamedDigestVec::new();
    // we have a reference file
    if let Some(digest_handle_inner) = ref_digest_handle {
        let digest_res = digest_handle_inner.join().expect("valid digest");
        let digest = digest_res?;
        digests.push(("reference_digest".to_string(), digest));
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

fn get_aligner_from_index(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    let idx_file = args.index.clone().expect("index file should exist");

    if let SourceType::ExistingMM2Index(_idx) = SourceType::from_path(&idx_file) {
        warn!(
            "You are using an existing minimap2 index (constructed outside of oarfish). This means that the parameters provided at index construction time will be applied."
        );
        warn!(
            "Thus, the parameters implied by your `--seq-tech` option will be ignored. If you have not done this intentionally, please make sure the proper parameters were used when building the index."
        );
    }
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
    let mut aligner = match args.seq_tech {
        Some(SequencingTech::OntCDNA) | Some(SequencingTech::OntDRNA) => {
            minimap2::Aligner::builder()
                .map_ont()
                .with_index_threads(*idx_threads)
                .with_cigar()
                .with_index(idx_file.clone(), idx_output)
                .expect("could not construct minimap2 index")
        }
        Some(SequencingTech::PacBio) => minimap2::Aligner::builder()
            .map_pb()
            .with_index_threads(*idx_threads)
            .with_cigar()
            .with_index(idx_file.clone(), idx_output)
            .expect("could not construct minimap2 index"),
        Some(SequencingTech::PacBioHifi) => minimap2::Aligner::builder()
            .map_hifi()
            .with_index_threads(*idx_threads)
            .with_cigar()
            .with_index(idx_file.clone(), idx_output)
            .expect("could not construct minimap2 index"),
        None => {
            anyhow::bail!("sequencing tech must be provided in read mode, but it was not!");
        }
    };

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

    let mut header = noodles_sam::header::Header::builder();

    #[derive(Debug, PartialEq, Eq)]
    pub struct SeqMetaData {
        pub name: String,
        pub length: u32,
        pub is_alt: bool,
    }

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

    let header = header.build();

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

fn get_aligner_from_args(args: &mut Args) -> anyhow::Result<HeaderReaderAlignerDigest> {
    info!("oarfish is operating in read-based mode");
    if args.index.is_some() {
        get_aligner_from_index(args)
    } else {
        assert!(
            args.reference
                .as_ref()
                .is_none_or(|f| is_fasta(f).expect("couldn't read input file."))
        );
        assert!(
            args.novel_transcripts
                .as_ref()
                .is_none_or(|f| is_fasta(f).expect("couldn't read input file."))
        );
        get_aligner_from_fastas(args)
    }
}

fn get_filter_opts(args: &Args) -> anyhow::Result<AlignmentFilters> {
    // set all of the filter options that the user
    // wants to apply.
    match args.filter_group {
        Some(FilterGroup::NoFilters) => {
            info!("disabling alignment filters.");
            // override individual parameters if the user passed them in explicitly
            let fpc = args
                .five_prime_clip
                .provided_or_u32("overriding 5' clip with user-provided value", u32::MAX);
            let tpc = args
                .three_prime_clip
                .provided_or_i64("overriding 3' clip with user-provided value", i64::MAX);
            let st = args
                .score_threshold
                .provided_or_f32("overriding score threshold with user-provided value", 0_f32);
            let maf = args.min_aligned_fraction.provided_or_f32(
                "overriding min aligned fraction with user-provided value",
                0_f32,
            );
            let mal = args.min_aligned_len.provided_or_u32(
                "overriding min aligned length with user-provided value",
                1_u32,
            );

            Ok(AlignmentFilters::builder()
                .five_prime_clip(fpc)
                .three_prime_clip(tpc)
                .score_threshold(st)
                .min_aligned_fraction(maf)
                .min_aligned_len(mal)
                .which_strand(args.strand_filter)
                .model_coverage(args.model_coverage)
                .logistic_growth_rate(args.growth_rate)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .build())
        }
        Some(FilterGroup::NanocountFilters) => {
            info!("setting filters to nanocount defaults.");
            // override individual parameters if the user passed them in explicitly
            let fpc = args
                .five_prime_clip
                .provided_or_u32("overriding 5' clip with user-provided value", u32::MAX);
            let tpc = args
                .three_prime_clip
                .provided_or_i64("overriding 3' clip with user-provided value", 50_i64);
            let st = args.score_threshold.provided_or_f32(
                "overriding score threshold with user-provided value",
                0.95_f32,
            );
            let maf = args.min_aligned_fraction.provided_or_f32(
                "overriding min aligned fraction with user-provided value",
                0.5_f32,
            );
            let mal = args.min_aligned_len.provided_or_u32(
                "overriding min aligned length with user-provided value",
                50_u32,
            );

            Ok(AlignmentFilters::builder()
                .five_prime_clip(fpc)
                .three_prime_clip(tpc)
                .score_threshold(st)
                .min_aligned_fraction(maf)
                .min_aligned_len(mal)
                .which_strand(bio_types::strand::Strand::Forward)
                .model_coverage(args.model_coverage)
                .logistic_growth_rate(args.growth_rate)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .build())
        }
        None => {
            info!("setting user-provided filter parameters.");
            Ok(AlignmentFilters::builder()
                .five_prime_clip(args.five_prime_clip.try_as_u32()?)
                .three_prime_clip(args.three_prime_clip.try_as_i64()?)
                .score_threshold(args.score_threshold.try_as_f32()?)
                .min_aligned_fraction(args.min_aligned_fraction.try_as_f32()?)
                .min_aligned_len(args.min_aligned_len.try_as_u32()?)
                .which_strand(args.strand_filter)
                .model_coverage(args.model_coverage)
                .logistic_growth_rate(args.growth_rate)
                .write_assignment_probs(args.write_assignment_probs.is_some())
                .write_assignment_probs_type(args.write_assignment_probs.clone())
                .build())
        }
    }
}

fn main() -> anyhow::Result<()> {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();
    let (filtered_layer, reload_handle) = tracing_subscriber::reload::Layer::new(env_filter);

    // set up the logging.  Here we will take the
    // logging level from the environment variable if
    // it is set.  Otherwise, we'll set the default
    tracing_subscriber::registry()
        // log level to INFO.
        .with(fmt::layer().with_writer(io::stderr))
        .with(filtered_layer)
        .init();

    let mut args = Args::parse();

    // change the logging filter if the user specified quiet or
    // verbose.
    if args.quiet {
        reload_handle.modify(|filter| *filter = EnvFilter::new("WARN"))?;
    }
    if args.verbose {
        reload_handle.modify(|filter| *filter = EnvFilter::new("TRACE"))?;
    }

    let filter_opts = get_filter_opts(&args)?;

    let (header, reader, aligner, digest) = if args.alignments.is_none() {
        get_aligner_from_args(&mut args)?
    } else {
        let alignments = args.alignments.clone().unwrap();
        let afile = File::open(&alignments)?;

        let decomp_threads = if args.single_cell {
            // we will overlap quantification with parsing, so don't try to use too many
            // parser threads, and adjust the worker threads accordingly.

            // is there a better heuristic than this?
            // <= 6 threads, use only 1 for decompression
            // 6-8 threads, use 2 for decompression
            // > 8 threads, use 3 for decompression
            match args.threads {
                1..=6 => 1,
                7 | 8 => 2,
                _ => 3,
            }
        } else {
            // try to use all but 1 thread, and assume we have at least 2.
            1.max(args.threads.saturating_sub(1))
        };

        let worker_count = NonZeroUsize::new(decomp_threads).expect("decompression threads >= 1");
        if args.single_cell {
            args.threads = 1.max(args.threads.saturating_sub(decomp_threads));
        }

        let decoder = bgzf::MultithreadedReader::with_worker_count(worker_count, afile);
        let mut reader = bam::io::Reader::from(decoder);
        // parse the header, and ensure that the reads were mapped with minimap2 (as far as we
        // can tell).
        let header = alignment_parser::read_and_verify_header(&mut reader, &alignments)?;
        let seqcol_digest = digest_utils::digest_from_header(&header)?;
        (
            header,
            Some(reader),
            None,
            vec![("bam_digest".to_string(), seqcol_digest)].into(),
        )
    };

    let num_ref_seqs = header.reference_sequences().len();

    // where we'll write down the per-transcript information we need
    // to track.
    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(num_ref_seqs);
    let mut txps_name: Vec<String> = Vec::with_capacity(num_ref_seqs);

    // loop over the transcripts in the header and fill in the relevant
    // information here.
    if args.model_coverage {
        for (rseq, rmap) in header.reference_sequences().iter() {
            txps.push(TranscriptInfo::with_len_and_bin_width(
                rmap.length(),
                args.bin_width,
            ));
            txps_name.push(rseq.to_string());
        }
    } else {
        for (rseq, rmap) in header.reference_sequences().iter() {
            txps.push(TranscriptInfo::with_len(rmap.length()));
            txps_name.push(rseq.to_string());
        }
    }
    info!(
        "parsed reference information for {} transcripts.",
        txps.len().to_formatted_string(&Locale::en)
    );

    if args.single_cell {
        // TODO: do this better (quiet the EM during single-cell quant)
        reload_handle.modify(|filter| {
            *filter = if args.quiet {
                EnvFilter::new("WARN")
                    .add_directive("oarfish=warn".parse().unwrap())
                    .add_directive("oarfish::single_cell=warn".parse().unwrap())
            } else if args.verbose {
                EnvFilter::new("TRACE")
                    .add_directive("oarfish=info".parse().unwrap())
                    .add_directive("oarfish::single_cell=trace".parse().unwrap())
            } else {
                // be quiet about normal things in single-cell mode
                // e.g. EM iterations, and only print out info for
                // oarfish::single_cell events.
                EnvFilter::new("INFO")
                    .add_directive("oarfish=warn".parse().unwrap())
                    .add_directive("oarfish::single_cell=info".parse().unwrap())
            }
        })?;

        single_cell::quantify_single_cell_from_collated_bam(
            &header,
            &filter_opts,
            &mut reader.unwrap(),
            &mut txps,
            &args,
            digest,
        )?;
    } else if args.alignments.is_some() {
        bulk::quantify_bulk_alignments_from_bam(
            &header,
            filter_opts,
            &mut reader.unwrap(),
            &mut txps,
            &txps_name,
            &args,
            digest,
        )?;
    } else {
        bulk::quantify_bulk_alignments_raw_reads(
            &header,
            aligner.expect("need valid alinger to align reads"),
            filter_opts,
            &args.reads.clone().expect("expected read file(s)"),
            &mut txps,
            &txps_name,
            &args,
            digest,
        )?;
    }

    info!("oarfish completed successfully.");
    Ok(())
}
