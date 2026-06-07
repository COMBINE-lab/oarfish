use clap::{Parser, builder::ArgPredicate};
use serde::Serialize;
use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;
use tracing::info;

/// These represent different "meta-options", specific settings
/// for all of the different filters that should be applied in
/// different cases.
#[derive(Clone, Debug, clap::ValueEnum, Serialize)]
pub enum FilterGroup {
    NoFilters,
    NanocountFilters,
}

fn parse_strand(arg: &str) -> anyhow::Result<bio_types::strand::Strand> {
    match arg {
        "+" | "fw" | "FW" | "f" | "F" => Ok(bio_types::strand::Strand::Forward),
        "-" | "rc" | "RC" | "r" | "R" => Ok(bio_types::strand::Strand::Reverse),
        "." | "both" | "either" => Ok(bio_types::strand::Strand::Unknown),
        _ => anyhow::bail!("Cannot parse {} as a valid strand type", arg),
    }
}

/// Parse a strictly-positive `f32` (used for parameters that divide, so 0 and
/// negatives are rejected to avoid degenerate probabilities).
fn parse_pos_f32(arg: &str) -> anyhow::Result<f32> {
    let v: f32 = arg
        .parse()
        .map_err(|_| anyhow::anyhow!("`{}` is not a valid number", arg))?;
    if v > 0.0 {
        Ok(v)
    } else {
        anyhow::bail!("value must be > 0, but got {}", v)
    }
}

#[derive(Debug, Clone, clap::ValueEnum, Serialize)]
pub enum ReadAssignmentProbOut {
    Uncompressed,
    Compressed,
}

/// Which signal to weight a read's projected (genome-mode) alignments by when
/// forming conditional probabilities for the EM.
#[derive(Debug, Clone, Copy, PartialEq, clap::ValueEnum, Serialize)]
pub enum ProjProbSource {
    /// bramble exonic-coverage similarity only: `exp((sim - best_sim) * beta)`.
    Similarity,
    /// aligner genomic-alignment score only: `exp((score - best_score) / 5)`.
    /// Discriminates paralogous loci better than similarity.
    Score,
    /// combine both: `exp((score - best_score)/5 + beta*(sim - best_sim))`.
    /// Score separates loci; similarity separates isoforms within a locus.
    Combined,
}

fn parse_assign_prob_out_value(s: &str) -> anyhow::Result<ReadAssignmentProbOut> {
    match s.to_lowercase().as_str() {
        "raw" => Ok(ReadAssignmentProbOut::Uncompressed),
        "uncompressed" => Ok(ReadAssignmentProbOut::Uncompressed),
        "compressed" => Ok(ReadAssignmentProbOut::Compressed),
        "lz4" => Ok(ReadAssignmentProbOut::Compressed),
        x => anyhow::bail!(
            "Cannot parse {} as a valid option for read assignment probability output",
            x
        ),
    }
}

fn parse_display_thresh(s: &str) -> anyhow::Result<f64> {
    match s.to_lowercase().as_str() {
        "none" => Ok(f64::MIN_POSITIVE),
        _ => {
            let val = s.parse::<f64>()?;
            if (0.0..=1.0).contains(&val) {
                Ok(val)
            } else {
                anyhow::bail!("display-thresh must be between 0.0 and 1.0, got {}", val)
            }
        }
    }
}

#[derive(Debug, Clone, clap::ValueEnum, Serialize)]
pub enum SequencingTech {
    OntCDNA,
    OntDRNA,
    PacBio,
    PacBioHifi,
}

impl FromStr for SequencingTech {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "ont" => Ok(SequencingTech::OntCDNA),
            "ont-cdna" => Ok(SequencingTech::OntCDNA),
            "ont-drna" => Ok(SequencingTech::OntDRNA),
            "pb" => Ok(SequencingTech::PacBio),
            "pacbio" => Ok(SequencingTech::PacBio),
            "pb-hifi" => Ok(SequencingTech::PacBioHifi),
            "pacbio-hifi" => Ok(SequencingTech::PacBioHifi),
            x => Err(format!("Unknown protocol type {:}", x)),
        }
    }
}

/// This tells us the value of the filter argument and
/// the type remembers if it was the default or if the
/// user provided it explicltiy.
/// TODO: see if there is some built-in clap functionality
/// to avoid this song and dance.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum FilterArg {
    DefaultI64(i64),
    ProvidedI64(i64),
    DefaultU32(u32),
    ProvidedU32(u32),
    DefaultF32(f32),
    ProvidedF32(f32),
}

/// because we have to (in the derive approach at least) round trip
/// the default values through strings, we need some way to designate
/// a default value from a provided one. Default values will all start
/// with '*'
const DEFAULT_FILTER_PREFIX: &str = "*";

/// How to convert a FilterArg to a string
impl fmt::Display for FilterArg {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FilterArg::DefaultI64(x) => write!(f, "{}{}", DEFAULT_FILTER_PREFIX, x),
            FilterArg::DefaultU32(x) => write!(f, "{}{}", DEFAULT_FILTER_PREFIX, x),
            FilterArg::DefaultF32(x) => write!(f, "{}{}", DEFAULT_FILTER_PREFIX, x),
            FilterArg::ProvidedI64(x) => write!(f, "{}", x),
            FilterArg::ProvidedU32(x) => write!(f, "{}", x),
            FilterArg::ProvidedF32(x) => write!(f, "{}", x),
        }
    }
}

/// Necessary functions on [FilterArg]
impl FilterArg {
    /// If it is an i64 type, get the i64
    pub fn try_as_i64(&self) -> anyhow::Result<i64> {
        match self {
            FilterArg::DefaultI64(x) => Ok(*x),
            FilterArg::ProvidedI64(x) => Ok(*x),
            _ => anyhow::bail!("Could not provide FilterArg variant as an i64"),
        }
    }

    /// If it is an u32 type, get the u32
    pub fn try_as_u32(&self) -> anyhow::Result<u32> {
        match self {
            FilterArg::DefaultU32(x) => Ok(*x),
            FilterArg::ProvidedU32(x) => Ok(*x),
            _ => anyhow::bail!("Could not provide FilterArg variant as a u32"),
        }
    }

    /// If it is an f32 type, get the f32
    pub fn try_as_f32(&self) -> anyhow::Result<f32> {
        match self {
            FilterArg::DefaultF32(x) => Ok(*x),
            FilterArg::ProvidedF32(x) => Ok(*x),
            _ => anyhow::bail!("Could not provide FilterArg variant as an f32"),
        }
    }

    /// If the value is user provided, return the value, otherwise
    /// return `other` and print the message provided in `msg` to the
    /// logger.
    pub fn provided_or_u32(&self, msg: &str, other: u32) -> u32 {
        match self {
            FilterArg::ProvidedU32(x) => {
                info!("{} {}", msg, x);
                *x
            }
            _ => other,
        }
    }

    /// If the value is user provided, return the value, otherwise
    /// return `other` and print the message provided in `msg` to the
    /// logger.
    pub fn provided_or_i64(&self, msg: &str, other: i64) -> i64 {
        match self {
            FilterArg::ProvidedI64(x) => {
                info!("{} {}", msg, x);
                *x
            }
            _ => other,
        }
    }

    /// If the value is user provided, return the value, otherwise
    /// return `other` and print the message provided in `msg` to the
    /// logger.
    pub fn provided_or_f32(&self, msg: &str, other: f32) -> f32 {
        match self {
            FilterArg::ProvidedF32(x) => {
                info!("{} {}", msg, x);
                *x
            }
            _ => other,
        }
    }
}

/// Parse a string as a [FilterArg] with an i64 inner type
fn parse_filter_i64(arg: &str) -> anyhow::Result<FilterArg> {
    if let Some(val) = arg.strip_prefix(DEFAULT_FILTER_PREFIX) {
        let v = val.parse::<i64>()?;
        Ok(FilterArg::DefaultI64(v))
    } else {
        let v = arg.parse::<i64>()?;
        Ok(FilterArg::ProvidedI64(v))
    }
}

/// Parse a string as a [FilterArg] with a u32 inner type
fn parse_filter_u32(arg: &str) -> anyhow::Result<FilterArg> {
    if let Some(val) = arg.strip_prefix(DEFAULT_FILTER_PREFIX) {
        let v = val.parse::<u32>()?;
        Ok(FilterArg::DefaultU32(v))
    } else {
        let v = arg.parse::<u32>()?;
        Ok(FilterArg::ProvidedU32(v))
    }
}

/// Parse a string as a [FilterArg] with an f32 inner type
fn parse_filter_f32(arg: &str) -> anyhow::Result<FilterArg> {
    if let Some(val) = arg.strip_prefix(DEFAULT_FILTER_PREFIX) {
        let v = val.parse::<f32>()?;
        Ok(FilterArg::DefaultF32(v))
    } else {
        let v = arg.parse::<f32>()?;
        Ok(FilterArg::ProvidedF32(v))
    }
}

/// accurate transcript quantification from long-read RNA-seq data
#[derive(Parser, Debug, Serialize)]
#[clap(author, version, about, long_about = None)]
#[command(group(
    clap::ArgGroup::new("input")
    .required(true)
    .args(["alignments", "reads", "only_index", "genome_alignments"])
))]
#[command(group(
    clap::ArgGroup::new("raw_input_type")
    .multiple(true)
    .args(["annotated", "novel"])
))]
#[command(group(
    clap::ArgGroup::new("raw_ref_type")
    .multiple(true)
    .args(["annotated", "novel", "index"])
))]
// Any reference that `--reads` can be mapped against. In transcriptome mode
// this is one of annotated/novel/index; in genome (spliced) mode it is
// `--genome`. `--reads` requires at least one member of this group.
#[command(group(
    clap::ArgGroup::new("read_ref")
    .multiple(true)
    .args(["annotated", "novel", "index", "genome"])
))]
pub struct Args {
    /// be quiet (i.e. don't output log messages that aren't at least warnings)
    #[arg(long, conflicts_with = "verbose")]
    pub quiet: bool,

    /// be verbose (i.e. output all non-developer logging messages)
    #[arg(long)]
    pub verbose: bool,

    /// path to the file containing the input alignments
    #[arg(short, long, help_heading = "alignment mode")]
    pub alignments: Option<PathBuf>,

    /// path to the file containing the input reads; these can be
    /// in FASTA/Q format (possibly gzipped), or provided in
    /// uBAM (unaligned BAM) format. The format will be inferred from
    /// the file suffixes, and if a format cannot be inferred, it will
    /// be assumed to be (possibly gzipped) FASTA/Q
    #[arg(
        long,
        help_heading = "raw read mode",
        value_delimiter = ',',
        requires_ifs([
            (ArgPredicate::IsPresent, "read_ref"),
            (ArgPredicate::IsPresent, "seq_tech"),
        ])
    )]
    pub reads: Option<Vec<PathBuf>>,

    /// path to the file containing the annotated transcriptome (e.g. GENCODE) against which
    /// to map.
    #[arg(long, conflicts_with = "alignments", help_heading = "raw read mode")]
    pub annotated: Option<PathBuf>,

    /// path to the file containing novel (de novo, or reference-guided assembled) transcripts against which
    /// to map. These are ultimately indexed together with reference transcripts, but passed in
    /// separately for the purposes of provenance tracking.
    #[arg(long, conflicts_with = "alignments", help_heading = "raw read mode")]
    pub novel: Option<PathBuf>,

    /// path to an existing index (an oarfish-created index, which is preferred, or a
    /// compatible prebuilt index)
    #[arg(
        long,
        conflicts_with_all = ["alignments", "raw_input_type"],
        help_heading = "raw read mode"
    )]
    pub index: Option<PathBuf>,

    /// path to a spliced, genome-aligned BAM file. The alignments are projected
    /// onto the transcripts in `--annotation` (using bramble) and then
    /// quantified. The BAM must be collated by read name (e.g. `samtools collate`).
    #[arg(
        long,
        help_heading = "genome mode",
        requires = "annotation",
        conflicts_with_all = ["alignments", "reads", "only_index", "annotated", "novel", "index", "genome"]
    )]
    pub genome_alignments: Option<PathBuf>,

    /// path to a genome FASTA or prebuilt genome index. Used with `--reads` to perform
    /// spliced alignment of the reads to the genome, followed by projection onto
    /// the transcripts in `--annotation`.
    #[arg(
        long,
        help_heading = "genome mode",
        requires = "annotation",
        conflicts_with_all = ["alignments", "only_index", "raw_ref_type"]
    )]
    pub genome: Option<PathBuf>,

    /// path to a transcript annotation (GTF/GFF). Required in genome mode; the
    /// genomic alignments are projected onto the transcripts it defines.
    #[arg(long, help_heading = "genome mode")]
    pub annotation: Option<PathBuf>,

    /// genome FASTA providing the reference sequence for bramble's soft-clip
    /// rescue during projection. Soft-clip rescue is ON by default in genome
    /// mode; in read mode the sequence is taken from `--genome` automatically
    /// when it is a FASTA, so this is only needed to override that or to supply
    /// the reference in genome-alignments (BAM) mode. See `--no-rescue` to disable.
    #[arg(long, help_heading = "genome mode")]
    pub genome_fasta: Option<PathBuf>,

    /// disable bramble's soft-clip rescue during projection (genome mode). Rescue
    /// re-aligns soft-clipped read ends against the transcript's neighboring exons
    /// to recover discriminating sequence; it is on by default (improves isoform
    /// accuracy at negligible cost). Use this to turn it off.
    #[arg(long, help_heading = "genome mode")]
    pub no_rescue: bool,

    /// optional BED12 file of splice junctions / transcript models used to hint
    /// the spliced genome alignment (genome read mode). If omitted, junctions are
    /// derived automatically from `--annotation`; if provided, this file is used
    /// instead.
    #[arg(long, help_heading = "genome mode", requires = "genome")]
    pub junctions: Option<PathBuf>,

    /// do NOT derive splice junctions from `--annotation` for spliced genome
    /// alignment (genome read mode). By default annotated junctions are used, as
    /// they improve alignment; this flag disables that (e.g. for ablation).
    #[arg(long, help_heading = "genome mode", hide = true)]
    pub ignore_annotation_junctions: bool,

    /// spread (`beta`) used to convert bramble similarity scores into alignment
    /// probabilities in genome mode: `prob = exp((sim - best_sim) * beta)`.
    #[arg(long, hide = true, default_value_t = 10.0)]
    pub projected_prob_beta: f32,

    /// which signal weights a read's projected alignments in genome mode:
    /// `similarity` (bramble exonic coverage), `score` (aligner alignment
    /// score), or `combined`.
    #[arg(long, hide = true, value_enum, default_value_t = ProjProbSource::Similarity)]
    pub projected_prob_source: ProjProbSource,

    /// denominator `D` in the score→probability conversion `exp((score - best)/D)`
    /// used to weight a read's alignments in the EM (default 5). Larger `D`
    /// flattens the weighting across alignments of differing score; smaller `D`
    /// sharpens it toward the best-scoring alignment. Transcriptome mode only —
    /// in genome mode projected alignments are weighted by bramble similarity, so
    /// passing this with `--genome`/`--genome-alignments` is an error.
    #[arg(long, value_parser = parse_pos_f32, help_heading = "filters")]
    pub score_prob_denom: Option<f32>,

    /// genome mode: per-internal-junction-mismatch discount in (0,1] applied to a
    /// transcript's projection similarity (sharpens isoform discrimination).
    /// 1.0 = off (default).
    #[arg(long, hide = true, default_value_t = 1.0)]
    pub junc_miss_discount: f64,

    /// If this flag is passed, oarfish only performs indexing and not quantification.
    /// Designed primarily for workflow management systems.
    /// Note: A prebuilt index is not needed to quantify with oarfish; an index can be
    /// written concurrently with quantification using the `--index-out` parameter.
    #[arg(long, help_heading = "indexing")]
    pub only_index: bool,

    /// path where the index will be written (if provided)
    #[arg(long, conflicts_with_all = ["alignments", "index"], requires_ifs([(ArgPredicate::IsPresent, "only_index")]), help_heading = "indexing")]
    pub index_out: Option<PathBuf>,

    /// sequencing technology in which to expect reads if using mapping based mode
    #[arg(
        long,
        help_heading = "raw read mode",
        required_unless_present_any = ["alignments", "genome_alignments"],
        value_parser = clap::value_parser!(SequencingTech)
    )]
    pub seq_tech: Option<SequencingTech>,

    /// maximum number of secondary mappings to consider when mapping reads to the transcriptome
    #[arg(
        long,
        default_value_t = 100,
        requires = "reads",
        help_heading = "raw read mode"
    )]
    pub best_n: usize,

    /// cap (in MB) on the per-thread DP alignment-scratch buffer the rammap mapper
    /// retains between reads; bounds peak RSS at high thread counts. Unset = use the
    /// default (128 MB) or the `RAMMAP_DP_CACHE_CAP_MB` env var; `0` = unbounded
    /// (largest buffer pinned per thread, original behavior).
    #[arg(long, conflicts_with = "alignments", help_heading = "raw read mode")]
    pub dp_cache_cap_mb: Option<usize>,

    /// location where output quantification file should be written
    #[arg(short, long, required_unless_present = "only_index")]
    pub output: Option<PathBuf>,

    #[arg(long, help_heading = "filters", value_enum)]
    pub filter_group: Option<FilterGroup>,

    /// maximum allowable distance of the right-most end of an alignment from the 3' transcript end
    #[arg(short, long, help_heading="filters", default_value_t = FilterArg::DefaultI64(u32::MAX as i64), value_parser = parse_filter_i64)]
    pub three_prime_clip: FilterArg,

    /// maximum allowable distance of the left-most end of an alignment from the 5' transcript end
    #[arg(short, long, help_heading="filters", default_value_t = FilterArg::DefaultU32(u32::MAX), value_parser = parse_filter_u32)]
    pub five_prime_clip: FilterArg,

    /// fraction of the best possible alignment score that a secondary alignment must have for
    /// consideration
    #[arg(short, long, help_heading = "filters", default_value_t = FilterArg::DefaultF32(0.95), value_parser = parse_filter_f32)]
    pub score_threshold: FilterArg,

    /// fraction of a query that must be mapped within an alignemnt to consider the alignemnt
    /// valid
    #[arg(short, long, help_heading = "filters", default_value_t = FilterArg::DefaultF32(0.5), value_parser = parse_filter_f32)]
    pub min_aligned_fraction: FilterArg,

    /// minimum number of nucleotides in the aligned portion of a read
    #[arg(short = 'l', long, help_heading = "filters", default_value_t = FilterArg::DefaultU32(50), value_parser = parse_filter_u32)]
    pub min_aligned_len: FilterArg,

    /// only alignments to this strand will be allowed; options are (fw /+, rc/-, or both/.)
    #[arg(
        short = 'd',
        long,
        help_heading = "filters",
        default_value_t = bio_types::strand::Strand::Unknown,
        value_parser = parse_strand
    )]
    pub strand_filter: bio_types::strand::Strand,

    /// input is assumed to be a single-cell, transcriptome-aligned BAM (passed
    /// with `--alignments`) and to have the `CB:z` tag for all read records.
    /// Single-cell mode is not supported for raw reads or the genome (projection)
    /// modes, so it requires `--alignments` and conflicts with `--reads`,
    /// `--genome`, and `--genome-alignments`.
    #[arg(
        long,
        requires = "alignments",
        conflicts_with_all = ["reads", "genome", "genome_alignments"]
    )]
    pub single_cell: bool,

    /// apply the coverage model
    #[arg(long, help_heading = "coverage model", value_parser)]
    pub model_coverage: bool,

    /// if using the coverage model, use this as the value of `k` in the logistic equation
    #[arg(
        short = 'k',
        long,
        help_heading = "coverage model",
        value_parser,
        default_value_t = 2.0
    )]
    pub growth_rate: f64,

    /// write output alignment probabilites (optionally compressed) for each mapped read.
    /// If <WRITE_ASSIGNMENT_PROBS> is present, it must be one of `uncompressed` (default) or
    /// `compressed`, which will cause the output file to be lz4 compressed.
    #[arg(
        long,
        help_heading = "output read-txps probabilities",
        conflicts_with = "single_cell",
        default_missing_value = "uncompressed",
        num_args = 0..=1,
        require_equals = true,
        value_parser = parse_assign_prob_out_value
    )]
    pub write_assignment_probs: Option<ReadAssignmentProbOut>,

    /// minimum posterior probability threshold for a read-transcript assignment to be written
    /// to the .prob file. Accepts a float value between 0.0 and 1.0, or the sentinel value
    /// 'none' to use the minimum machine precision (f64::MIN_POSITIVE).
    #[arg(
        long,
        help_heading = "output read-txps probabilities",
        default_value_t = 1e-6,
        value_parser = parse_display_thresh
    )]
    pub display_thresh: f64,

    /// maximum number of iterations for which to run the EM algorithm
    #[arg(long, help_heading = "EM", default_value_t = 1000)]
    pub max_em_iter: u32,

    /// maximum number of iterations for which to run the EM algorithm
    #[arg(long, help_heading = "EM", default_value_t = 1e-3)]
    pub convergence_thresh: f64,

    /// number of cores that oarfish will use during different phases
    /// of quantification. Note: This value will be at least 2 for bulk
    /// quantification and at least 3 for single-cell quantification due to
    /// the use of dedicated parsing threads.
    #[arg(short = 'j', long, default_value_t = 3)]
    pub threads: usize,

    /// location of short read quantification (if provided)
    #[arg(short = 'q', long, help_heading = "EM")]
    pub short_quant: Option<String>,

    /// number of bootstrap replicates to produce to assess quantification uncertainty
    #[arg(long, default_value_t = 0)]
    pub num_bootstraps: u32,

    /// width of the bins used in the coverage model
    #[arg(short, long, help_heading = "coverage model", default_value_t = 100)]
    pub bin_width: u32,

    /// Number of alignment records to check for name collation when attempting
    /// to validate that the input BAM is name collated.
    #[arg(long, hide = true, default_value_t = 100_000)]
    pub sort_check_num: usize,

    /// use a KDE model of the observed fragment length distribution
    #[arg(short, long, hide = true)]
    pub use_kde: bool,
}

#[cfg(test)]
mod tests {
    use super::Args;
    use clap::Parser;

    #[test]
    fn allows_annotated_and_novel_together() {
        let parsed = Args::try_parse_from([
            "oarfish",
            "--reads",
            "reads.fq.gz",
            "--annotated",
            "annotated.fa",
            "--novel",
            "novel.fa",
            "--seq-tech",
            "ont-cdna",
            "-o",
            "out",
        ]);

        assert!(parsed.is_ok(), "expected parse success, got {parsed:?}");
        let parsed = parsed.expect("valid arguments");
        assert_eq!(
            parsed.annotated.as_deref(),
            Some(std::path::Path::new("annotated.fa"))
        );
        assert_eq!(
            parsed.novel.as_deref(),
            Some(std::path::Path::new("novel.fa"))
        );
        assert!(parsed.index.is_none());
    }

    #[test]
    fn rejects_index_with_raw_reference_fastas() {
        let parsed = Args::try_parse_from([
            "oarfish",
            "--reads",
            "reads.fq.gz",
            "--annotated",
            "annotated.fa",
            "--index",
            "transcripts.mmi",
            "--seq-tech",
            "ont-cdna",
            "-o",
            "out",
        ]);

        assert!(parsed.is_err(), "expected parse failure");
    }

    #[test]
    fn genome_alignments_requires_annotation() {
        // genome-alignments mode without --annotation must fail
        let missing = Args::try_parse_from([
            "oarfish",
            "--genome-alignments",
            "aln.genome.bam",
            "-o",
            "out",
        ]);
        assert!(missing.is_err(), "expected failure without --annotation");

        // with --annotation it should parse (and not require --seq-tech)
        let ok = Args::try_parse_from([
            "oarfish",
            "--genome-alignments",
            "aln.genome.bam",
            "--annotation",
            "anno.gtf",
            "-o",
            "out",
        ]);
        assert!(ok.is_ok(), "expected parse success, got {ok:?}");
        let ok = ok.expect("valid arguments");
        assert_eq!(
            ok.genome_alignments.as_deref(),
            Some(std::path::Path::new("aln.genome.bam"))
        );
        assert_eq!(
            ok.annotation.as_deref(),
            Some(std::path::Path::new("anno.gtf"))
        );
    }

    #[test]
    fn genome_reads_mode_parses_and_conflicts() {
        // --reads + --genome + --annotation + --seq-tech is the genome-read mode
        let ok = Args::try_parse_from([
            "oarfish",
            "--reads",
            "reads.fq.gz",
            "--genome",
            "genome.fa",
            "--annotation",
            "anno.gtf",
            "--seq-tech",
            "ont-cdna",
            "-o",
            "out",
        ]);
        assert!(ok.is_ok(), "expected parse success, got {ok:?}");

        // --genome conflicts with a transcriptome --index
        let conflict = Args::try_parse_from([
            "oarfish",
            "--reads",
            "reads.fq.gz",
            "--genome",
            "genome.fa",
            "--index",
            "txps.mmi",
            "--annotation",
            "anno.gtf",
            "--seq-tech",
            "ont-cdna",
            "-o",
            "out",
        ]);
        assert!(conflict.is_err(), "expected --genome/--index conflict");
    }
}
