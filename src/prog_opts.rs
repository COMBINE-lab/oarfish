use clap::{builder::ArgPredicate, Parser};
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

#[derive(Debug, Clone, clap::ValueEnum, Serialize)]
pub enum ReadAssignmentProbOut {
    Uncompressed,
    Compressed,
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
    .args(["alignments", "reads"])
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

    /// path to the file containing the input reads
    #[arg(
        long,
        help_heading = "raw read mode",
        value_delimiter = ',',
        requires_ifs([
            (ArgPredicate::IsPresent, "reference"),
            (ArgPredicate::IsPresent, "seq_tech")
        ])
    )]
    pub reads: Option<Vec<PathBuf>>,

    /// path to the file containing the reference transcriptome (or existing index) against which
    /// to map
    #[arg(long, conflicts_with = "alignments", help_heading = "raw read mode")]
    pub reference: Option<PathBuf>,

    /// path where minimap2 index will be written (if provided)
    #[arg(long, conflicts_with = "alignments", help_heading = "raw read mode")]
    pub index_out: Option<PathBuf>,

    /// sequencing technology in which to expect reads if using mapping based mode
    #[arg(
        long,
        help_heading = "raw read mode",
        required_unless_present = "alignments",
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

    /// location where output quantification file should be written
    #[arg(short, long, required = true)]
    pub output: PathBuf,

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

    /// input is assumed to be a single-cell BAM and to have the `CB:z` tag for all read records
    #[arg(long, conflicts_with = "reads")]
    pub single_cell: bool,

    /// apply the coverage model
    #[arg(long, help_heading = "coverage model", value_parser)]
    pub model_coverage: bool,

    /// write output alignment probabilites (optionally compressed) for each mapped read.
    /// If <WRITE_ASSIGNMENT_PROBS> is present, it must be one of `uncompressed` (default) or
    /// `compressed`, which will cause the output file to be lz4 compressed.
    #[arg(
        long,
        help_heading = "output read-txps probabilities",
        conflicts_with = "single-cell",
        default_missing_value = "uncompressed",
        num_args = 0..=1,
        require_equals = true,
        value_parser = parse_assign_prob_out_value
    )]
    pub write_assignment_probs: Option<ReadAssignmentProbOut>,

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
