use clap::Parser;
use serde::Serialize;
use std::path::PathBuf;

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

/// accurate transcript quantification from long-read RNA-seq data
#[derive(Parser, Debug, Serialize)]
#[clap(author, version, about, long_about = None)]
pub struct Args {
    /// be quiet (i.e. don't output log messages that aren't at least warnings)
    #[arg(long, conflicts_with = "verbose")]
    pub quiet: bool,

    /// be verbose (i.e. output all non-developer logging messages)
    #[arg(long)]
    pub verbose: bool,

    /// path to the file containing the input alignments
    #[arg(short, long, required = true)]
    pub alignments: PathBuf,
    /// location where output quantification file should be written
    #[arg(short, long, required = true)]
    pub output: PathBuf,

    #[arg(long, help_heading = "filters", value_enum)]
    pub filter_group: Option<FilterGroup>,

    /// maximum allowable distance of the right-most end of an alignment from the 3' transcript end
    #[arg(short, long, conflicts_with = "filter-group", help_heading="filters", default_value_t = u32::MAX as i64)]
    pub three_prime_clip: i64,

    /// maximum allowable distance of the left-most end of an alignment from the 5' transcript end
    #[arg(short, long, conflicts_with = "filter-group", help_heading="filters", default_value_t = u32::MAX)]
    pub five_prime_clip: u32,

    /// fraction of the best possible alignment score that a secondary alignment must have for
    /// consideration
    #[arg(
        short,
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = 0.95
    )]
    pub score_threshold: f32,

    /// fraction of a query that must be mapped within an alignemnt to consider the alignemnt
    /// valid
    #[arg(
        short,
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = 0.5
    )]
    pub min_aligned_fraction: f32,

    /// minimum number of nucleotides in the aligned portion of a read
    #[arg(
        short = 'l',
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = 50
    )]
    pub min_aligned_len: u32,

    /// only alignments to this strand will be allowed; options are (fw /+, rc/-, or both/.)
    #[arg(
        short = 'd',
        long,
        conflicts_with = "filter-group",
        help_heading = "filters",
        default_value_t = bio_types::strand::Strand::Unknown,
        value_parser = parse_strand
    )]
    pub strand_filter: bio_types::strand::Strand,

    /// input is assumed to be a single-cell BAM and to have the `CB:z` tag for all read records
    #[arg(long)]
    pub single_cell: bool,

    /// apply the coverage model
    #[arg(long, help_heading = "coverage model", value_parser)]
    pub model_coverage: bool,

    /// maximum number of iterations for which to run the EM algorithm
    #[arg(long, help_heading = "EM", default_value_t = 1000)]
    pub max_em_iter: u32,

    /// maximum number of iterations for which to run the EM algorithm
    #[arg(long, help_heading = "EM", default_value_t = 1e-3)]
    pub convergence_thresh: f64,

    /// maximum number of cores that the oarfish can use to obtain binomial probability
    #[arg(short = 'j', long, default_value_t = 1)]
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

    /// use a KDE model of the observed fragment length distribution
    #[arg(short, long, hide = true)]
    pub use_kde: bool,
}
