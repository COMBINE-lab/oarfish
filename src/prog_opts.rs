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

/// `--bump-alloc[=<size>]`: the largest block, in bytes, that the per-thread arenas serve,
/// and 0 (`off`, `none`) for no arenas. Without the option, `--mapper spearfish` has the
/// arenas on with the 8 KB limit.
fn parse_bump_alloc(arg: &str) -> anyhow::Result<usize> {
    if arg.eq_ignore_ascii_case("off") || arg.eq_ignore_ascii_case("none") {
        return Ok(0);
    }
    Ok(arg.parse::<nudgy::Config>()?.max_block())
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

/// Parse a non-negative, finite `f32` (`--end-fit`, where 0 turns the end fit off).
fn parse_nonneg_f32(arg: &str) -> anyhow::Result<f32> {
    let v: f32 = arg
        .parse()
        .map_err(|_| anyhow::anyhow!("`{}` is not a valid number", arg))?;
    if v >= 0.0 && v.is_finite() {
        Ok(v)
    } else {
        anyhow::bail!("value must be >= 0, but got {}", v)
    }
}

/// How spearfish's end fit measures where a read ends on each candidate transcript
/// (`--end-fit-from`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum EndFitFrom {
    /// The read taken to continue along the transcript from its last seed, without gaps (no
    /// alignment is run).
    Chain,
    /// The alignment of the read's right end.
    Alignment,
}

/// What spearfish does with its orientation evidence for antisense pairs (`--orient-evidence`):
/// the exact comparison of read orientations whose best hits are one sequence on opposite
/// strands, and the longer allowance after a poly(A) tail for cDNA reads that kept their
/// adapters (spearfish's `orient_same_exact` and `orient_tail_after`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum OrientEvidence {
    /// The evidence moves only the quantification scores, which alignments are weighed by;
    /// the mapping scores, the orientation and the hits reported are as without it
    /// (spearfish's `orient_quant`, on in its 0.3.0 presets).
    Quant,
    /// The evidence moves the mapping scores and decides the orientation, and so which hits
    /// are reported (spearfish's `--no-orient-quant`).
    Mapping,
    /// No orientation evidence: reads map and are weighed as with spearfish 0.2.1 (spearfish's
    /// `--no-orient-same-exact --orient-tail-after=80`).
    Off,
}

/// The default of `--orient-evidence`. `None` keeps the spearfish preset's setting, which is
/// `quant` in spearfish 0.3.0. To map as spearfish 0.2.1 did by default, set this to
/// `Some(OrientEvidence::Off)` (and update the default named in the help of
/// `--orient-evidence` and in the docs).
pub const ORIENT_EVIDENCE_DEFAULT: Option<OrientEvidence> = None;

/// spearfish's `orient_tail_after` with the orientation evidence off: the allowance any A-run
/// has after it (spearfish's `TAIL_AFTER`, the value of its 0.2.1 presets and of its `-fast`
/// presets).
pub const ORIENT_TAIL_AFTER_OFF: u32 = 80;

impl OrientEvidence {
    /// The setting spearfish's tie-break settings amount to: `off` when neither part of the
    /// evidence is on, otherwise where `orient_quant` sends it.
    pub fn of(orient_same_exact: bool, orient_tail_after: u32, orient_quant: bool) -> Self {
        if !orient_same_exact && orient_tail_after <= ORIENT_TAIL_AFTER_OFF {
            OrientEvidence::Off
        } else if orient_quant {
            OrientEvidence::Quant
        } else {
            OrientEvidence::Mapping
        }
    }
}

/// Default denominator `D` of the score→probability conversion `exp((score - best) / D)`
/// (`--score-prob-denom`). 3 rather than the earlier 5: sharper weighting toward the
/// best-scoring alignment improved accuracy on simulations with exact truth and on LongBench
/// (see the CHANGELOG); `--score-prob-denom 5` restores the old behaviour.
pub const DEFAULT_SCORE_PROB_DENOM: f32 = 3.0;

/// Default `--max-em-iter`.
pub const DEFAULT_MAX_EM_ITER: u32 = 1000;
/// Default `--max-em-iter` with `--model-truncation`. With its convergence test (see
/// `crate::util::constants::TRUNC_CONV_MIN_READS`) the truncation EM typically converges after
/// 1,000-1,300 iterations, just past [`DEFAULT_MAX_EM_ITER`]; iterations on its component
/// layout are cheap.
pub const DEFAULT_MAX_EM_ITER_TRUNCATION: u32 = 5000;

/// Model of where reads land on transcripts (`--coverage-model`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, Serialize)]
pub enum CoverageModel {
    /// No positional term.
    #[default]
    None,
    /// The coverage model (`--model-coverage`).
    Logistic,
    /// The learned placement mixture (`--model-truncation`).
    Truncation,
}

impl Args {
    /// Resolves `--coverage-model` and its shorthands into the model in effect, and sets the
    /// shorthand flags to match, which is what the rest of the program reads. Contradictory
    /// choices, and models a mode does not support, are errors.
    pub fn resolve_coverage_model(&mut self) -> anyhow::Result<CoverageModel> {
        let implied = if self.model_coverage {
            Some(CoverageModel::Logistic)
        } else if self.model_truncation {
            Some(CoverageModel::Truncation)
        } else {
            None
        };
        let model = match (self.coverage_model, implied) {
            (Some(m), Some(i)) if m != i => anyhow::bail!(
                "--coverage-model {} contradicts --model-{}",
                clap::ValueEnum::to_possible_value(&m)
                    .expect("named")
                    .get_name(),
                if i == CoverageModel::Logistic {
                    "coverage"
                } else {
                    "truncation"
                }
            ),
            (Some(m), _) => m,
            (None, Some(i)) => i,
            (None, None) => CoverageModel::None,
        };
        if self.single_cell && model == CoverageModel::Truncation {
            anyhow::bail!("the truncation model is not supported in single-cell mode");
        }
        self.coverage_model = Some(model);
        self.model_coverage = model == CoverageModel::Logistic;
        self.model_truncation = model == CoverageModel::Truncation;
        Ok(model)
    }

    /// The EM's iteration limit: `--max-em-iter` if given, else [`DEFAULT_MAX_EM_ITER`], or
    /// [`DEFAULT_MAX_EM_ITER_TRUNCATION`] with the truncation model (after
    /// [`Args::resolve_coverage_model`]).
    pub fn em_iter_limit(&self) -> u32 {
        self.max_em_iter.unwrap_or(if self.model_truncation {
            DEFAULT_MAX_EM_ITER_TRUNCATION
        } else {
            DEFAULT_MAX_EM_ITER
        })
    }
}

/// Presence/absence component for the bulk EM.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, Serialize)]
pub enum PresenceModel {
    /// No presence modelling.
    #[default]
    None,
    /// Spike-and-slab: per-transcript Bernoulli presence posterior with a fixed prior; the EM
    /// sees presence-discounted abundances (zero-inflated EM).
    SpikeSlab,
    /// Spike-and-slab with a marginal-likelihood Bayes factor per transcript and a learned
    /// base rate (see `crate::presence`).
    Bayes,
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, Serialize)]
pub enum MapperBackend {
    Rammap,
    Spearfish,
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
    /// used to weight a read's alignments in the EM (default 3). spearfish
    /// mappings are scored in tenths of a spearfish score point, and the default
    /// is the scale spearfish suggests for them: 20 with `--mapper spearfish`,
    /// or the one a RAD file records. Larger `D` flattens the weighting across
    /// alignments of differing score; smaller `D` sharpens it toward the
    /// best-scoring alignment. Transcriptome mode only — in genome mode projected
    /// alignments are weighted by bramble similarity, so passing this with
    /// `--genome`/`--genome-alignments` is an error.
    #[arg(long, value_parser = parse_pos_f32, help_heading = "filters")]
    pub score_prob_denom: Option<f32>,

    /// per-base match score of the aligner, used to remove the score credit a read's poly(A)
    /// tail earns by aligning to a transcript whose sequence ends in an A-run. By default it is
    /// taken from the aligner (a minimap2 `@PG` command line, or the rammap preset); with an
    /// unrecognized aligner the adjustment is off. 0 turns it off.
    #[arg(long, hide = true)]
    pub tail_match_score: Option<i32>,

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

    /// path where the index will be written (if provided): a rammap index, or with
    /// `--mapper spearfish` a spearfish `.sfi` (with its build manifest as `<path>.json`)
    #[arg(long, conflicts_with_all = ["alignments", "index"], requires_ifs([(ArgPredicate::IsPresent, "only_index")]), help_heading = "indexing")]
    pub index_out: Option<PathBuf>,

    /// with `--mapper spearfish`, when oarfish builds the index from `--annotated`/`--novel`:
    /// a transcript-to-gene table (whitespace-separated `transcript gene` per line; further
    /// columns are ignored) whose genes replace those in the FASTA headers. spearfish groups
    /// isoforms by gene; without gene labels every transcript is its own gene. Transcripts
    /// the table lacks keep the gene from their header
    #[arg(
        long,
        value_name = "FILE",
        conflicts_with_all = ["alignments", "index", "infer_genes"],
        help_heading = "indexing"
    )]
    pub t2g: Option<PathBuf>,

    /// with `--mapper spearfish`, when oarfish builds the index from `--annotated`/`--novel`:
    /// group transcripts into genes from shared sequence (spearfish's gene inference, with its
    /// default settings), for transcriptomes whose FASTA headers carry no gene IDs and that
    /// have no transcript-to-gene table
    #[arg(
        long,
        conflicts_with_all = ["alignments", "index"],
        help_heading = "indexing"
    )]
    pub infer_genes: bool,

    /// sequencing technology in which to expect reads if using mapping based mode
    #[arg(
        long,
        help_heading = "raw read mode",
        required_unless_present_any = ["alignments", "genome_alignments"],
        value_parser = clap::value_parser!(SequencingTech)
    )]
    pub seq_tech: Option<SequencingTech>,

    /// read-mode mapping backend for transcriptome quantification: `rammap` (default) or
    /// `spearfish`. spearfish builds its own index from `--annotated`/`--novel` (and writes it
    /// with `--index-out`), or opens a spearfish `.sfi` given as `--index`; its preset follows
    /// `--seq-tech`.
    #[arg(long, value_enum, default_value_t = MapperBackend::Rammap, help_heading = "raw read mode")]
    pub mapper: MapperBackend,

    /// with `--mapper spearfish`, which orientation of each read to map: `fw` (the read as
    /// given), `rc` (its reverse complement) or `both`. The default is the preset's: `both`,
    /// and `fw` for `ont-drna`. For reads already oriented to the transcript strand (e.g.
    /// Iso-Seq FLNC reads), `fw` maps the one orientation, which takes less time. This is not
    /// a filter (see `--strand-filter`): a read of the other strand can still have a weak
    /// hit in the mapped orientation, so use it only for oriented reads.
    #[arg(long, value_parser = parse_strand, requires = "reads", help_heading = "raw read mode")]
    pub map_strand: Option<bio_types::strand::Strand>,

    /// with `--mapper spearfish`, the mapper's small allocations (up to 8 KB) are served from
    /// memory that belongs to each thread instead of from the general allocator. Mapping a
    /// read makes a few hundred of them and frees nearly all before the next read. A little
    /// less mapping time, for a few MB of memory a thread; what is mapped does not change.
    /// This is on by default: `--bump-alloc=off` (or `none`) turns the arenas off, and a
    /// size sets the largest allocation served this way in place of 8 KB
    /// (`--bump-alloc=1KB`, at most 32KB).
    #[arg(
        long,
        value_name = "SIZE",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "8KB",
        value_parser = parse_bump_alloc,
        requires = "reads",
        help_heading = "raw read mode"
    )]
    pub bump_alloc: Option<usize>,

    /// with `--mapper spearfish`, the end fit's weight: points added to the score of the hits
    /// whose annotated 3' end a read without a poly(A) tail reaches, when the read's candidate
    /// transcripts end at distinguishable places. It only changes the scores alignments are
    /// weighed by, not which hits a read has. 0 turns it off; alignments are then weighed by
    /// their mapping scores, but for the reads spearfish's orientation evidence moves (see
    /// `--orient-evidence`; `--orient-evidence off` too gives mapping-score weights). The
    /// default is the preset's (5).
    #[arg(
        long,
        value_name = "SCORE",
        value_parser = parse_nonneg_f32,
        requires = "reads",
        help_heading = "raw read mode"
    )]
    pub end_fit: Option<f32>,

    /// with `--mapper spearfish`, the end fit's tolerance: how far (in homopolymer-compressed
    /// bases) where a read ends may lie from a transcript's annotated 3' end and still fit
    /// it. Transcripts whose 3' ends are closer to each other than this are not told apart.
    /// The default is the preset's (4).
    #[arg(
        long,
        value_name = "BASES",
        requires = "reads",
        help_heading = "raw read mode"
    )]
    pub end_fit_tol: Option<u32>,

    /// with `--mapper spearfish`, how the end fit measures where a read ends on each
    /// candidate transcript: `chain` (the read taken to continue along the transcript from
    /// its last seed, without gaps; no alignment is run) or `alignment` (the alignment of the
    /// read's right end; more exact, and more mapping time where spans are not extended).
    /// The default is the preset's: `chain` for `ont-drna` and `alignment` otherwise, but
    /// `alignment` for every preset with `--model-truncation`, which aligns the read ends
    /// anyway. Hits whose spans were extended are measured from their alignment either way.
    #[arg(long, value_enum, requires = "reads", help_heading = "raw read mode")]
    pub end_fit_from: Option<EndFitFrom>,

    /// with `--mapper spearfish`, what spearfish's orientation evidence for antisense pairs
    /// does. Where a read's best hits in its two orientations are one sequence on opposite
    /// strands (antisense isoforms, a transcript inside another's intron on the other strand),
    /// spearfish compares the two by exact alignment of the whole read, and gives cDNA reads
    /// that kept their adapters a longer allowance after a poly(A) tail. `quant`: the evidence
    /// moves only the quantification scores alignments are weighed by; which hits a read has
    /// does not change. `mapping`: it moves the mapping scores and decides the orientation.
    /// `off`: no orientation evidence, as spearfish 0.2.1 mapped; with `--end-fit 0` as well,
    /// alignments are weighed by their mapping scores alone. The default is the preset's
    /// (`quant`).
    #[arg(long, value_enum, requires = "reads", help_heading = "raw read mode")]
    pub orient_evidence: Option<OrientEvidence>,

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

    /// fraction of a query that must be mapped within an alignment to consider the alignment
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

    /// which model of where reads land on transcripts to use: `none`, `logistic` (the coverage
    /// model) or `truncation` (a learned mixture of 3'-anchored, 5'-anchored, full-length and
    /// uniformly placed reads); default `none`. `--model-coverage` and `--model-truncation` are
    /// shorthands for the last two.
    #[arg(long, help_heading = "coverage model", value_enum)]
    pub coverage_model: Option<CoverageModel>,

    /// apply the coverage model (same as `--coverage-model logistic`)
    #[arg(long, help_heading = "coverage model", value_parser)]
    pub model_coverage: bool,

    /// model where reads land on transcripts: a mixture of 3'-anchored, 5'-anchored, full-length
    /// and uniformly placed reads, whose weights, read-length laws and end-gap laws are learned
    /// from the data in the EM; the fitted mixture is written to the log (same as
    /// `--coverage-model truncation`). With the spearfish backend (`--mapper spearfish`), this
    /// also extends alignment spans to the read ends where the preset leaves them at the
    /// outermost seeds (`--seq-tech ont-drna` and `pac-bio-hifi`). Longer spans change which reads
    /// `--min-aligned-fraction` drops, so the reads kept differ from a run without
    /// `--model-truncation`; the log and `meta_info.json` (`spearfish_settings`) record
    /// whether spans were extended and why
    #[arg(
        long,
        help_heading = "coverage model",
        conflicts_with = "model_coverage"
    )]
    pub model_truncation: bool,

    /// if using the coverage model, use this as the value of `k` in the logistic equation
    #[arg(
        short = 'k',
        long,
        help_heading = "coverage model",
        value_parser,
        default_value_t = 2.0
    )]
    pub growth_rate: f64,

    /// write output alignment probabilities (optionally compressed) for each mapped read.
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

    /// maximum number of iterations for which to run the EM algorithm (default 1000; 5000 with
    /// `--model-truncation`). A run that reaches the limit before converging gets a warning in
    /// the log and in `meta_info.json`
    #[arg(long, help_heading = "EM")]
    pub max_em_iter: Option<u32>,

    /// the EM stops once the largest relative change in a transcript's abundance between two
    /// rounds is below this value
    #[arg(long, help_heading = "EM", default_value_t = 1e-3)]
    pub convergence_thresh: f64,

    /// run the EM over individual reads even when no per-read model term (coverage, KDE) is
    /// active; by default such reads are collapsed into equivalence classes
    #[arg(long, help_heading = "EM")]
    pub per_read_em: bool,

    /// presence/absence model (experimental): each transcript carries a Bernoulli presence
    /// posterior that discounts its effective abundance during the EM (zero-inflated EM). It
    /// suppresses low-evidence false-positive transcripts at some cost to true low-abundance
    /// ones. Supported by the equivalence-class EM and `--model-truncation`.
    #[arg(long, help_heading = "EM", value_enum, default_value_t = PresenceModel::None)]
    pub presence_model: PresenceModel,

    /// presence prior probability, a fixed penalty: roughly `-logit(rho)` effective exclusive
    /// reads are needed for a transcript to be called present
    #[arg(long, hide = true, default_value_t = 0.05)]
    pub presence_rho: f64,

    /// EM iterations before the first presence update
    #[arg(long, hide = true, default_value_t = 50)]
    pub presence_warmup: u32,

    /// EM iterations between presence updates
    #[arg(long, hide = true, default_value_t = 10)]
    pub presence_period: u32,

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

    /// The spearfish settings the alignments were made with: span extension and why, the
    /// orientations mapped and the quantification score's settings (set by the spearfish
    /// backend and by RAD input; recorded in `meta_info.json` as `spearfish_settings`).
    #[arg(skip)]
    #[serde(skip)]
    pub spearfish_settings: Option<SpearfishSettings>,
}

/// The spearfish settings a run's alignments were made with, as `meta_info.json` records them
/// (`spearfish_settings`). In process they are resolved from the preset and oarfish's options;
/// for a RAD file they are read from its header, and a setting the header does not carry is
/// `null`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SpearfishSettings {
    /// `in_process` (oarfish mapped the reads) or `rad_header` (read from a RAD file).
    pub source: &'static str,
    /// The mapper and its version, as the RAD file records them (RAD input only).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mapper: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mapper_version: Option<String>,
    /// The mapping preset.
    pub preset: String,
    /// Whether alignment spans were extended to the read ends.
    pub extend_ends: bool,
    /// Why (see [`SpanExtension::reason`]).
    pub extend_ends_reason: &'static str,
    /// The orientations mapped: `both`, `fwd` or `rev`.
    pub strands: String,
    /// Whether the alignments carry quantification scores that may differ from their mapping
    /// scores (otherwise the two are equal).
    pub quant_scores: bool,
    /// The end fit's points (0: off), tolerance and measurement.
    pub end_fit: Option<f32>,
    pub end_fit_tol: Option<u32>,
    pub end_fit_from: Option<EndFitFrom>,
    /// Whether spearfish's orientation evidence acts through the quantification score.
    pub orient_quant: Option<bool>,
    /// What the orientation evidence does: `quant`, `mapping` or `off` (see
    /// [`OrientEvidence`]).
    pub orient_evidence: Option<OrientEvidence>,
    /// Where the index's gene labels came from, and whether it groups isoforms into genes.
    pub index_genes: Option<IndexGenes>,
}

/// The gene labels of the spearfish index the alignments were made against, as
/// `meta_info.json` records them (`spearfish_settings.index_genes`).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct IndexGenes {
    /// Where the genes came from, as the index manifest's `genes_from` (or, for RAD input, the
    /// file tag `index_genes_from`): `fasta` (gene IDs in the FASTA headers) or `t2g` (a
    /// transcript-to-gene table), both provided; `provided` (genes given through spearfish's
    /// Python API, `given` in the RAD file); `inferred` (spearfish's gene inference); `none`
    /// (every transcript its own gene); or `unknown` (an index without a manifest, or a RAD
    /// file written before spearfish 0.3.1, which does not record it).
    pub source: String,
    /// Genes in the index (`null` when unknown).
    pub genes: Option<usize>,
    /// Transcripts in the index.
    pub transcripts: Option<usize>,
    /// Whether the index groups isoforms into genes (spearfish's `Index::has_gene_groups`;
    /// `null` when unknown).
    pub has_gene_groups: Option<bool>,
}

impl IndexGenes {
    /// The record for an index whose gene labels are not known (`source: unknown`).
    pub fn unknown(transcripts: Option<usize>) -> Self {
        IndexGenes {
            source: "unknown".to_string(),
            genes: None,
            transcripts,
            has_gene_groups: None,
        }
    }
}

impl SpearfishSettings {
    /// The settings from a RAD file's header. The end fit is read from the mapper's parameters
    /// (spearfish records its options there as JSON); it is `null` where they do not have it.
    pub fn from_rad(m: &long_read_rad::Meta) -> Self {
        let params: serde_json::Value = serde_json::from_str(&m.params).unwrap_or_default();
        let tb = &params["parameters"]["params"]["tiebreak"];
        SpearfishSettings {
            source: "rad_header",
            mapper: Some(m.mapper.clone()),
            mapper_version: Some(m.mapper_version.clone()),
            preset: m.preset.clone(),
            extend_ends: m.ends_extended,
            extend_ends_reason: "rad_header",
            strands: m.strands.clone(),
            quant_scores: m.quant_scores,
            end_fit: tb["end_fit"].as_f64().map(|v| v as f32),
            end_fit_tol: tb["end_fit_tol"].as_u64().map(|v| v as u32),
            end_fit_from: tb["end_fit_chain"].as_bool().map(|chain| {
                if chain {
                    EndFitFrom::Chain
                } else {
                    EndFitFrom::Alignment
                }
            }),
            orient_quant: tb["orient_quant"].as_bool(),
            index_genes: None,
            orient_evidence: match (
                tb["orient_same_exact"].as_bool(),
                tb["orient_tail_after"].as_u64(),
                tb["orient_quant"].as_bool(),
            ) {
                (Some(exact), Some(after), Some(quant)) => {
                    Some(OrientEvidence::of(exact, after as u32, quant))
                }
                _ => None,
            },
        }
    }

    /// Whether alignment spans were extended, and why, for the log line.
    pub fn span_extension(&self) -> SpanExtension {
        SpanExtension {
            extended: self.extend_ends,
            reason: self.extend_ends_reason,
            preset: self.preset.clone(),
        }
    }
}

/// Whether alignment spans reach the read ends, or stop at the outermost seeds, and why. Span
/// length decides the aligned fraction of a read, and so which reads `--min-aligned-fraction`
/// drops. Recorded in `meta_info.json` as part of [`SpearfishSettings`].
#[derive(Debug, Clone, PartialEq)]
pub struct SpanExtension {
    /// Whether alignment spans were extended to the read ends.
    pub extended: bool,
    /// Why: `preset` (the mapping preset extends them), `model_truncation` (the preset leaves
    /// them at the outermost seeds and `--model-truncation` turned extension on), `preset_off`
    /// (the preset leaves them at the outermost seeds), or `rad_header` (RAD input: as the
    /// mapper recorded in the file).
    pub reason: &'static str,
    /// The mapper's preset.
    pub preset: String,
}

impl SpanExtension {
    /// The spearfish backend's setting: extended if the preset extends spans (`preset_extends`)
    /// or if the truncation model is on.
    pub fn in_process(preset: String, preset_extends: bool, model_truncation: bool) -> Self {
        let (extended, reason) = match (preset_extends, model_truncation) {
            (true, _) => (true, "preset"),
            (false, true) => (true, "model_truncation"),
            (false, false) => (false, "preset_off"),
        };
        SpanExtension {
            extended,
            reason,
            preset,
        }
    }

    /// One line for the log.
    pub fn describe(&self) -> String {
        match self.reason {
            "preset" => format!(
                "alignment spans extended to the read ends (as preset {} does)",
                self.preset
            ),
            "model_truncation" => format!(
                "alignment spans extended to the read ends because of --model-truncation \
                 (preset {} leaves them at the outermost seeds). Longer spans change which \
                 reads --min-aligned-fraction drops, so the reads kept differ from a run \
                 without --model-truncation",
                self.preset
            ),
            "rad_header" if self.extended => format!(
                "alignment spans extended to the read ends by the mapper (preset {}, from the \
                 RAD header)",
                self.preset
            ),
            "rad_header" => format!(
                "alignment spans end at the outermost seeds (preset {}, from the RAD header)",
                self.preset
            ),
            _ => format!(
                "alignment spans end at the outermost seeds (preset {}; --model-truncation \
                 would extend them)",
                self.preset
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Args, DEFAULT_MAX_EM_ITER, DEFAULT_MAX_EM_ITER_TRUNCATION};
    use clap::Parser;

    /// `--coverage-model` and its shorthands resolve to one model; contradictions and the
    /// truncation model in single-cell mode are errors.
    #[test]
    fn coverage_model_selector_and_shorthands_agree() {
        use super::CoverageModel as M;
        let resolve = |extra: &[&str]| {
            let mut v = vec![
                "oarfish",
                "-a",
                "x.bam",
                "--seq-tech",
                "ont-cdna",
                "-o",
                "out",
            ];
            v.extend_from_slice(extra);
            let mut a = Args::try_parse_from(v).expect("valid arguments");
            a.resolve_coverage_model()
                .map(|m| (m, a.model_coverage, a.model_truncation))
        };
        assert_eq!(resolve(&[]).unwrap(), (M::None, false, false));
        assert_eq!(
            resolve(&["--model-coverage"]).unwrap(),
            (M::Logistic, true, false)
        );
        assert_eq!(
            resolve(&["--model-truncation"]).unwrap(),
            (M::Truncation, false, true)
        );
        assert_eq!(
            resolve(&["--coverage-model", "logistic"]).unwrap(),
            (M::Logistic, true, false)
        );
        assert_eq!(
            resolve(&["--coverage-model", "truncation"]).unwrap(),
            (M::Truncation, false, true)
        );
        assert_eq!(
            resolve(&["--coverage-model", "truncation", "--model-truncation"]).unwrap(),
            (M::Truncation, false, true)
        );
        assert!(resolve(&["--coverage-model", "none", "--model-coverage"]).is_err());
        assert!(resolve(&["--coverage-model", "logistic", "--model-truncation"]).is_err());
        assert!(resolve(&["--single-cell", "--coverage-model", "truncation"]).is_err());
        assert_eq!(
            resolve(&["--single-cell", "--model-coverage"]).unwrap(),
            (M::Logistic, true, false)
        );
    }

    #[test]
    fn span_extension_records_why() {
        use super::SpanExtension;
        let s = SpanExtension::in_process("ont".into(), true, false);
        assert_eq!((s.extended, s.reason), (true, "preset"));
        let s = SpanExtension::in_process("ont".into(), true, true);
        assert_eq!((s.extended, s.reason), (true, "preset"));
        let s = SpanExtension::in_process("drna".into(), false, true);
        assert_eq!((s.extended, s.reason), (true, "model_truncation"));
        assert!(s.describe().contains("--min-aligned-fraction"));
        let s = SpanExtension::in_process("hifi".into(), false, false);
        assert_eq!((s.extended, s.reason), (false, "preset_off"));
    }

    /// A RAD file's settings come from its header, with the end fit from the mapper's
    /// parameters where they have it.
    #[test]
    fn spearfish_settings_from_a_rad_header() {
        let mut m = long_read_rad::Meta {
            mapper: "spearfish".into(),
            mapper_version: "0.3.0".into(),
            preset: "drna".into(),
            strands: "fwd".into(),
            quant_scores: true,
            params: r#"{"parameters":{"params":{"tiebreak":{"end_fit":5.0,"end_fit_tol":4,
                "end_fit_chain":true,"orient_quant":true,"orient_same_exact":false,
                "orient_tail_after":80}}}}"#
                .into(),
            ..Default::default()
        };
        let v = serde_json::to_value(super::SpearfishSettings::from_rad(&m)).unwrap();
        assert_eq!(v["source"], "rad_header");
        assert_eq!(v["mapper_version"], "0.3.0");
        assert_eq!(v["extend_ends"], false);
        assert_eq!(v["extend_ends_reason"], "rad_header");
        assert_eq!(v["quant_scores"], true);
        assert_eq!(v["end_fit"], 5.0);
        assert_eq!(v["end_fit_tol"], 4);
        assert_eq!(v["end_fit_from"], "chain");
        assert_eq!(v["orient_quant"], true);
        assert_eq!(v["orient_evidence"], "off");
        // A file from another mapper, or an older one, has no end fit in its parameters.
        m.params = "not json".into();
        m.quant_scores = false;
        m.ends_extended = true;
        let s = super::SpearfishSettings::from_rad(&m);
        assert!(s.end_fit.is_none() && s.end_fit_from.is_none() && !s.quant_scores);
        assert!(s.span_extension().describe().contains("extended"));
    }

    #[test]
    fn em_iter_limit_depends_on_the_model_unless_given() {
        let limit = |extra: &[&str]| {
            let mut v = vec![
                "oarfish",
                "-a",
                "x.bam",
                "--seq-tech",
                "ont-cdna",
                "-o",
                "out",
            ];
            v.extend_from_slice(extra);
            let mut a = Args::try_parse_from(v).expect("valid arguments");
            a.resolve_coverage_model().expect("a valid model");
            a.em_iter_limit()
        };
        assert_eq!(limit(&[]), DEFAULT_MAX_EM_ITER);
        assert_eq!(limit(&["--model-coverage"]), DEFAULT_MAX_EM_ITER);
        assert_eq!(
            limit(&["--model-truncation"]),
            DEFAULT_MAX_EM_ITER_TRUNCATION
        );
        assert_eq!(
            limit(&["--coverage-model", "truncation"]),
            DEFAULT_MAX_EM_ITER_TRUNCATION
        );
        assert_eq!(limit(&["--max-em-iter", "700"]), 700);
        assert_eq!(limit(&["--model-truncation", "--max-em-iter", "700"]), 700);
    }

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

    /// `--bump-alloc` is unset without the option (the spearfish backend then turns the
    /// arenas on with their default limit), 8 KB without a value, and takes a size or
    /// `off` / `none` joined by `=`; a following argument is never its value.
    #[test]
    fn bump_alloc_takes_an_optional_size() {
        let parse = |extra: &[&str]| {
            let mut v = vec!["oarfish", "--mapper", "spearfish"];
            v.extend_from_slice(extra);
            v.extend_from_slice(&["--reads", "r.fq", "--index", "t.sfi", "-o", "out"]);
            v.extend_from_slice(&["--seq-tech", "ont-cdna"]);
            Args::try_parse_from(v).map(|a| a.bump_alloc)
        };
        assert_eq!(parse(&[]).unwrap(), None);
        assert_eq!(parse(&["--bump-alloc"]).unwrap(), Some(8192));
        assert_eq!(parse(&["--bump-alloc=8KB"]).unwrap(), Some(8192));
        assert_eq!(parse(&["--bump-alloc=1KB"]).unwrap(), Some(1024));
        assert_eq!(parse(&["--bump-alloc=500"]).unwrap(), Some(512));
        assert_eq!(parse(&["--bump-alloc=32kb"]).unwrap(), Some(32768));
        for off in ["--bump-alloc=off", "--bump-alloc=none", "--bump-alloc=0"] {
            assert_eq!(parse(&[off]).unwrap(), Some(0), "{off}");
        }
        for bad in [
            "--bump-alloc=64KB",
            "--bump-alloc=8MB",
            "--bump-alloc=on",
            "--bump-alloc=",
        ] {
            assert!(parse(&[bad]).is_err(), "{bad}");
        }
        // The option needs reads to map, and its default does not.
        let no_reads = |extra: &[&str]| {
            let mut v = vec!["oarfish", "-a", "x.bam", "-o", "out"];
            v.extend_from_slice(extra);
            Args::try_parse_from(v).map(|a| a.bump_alloc)
        };
        assert_eq!(no_reads(&[]).unwrap(), None);
        assert!(no_reads(&["--bump-alloc"]).is_err());
    }
}
