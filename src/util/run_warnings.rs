//! Warnings about a run as a whole: conditions that do not stop quantification, but that the
//! user should look at before trusting the result. Each is written to the log and recorded
//! in the `warnings` array of `meta_info.json`.
//!
//! All of them are computed from what the run keeps anyway (the discard table, and how the EM
//! stopped), so they cost nothing per read.

use bio_types::strand::Strand;
use serde::Serialize;

use crate::util::oarfish_types::DiscardTable;

/// The share of mapped reads the aligned-fraction filter may drop before a warning.
pub const MAX_ALIGNED_FRACTION_DROP: f64 = 0.10;
/// The share of mapped reads the strand filter may leave without an alignment before a warning.
pub const MAX_STRAND_FILTER_DROP: f64 = 0.10;
/// The share of mapped reads the other alignment filters may leave without an alignment
/// before a warning.
pub const MAX_ALIGNMENT_FILTER_DROP: f64 = 0.10;
/// The share of reads that may stay unmapped under `--map-strand fw|rc` before a warning.
pub const MAX_UNMAPPED_ONE_STRAND: f64 = 0.30;
/// The share of reads that may stay unmapped before a warning.
pub const MAX_UNMAPPED: f64 = 0.50;

/// One warning about the run.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RunWarning {
    /// A stable name for the condition, for programs that read `meta_info.json`.
    pub code: &'static str,
    /// The share of reads that triggered the warning (for `em_not_converged` and
    /// `truncation_em_not_converged`, the EM's last convergence statistic; for
    /// `index_without_gene_groups`, the index's genes per transcript).
    pub fraction: f64,
    /// The share above which the warning is issued (for the two EM codes, the convergence
    /// threshold; for `index_without_gene_groups`, 1).
    pub threshold: f64,
    /// What happened and what to check, for a person.
    pub message: String,
}

/// How an EM run stopped: after `niter` iterations, either because its convergence test held
/// or because it reached the iteration limit `max_iter`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmConvergence {
    /// Iterations run (a SQUAREM cycle counts as three).
    pub niter: u32,
    /// The iteration limit (`--max-em-iter`).
    pub max_iter: u32,
    /// Whether the convergence test held before the limit.
    pub converged: bool,
    /// The convergence statistic of the last iteration (the largest relative change of a
    /// transcript's abundance, as that EM measures it).
    pub rel_diff: f64,
    /// The threshold the statistic had to fall below (`--convergence-thresh`).
    pub thresh: f64,
}

/// Which EM ran: the truncation model's, or one of the others (the equivalence-class EM, the
/// per-read EMs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum EmKind {
    Standard,
    Truncation,
}

/// What `meta_info.json` records about the EM (`em`): which one ran, its iterations and limit,
/// and whether it converged.
pub fn em_info(em: &EmConvergence, kind: EmKind) -> serde_json::Value {
    serde_json::json!({
        "kind": kind,
        "iterations": em.niter,
        "max_iter": em.max_iter,
        "converged": em.converged,
        "last_rel_change": em.rel_diff,
        "convergence_thresh": em.thresh,
    })
}

/// The warning for an EM that stopped at its iteration limit before converging, if it did:
/// `em_not_converged` for the standard EMs, `truncation_em_not_converged` for the truncation
/// EM, each with its own advice.
pub fn check_em(em: &EmConvergence, kind: EmKind) -> Option<RunWarning> {
    let stat = format!(
        "the largest relative change in a transcript's abundance in the last iteration was \
         {:.3e}, above the threshold {:.0e} (--convergence-thresh)",
        em.rel_diff, em.thresh
    );
    (!em.converged).then(|| match kind {
        EmKind::Standard => RunWarning {
            code: "em_not_converged",
            fraction: em.rel_diff,
            threshold: em.thresh,
            message: format!(
                "The EM stopped at its iteration limit ({} iterations, --max-em-iter) before \
                 converging: {stat}. Estimates are usually stable by then; to let the EM \
                 converge, raise --max-em-iter (default 1,000).",
                em.max_iter
            ),
        },
        EmKind::Truncation => RunWarning {
            code: "truncation_em_not_converged",
            fraction: em.rel_diff,
            threshold: em.thresh,
            message: format!(
                "The truncation EM (--model-truncation) stopped at its iteration limit ({} \
                 iterations, --max-em-iter) before converging: {stat}. Raise --max-em-iter \
                 (default 5,000 with --model-truncation) to let it converge.",
                em.max_iter
            ),
        },
    })
}

/// The warning for a spearfish index without gene groups (`index_without_gene_groups`), if the
/// index has none. Its `fraction` is the index's genes per transcript (1: every transcript its
/// own gene) and its `threshold` is 1.
pub fn check_index_genes(genes: &crate::prog_opts::IndexGenes) -> Option<RunWarning> {
    if genes.has_gene_groups != Some(false) {
        return None;
    }
    let cause = match genes.source.as_str() {
        "t2g" => "the transcript-to-gene table puts each transcript in a gene of its own",
        "inferred" => "gene inference grouped no transcripts",
        "fasta" => "the gene IDs in the FASTA headers put each transcript in a gene of its own",
        "provided" => {
            "the gene labels given when the index was built put each transcript in a \
             gene of its own"
        }
        "none" => {
            "the FASTA headers carry no gene IDs, and the index was built without a \
             transcript-to-gene table or gene inference"
        }
        _ => "the index does not record where its genes came from",
    };
    let fix = if genes.source == "inferred" {
        "rebuild the index with a transcript-to-gene table (`--t2g FILE` when oarfish builds \
         it, `spearfish index --t2g FILE` otherwise)"
    } else {
        "rebuild the index with a transcript-to-gene table (`--t2g FILE` when oarfish builds \
         it, `spearfish index --t2g FILE` otherwise), or with gene inference (`--infer-genes`, \
         or `spearfish index --infer-genes`) when no gene annotation exists"
    };
    let (n_genes, n_tx) = (
        genes.genes.unwrap_or(0),
        genes.transcripts.unwrap_or(0).max(1),
    );
    Some(RunWarning {
        code: "index_without_gene_groups",
        fraction: n_genes as f64 / n_tx as f64,
        threshold: 1.0,
        message: format!(
            "The spearfish index has no gene groups: each of its {n_tx} transcripts is its own \
             gene ({cause}). spearfish keeps a few candidate places per gene for each read, so \
             isoforms then compete for those places one by one instead of as a gene, and ties \
             among them go by their order in the FASTA. This costs accuracy. To fix it, \
             {fix}."
        ),
    })
}

/// What the checks need to know about the run besides the discard table.
#[derive(Debug, Clone, Copy)]
pub struct RunContext {
    /// The strand whose alignments the strand filter keeps (`Unknown`: both).
    pub strand_filter: Strand,
    /// The orientation given as `--map-strand`, if any.
    pub map_strand: Option<Strand>,
    /// The minimum aligned fraction in force.
    pub min_aligned_fraction: f32,
}

fn share(part: u64, whole: u64) -> f64 {
    if whole == 0 {
        0.0
    } else {
        part as f64 / whole as f64
    }
}

fn pct(x: f64) -> String {
    format!("{:.1}%", 100.0 * x)
}

/// Whether `s` names one strand. (`Strand::Unknown` does not compare equal to itself, so
/// this cannot be written as `s != Strand::Unknown`.)
fn one_strand(s: Strand) -> bool {
    matches!(s, Strand::Forward | Strand::Reverse)
}

fn strand_name(s: Strand) -> &'static str {
    match s {
        Strand::Forward => "fw",
        Strand::Reverse => "rc",
        Strand::Unknown => "both",
    }
}

/// The warnings the discard table of a finished run gives rise to.
pub fn check(dt: &DiscardTable, ctx: &RunContext) -> Vec<RunWarning> {
    let mut warnings = Vec::new();
    let seen = dt.reads_seen();
    let mapped = dt.reads_mapped();

    // Reads dropped because too little of the read is in its best alignment.
    let low_frac = dt.reads_low_aligned_fraction();
    let f = share(low_frac, mapped);
    if f > MAX_ALIGNED_FRACTION_DROP {
        // under `--map-strand`, reads of the other strand that map in part end up here too
        let one_strand_note = match ctx.map_strand {
            Some(strand) if one_strand(strand) => {
                " With `--map-strand`, reads of the other strand that match only in part are \
                 dropped here too."
            }
            _ => "",
        };
        warnings.push(RunWarning {
            code: "aligned_fraction_filter",
            fraction: f,
            threshold: MAX_ALIGNED_FRACTION_DROP,
            message: format!(
                "{low_frac} of {mapped} mapped reads ({}) were dropped because less than {} of the \
                 read is in its best alignment (--min-aligned-fraction). Reads that still carry \
                 adapters or primers are the usual cause: trim them, keeping the poly(A) tails \
                 (see https://combine-lab.github.io/oarfish/guides/sequencing-technologies/\
                 #when-the-aligned-fraction-filter-drops-many-reads). For gene-level analyses, \
                 `--filter-group no-filters` or a lower `--min-aligned-fraction` keeps these \
                 reads; for transcript-level estimates, the default agreed better with \
                 short-read estimates in our tests.{one_strand_note}",
                pct(f),
                ctx.min_aligned_fraction
            ),
        });
    }

    // Reads all of whose alignments are on the strand the strand filter removes.
    let wrong_strand = dt.reads_filtered_strand();
    let f = share(wrong_strand, mapped);
    if one_strand(ctx.strand_filter) && f > MAX_STRAND_FILTER_DROP {
        warnings.push(RunWarning {
            code: "strand_filter",
            fraction: f,
            threshold: MAX_STRAND_FILTER_DROP,
            message: format!(
                "{wrong_strand} of {mapped} mapped reads ({}) align only to the strand the strand \
                 filter removes (the filter keeps `{}`). Keep one strand only for reads that are \
                 all in that orientation (direct RNA, or cDNA reads that were oriented \
                 beforehand); otherwise use `--strand-filter both`.",
                pct(f),
                strand_name(ctx.strand_filter)
            ),
        });
    }

    // Reads that the other per-alignment filters leave without an alignment.
    let other = dt.reads_no_valid_aln().saturating_sub(wrong_strand);
    let f = share(other, mapped);
    if f > MAX_ALIGNMENT_FILTER_DROP {
        warnings.push(RunWarning {
            code: "alignment_filters",
            fraction: f,
            threshold: MAX_ALIGNMENT_FILTER_DROP,
            message: format!(
                "{other} of {mapped} mapped reads ({}) have no alignment that passes the alignment \
                 filters (aligned length, distance from the transcript ends, a positive score). \
                 Check the filter settings recorded under `filter_options`.",
                pct(f)
            ),
        });
    }

    // Reads without any mapping.
    let unmapped = dt.reads_unmapped();
    let f = share(unmapped, seen);
    match ctx.map_strand {
        Some(strand) if one_strand(strand) && f > MAX_UNMAPPED_ONE_STRAND => {
            warnings.push(RunWarning {
                code: "map_strand_unmapped",
                fraction: f,
                threshold: MAX_UNMAPPED_ONE_STRAND,
                message: format!(
                    "{unmapped} of {seen} reads ({}) did not map with `--map-strand {}`. If the \
                     reads are not all oriented to the transcript strand, about half of them \
                     cannot map this way: drop `--map-strand` (or orient the reads first).",
                    pct(f),
                    strand_name(strand)
                ),
            });
        }
        _ if f > MAX_UNMAPPED => {
            warnings.push(RunWarning {
                code: "low_mapping_rate",
                fraction: f,
                threshold: MAX_UNMAPPED,
                message: format!(
                    "{unmapped} of {seen} reads ({}) did not map. Check that the reference matches \
                     the sample, and that `--seq-tech` matches the reads.",
                    pct(f)
                ),
            });
        }
        _ => {}
    }

    warnings
}

#[cfg(test)]
mod tests {
    use super::*;

    const CTX: RunContext = RunContext {
        strand_filter: Strand::Unknown,
        map_strand: None,
        min_aligned_fraction: 0.5,
    };

    fn codes(w: &[RunWarning]) -> Vec<&'static str> {
        w.iter().map(|w| w.code).collect()
    }

    #[test]
    fn index_without_gene_groups() {
        use crate::prog_opts::IndexGenes;
        let genes = |source: &str, genes, has| IndexGenes {
            source: source.to_string(),
            genes: Some(genes),
            transcripts: Some(10),
            has_gene_groups: has,
        };
        let w = check_index_genes(&genes("none", 10, Some(false))).unwrap();
        assert_eq!(w.code, "index_without_gene_groups");
        assert_eq!((w.fraction, w.threshold), (1.0, 1.0));
        assert!(w.message.contains("no gene IDs"));
        assert!(w.message.contains("--t2g") && w.message.contains("--infer-genes"));
        // Inference that grouped nothing: a table is the fix.
        let w = check_index_genes(&genes("inferred", 10, Some(false))).unwrap();
        assert!(w.message.contains("--t2g") && !w.message.contains("`--infer-genes`"));
        assert!(check_index_genes(&genes("fasta", 4, Some(true))).is_none());
        // RAD input: not known, no warning.
        assert!(check_index_genes(&IndexGenes::unknown(Some(10))).is_none());
        let v = serde_json::to_value(IndexGenes::unknown(Some(10))).unwrap();
        assert_eq!(v["source"], "unknown");
        assert!(v["genes"].is_null() && v["has_gene_groups"].is_null());
        assert_eq!(v["transcripts"], 10);
    }

    #[test]
    fn a_clean_run_has_no_warnings() {
        // 5% unmapped, 5% of mapped reads under the aligned fraction, 1% without a valid one
        let dt = DiscardTable::from_read_counts(50, 10, 0, 48, 892);
        assert!(check(&dt, &CTX).is_empty());
        assert!(check(&DiscardTable::new(), &CTX).is_empty());
    }

    #[test]
    fn aligned_fraction_drop_warns_above_the_threshold() {
        let dt = DiscardTable::from_read_counts(0, 0, 0, 200, 800);
        let w = check(&dt, &CTX);
        assert_eq!(codes(&w), ["aligned_fraction_filter"]);
        assert!((w[0].fraction - 0.2).abs() < 1e-12);
        assert!(w[0].message.contains("200 of 1000 mapped reads (20.0%)"));
        // the advice agrees with the guide it links to
        assert!(w[0].message.contains(
            "sequencing-technologies/#when-the-aligned-fraction-filter-drops-many-reads"
        ));
        assert!(w[0].message.contains("keeping the poly(A) tails"));
        assert!(w[0].message.contains("For gene-level analyses"));
        // exactly at the threshold: no warning
        let dt = DiscardTable::from_read_counts(0, 0, 0, 100, 900);
        assert!(check(&dt, &CTX).is_empty());
    }

    #[test]
    fn strand_filter_warns_only_when_it_is_on() {
        // 300 of 1000 mapped reads have alignments on the removed strand only
        let dt = DiscardTable::from_read_counts(0, 300, 300, 0, 700);
        assert!(check(&dt, &CTX).is_empty());
        let ctx = RunContext {
            strand_filter: Strand::Forward,
            ..CTX
        };
        let w = check(&dt, &ctx);
        assert_eq!(codes(&w), ["strand_filter"]);
        assert!(w[0].message.contains("keeps `fw`"));
        // the other alignment filters are reported apart from the strand filter
        let dt = DiscardTable::from_read_counts(0, 300, 100, 0, 700);
        assert_eq!(codes(&check(&dt, &CTX)), ["alignment_filters"]);
    }

    #[test]
    fn unmapped_reads_warn_sooner_under_map_strand() {
        // 40% unmapped
        let dt = DiscardTable::from_read_counts(400, 0, 0, 0, 600);
        assert!(check(&dt, &CTX).is_empty());
        let fw = RunContext {
            map_strand: Some(Strand::Forward),
            ..CTX
        };
        assert_eq!(codes(&check(&dt, &fw)), ["map_strand_unmapped"]);
        // `--map-strand both` is no restriction
        let both = RunContext {
            map_strand: Some(Strand::Unknown),
            ..CTX
        };
        assert!(check(&dt, &both).is_empty());
        // 60% unmapped
        let dt = DiscardTable::from_read_counts(600, 0, 0, 0, 400);
        assert_eq!(codes(&check(&dt, &CTX)), ["low_mapping_rate"]);
        assert_eq!(codes(&check(&dt, &both)), ["low_mapping_rate"]);
        assert_eq!(codes(&check(&dt, &fw)), ["map_strand_unmapped"]);
    }

    #[test]
    fn em_warns_only_when_it_stops_at_the_limit() {
        let mut em = EmConvergence {
            niter: 412,
            max_iter: 1000,
            converged: true,
            rel_diff: 9e-4,
            thresh: 1e-3,
        };
        assert_eq!(check_em(&em, EmKind::Standard), None);
        assert_eq!(check_em(&em, EmKind::Truncation), None);
        let v = em_info(&em, EmKind::Truncation);
        assert_eq!(v["kind"], "truncation");
        assert_eq!(
            (v["iterations"].as_u64(), v["converged"].as_bool()),
            (Some(412), Some(true))
        );
        em.niter = 1000;
        em.converged = false;
        em.rel_diff = 0.0123;
        let w = check_em(&em, EmKind::Standard).expect("a warning");
        assert_eq!(w.code, "em_not_converged");
        assert!(w.message.contains("default 1,000"));
        assert!(w.message.contains("usually stable"));
        let t = check_em(&em, EmKind::Truncation).expect("a warning");
        assert_eq!(t.code, "truncation_em_not_converged");
        assert_eq!((t.fraction, t.threshold), (0.0123, 1e-3));
        assert!(t.message.contains("default 5,000 with --model-truncation"));
        assert_eq!(em_info(&em, EmKind::Standard)["kind"], "standard");
        assert_eq!((w.fraction, w.threshold), (0.0123, 1e-3));
        assert!(w.message.contains("iteration limit (1000 iterations"));
        assert!(t.message.contains("truncation EM (--model-truncation)"));
        assert!(w.message.contains("1.230e-2"));
        assert!(w.message.contains("--max-em-iter"));
    }

    #[test]
    fn warnings_serialize_with_their_code_and_numbers() {
        let dt = DiscardTable::from_read_counts(0, 0, 0, 200, 800);
        let v = serde_json::to_value(check(&dt, &CTX)).unwrap();
        assert_eq!(v[0]["code"], "aligned_fraction_filter");
        assert_eq!(v[0]["threshold"], 0.1);
        assert!(v[0]["message"].is_string());
    }
}
