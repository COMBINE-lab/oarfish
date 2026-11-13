//! Warnings about a run as a whole: conditions that do not stop quantification, but that the
//! user should look at before trusting the result. Each is written to the log and recorded
//! in the `warnings` array of `meta_info.json`.
//!
//! All of them are computed from counts the run keeps anyway (the discard table), so they
//! cost nothing per read.

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
    /// The share of reads that triggered the warning.
    pub fraction: f64,
    /// The share above which the warning is issued.
    pub threshold: f64,
    /// What happened and what to check, for a person.
    pub message: String,
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
                 adapters or primers are the usual cause: trim them (see \
                 https://combine-lab.github.io/oarfish/guides/sequencing-technologies/), or \
                 keep these reads with `--min-aligned-fraction 0` or \
                 `--filter-group no-filters`.{one_strand_note}",
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
    fn warnings_serialize_with_their_code_and_numbers() {
        let dt = DiscardTable::from_read_counts(0, 0, 0, 200, 800);
        let v = serde_json::to_value(check(&dt, &CTX)).unwrap();
        assert_eq!(v[0]["code"], "aligned_fraction_filter");
        assert_eq!(v[0]["threshold"], 0.1);
        assert!(v[0]["message"].is_string());
    }
}
