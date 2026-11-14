pub(crate) const MIN_READ_THRESH: f64 = 1e-5;
/// The truncation EM's convergence test looks only at transcripts with more than this many
/// reads (salmon's cutoff for the same test). Below it, abundances decaying toward zero keep a
/// relative change above the threshold for thousands of iterations while moving a negligible
/// number of reads.
pub(crate) const TRUNC_CONV_MIN_READS: f64 = 0.01;
pub(crate) const EM_DENOM_THRESH: f64 = 1e-30_f64;
pub(crate) const EMPTY_READ_NAME: &str = "no_read_name_available";
