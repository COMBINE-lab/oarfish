pub(crate) const MIN_READ_THRESH: f64 = 1e-5;
// Transcripts below this estimated read count (in either of two successive
// iterates) are excluded from the convergence criterion: their relative change
// never decays for geometrically dying parameters, while their absolute effect
// on any reported quantity is below output precision.
pub(crate) const MIN_ACTIVE_COUNT: f64 = 1e-2;
pub(crate) const EM_DENOM_THRESH: f64 = 1e-30_f64;
pub(crate) const EMPTY_READ_NAME: &str = "no_read_name_available";
