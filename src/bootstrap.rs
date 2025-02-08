use rand::distr::{Distribution, Uniform};
use rand::Rng;

/// Get a random uniform sample of `n` numbers in the range [0,n).
/// Duplicates are explicitly allowed. The numbers are returned in
/// sorted order.
pub fn get_sample_inds<R: Rng + ?Sized>(n: usize, rng: &mut R) -> Vec<usize> {
    let dist = Uniform::new(0, n);
    let mut inds: Vec<usize> = dist
        .expect("could not create distribution")
        .sample_iter(rng)
        .take(n)
        .collect();
    inds.sort_unstable();
    inds
}
