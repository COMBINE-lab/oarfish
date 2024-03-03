use rand::distributions::{Distribution, Uniform};
use rand::Rng;

pub fn get_sample_inds<R: Rng + ?Sized>(n: usize, rng: &mut R) -> Vec<usize> {
    let dist = Uniform::new(0, n);
    let mut inds: Vec<usize> = dist.sample_iter(rng).take(n).collect();
    inds.sort_unstable();
    inds
}
