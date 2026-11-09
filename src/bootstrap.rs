use rand::Rng;
use rand::distr::{Distribution, Uniform};

/// Get a random uniform sample of `n` numbers in the range [0,n).
/// Duplicates are explicitly allowed. The numbers are returned in
/// sorted order. With `n == 0` the sample is empty.
pub fn get_sample_inds<R: Rng + ?Sized>(n: usize, rng: &mut R) -> Vec<usize> {
    if n == 0 {
        return Vec::new();
    }
    let dist = Uniform::new(0, n);
    let mut inds: Vec<usize> = dist
        .expect("could not create distribution")
        .sample_iter(rng)
        .take(n)
        .collect();
    inds.sort_unstable();
    inds
}

#[cfg(test)]
mod tests {
    use super::get_sample_inds;

    #[test]
    fn sample_of_nothing_is_empty() {
        assert!(get_sample_inds(0, &mut rand::rng()).is_empty());
        let s = get_sample_inds(5, &mut rand::rng());
        assert_eq!(s.len(), 5);
        assert!(s.windows(2).all(|w| w[0] <= w[1]) && s.iter().all(|&i| i < 5));
    }
}
