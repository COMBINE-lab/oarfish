use rand::distributions::{Distribution, Uniform};
use rand::Rng;

/// Wraps an iterator `I` to produce a new iterator that samples the original
/// with replacement.
pub struct SampleWithReplaceIter<'b, I: ExactSizeIterator> {
    pub inds: &'b [usize],
    pub iter: I,
    pub idx: usize,
    pub e_ptr: usize,
    pub prev: Option<I::Item>,
}

impl<'b, I: ExactSizeIterator> Iterator for SampleWithReplaceIter<'b, I>
where
    I::Item: Copy,
{
    type Item = I::Item;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next_ind) = self.inds.get(self.idx) {
            // always makes sure we have a previous
            // element. This will get the first
            // element but wont yet increment
            // the pointer (e_ptr)
            if self.prev.is_none() {
                self.prev = self.iter.next();
            }
            // next element is not a duplicate
            if self.e_ptr != *next_ind {
                // otherwise, call next until we have the
                // correct element
                loop {
                    self.prev = self.iter.next();
                    self.e_ptr += 1;
                    if self.e_ptr == *next_ind {
                        break;
                    }
                }
            }
            self.idx += 1;
            self.prev
        } else {
            None
        }
    }
}

pub fn get_sample_inds<R: Rng + ?Sized>(n: usize, rng: &mut R) -> Vec<usize> {
    let dist = Uniform::new(0, n);
    let mut inds: Vec<usize> = dist.sample_iter(rng).take(n).collect();
    inds.sort_unstable();
    inds
}
