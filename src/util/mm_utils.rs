use crate::util::mm_utils;
use minimap2_sys::{MmIdx, mm_idx_seq_t};
use std::ffi::CStr;
use std::sync::Arc;

//https://github.com/lh3/minimap2/blob/618d33515e5853c4576d5a3d126fdcda28f0e8a4/mmpriv.h#L32
// #define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)
const fn mm_seq4_get(seq: *const u32, i: isize) -> char {
    // SAFETY: i is obtained from a valid `mm_idx_seq_t` object
    let r = (unsafe { seq.offset(i >> 3).read_unaligned() } >> ((i) & 7) << 2) & 0xf;
    match r {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

pub(crate) fn get_sequence(idx: &Arc<MmIdx>, st: &mm_idx_seq_t) -> String {
    let offset = st.offset;
    let len = st.len;
    let end = offset + len as u64;
    let str_base: *const u32 = (idx).S;
    let mut st = String::with_capacity(len as usize);
    for i in offset..end {
        st.push(mm_seq4_get(str_base, i as isize));
    }
    st
}

pub struct MMIdxNameSeqIter {
    idx: Arc<MmIdx>,
    curr_seq: isize,
    nseq: isize,
}

impl MMIdxNameSeqIter {
    pub fn from_idx(idx_in: &Arc<MmIdx>) -> Self {
        Self {
            idx: idx_in.clone(),
            curr_seq: 0,
            nseq: idx_in.n_seq as isize,
        }
    }
}

impl Iterator for MMIdxNameSeqIter {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_seq < self.nseq {
            // SAFETY: curr_seq < nseq, so this access must be in bounds
            let seq = unsafe { *(self.idx).seq.offset(self.curr_seq) };
            let sequence = mm_utils::get_sequence(&self.idx, &seq);
            let c_str = unsafe { CStr::from_ptr(seq.name) };
            let rust_str = c_str.to_str().unwrap().to_string();
            self.curr_seq += 1;
            Some((rust_str, sequence))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let rem = (self.curr_seq - self.nseq) as usize;
        (rem, Some(rem))
    }
}

impl ExactSizeIterator for MMIdxNameSeqIter {}
