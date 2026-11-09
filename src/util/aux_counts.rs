use crate::TranscriptInfo;
use crate::util::oarfish_types::InMemoryAlignmentStore;
use anyhow;

pub struct CountInfo {
    pub unique_count: u32,
    pub total_count: u32,
    #[allow(dead_code)]
    pub expected_count: f64,
}

impl CountInfo {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            unique_count: 0,
            total_count: 0,
            expected_count: 0.0,
        }
    }
}

/// Alignments per chunk below which [`get_aux_counts`] adds no further chunk.
const MIN_CHUNK_ALIGNMENTS: usize = 1 << 18;

/// Total and unique alignment counts per transcript. The reads are counted in chunks in
/// parallel on the current rayon pool (one per thread, fewer for a small store), and each
/// transcript's per-chunk counts are summed afterwards; the counts are integers, so the result
/// does not depend on the chunking.
pub fn get_aux_counts(
    store: &InMemoryAlignmentStore<'_>,
    txps: &[TranscriptInfo],
) -> anyhow::Result<Vec<CountInfo>> {
    use rayon::prelude::*;

    let m = txps.len();
    let n = store.len();
    let bounds = store.boundaries();
    let alns = &store.alignments;
    let nchunks = rayon::current_num_threads()
        .min(alns.len() / MIN_CHUNK_ALIGNMENTS)
        .min(n)
        .max(1);

    // Counts `[total, unique]` per transcript over reads `lo..hi`.
    let count = |lo: usize, hi: usize| -> Vec<[u32; 2]> {
        let mut c = vec![[0u32; 2]; m];
        for r in lo..hi {
            let (start, end) = (bounds[r], bounds[r + 1]);
            let is_unique = (end - start == 1) as u32;
            for a in &alns[start..end] {
                if let Some(ci) = c.get_mut(a.ref_id as usize) {
                    ci[0] += 1;
                    ci[1] += is_unique;
                }
            }
        }
        c
    };

    let to_info = |[total_count, unique_count]: [u32; 2]| CountInfo {
        unique_count,
        total_count,
        expected_count: 0.0,
    };

    if nchunks == 1 {
        return Ok(count(0, n).into_iter().map(to_info).collect());
    }

    let parts: Vec<Vec<[u32; 2]>> = (0..nchunks)
        .into_par_iter()
        .map(|i| count(i * n / nchunks, (i + 1) * n / nchunks))
        .collect();

    let mut cinfo = Vec::with_capacity(m);
    (0..m)
        .into_par_iter()
        .map(|t| {
            let mut s = [0u32; 2];
            for p in &parts {
                s[0] += p[t][0];
                s[1] += p[t][1];
            }
            to_info(s)
        })
        .collect_into_vec(&mut cinfo);
    Ok(cinfo)
}
