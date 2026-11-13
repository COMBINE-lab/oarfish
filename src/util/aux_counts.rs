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
/// Alignments per transcript a chunk must count before [`get_aux_counts`] gives it its own
/// per-transcript counts (8 bytes per transcript): this keeps those buffers to at most 4 bytes
/// per alignment, so a large transcriptome with few alignments is counted in fewer chunks
/// instead of holding many mostly-zero arrays.
const MIN_ALIGNMENTS_PER_TRANSCRIPT: usize = 2;

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
        .min(alns.len() / (MIN_ALIGNMENTS_PER_TRANSCRIPT * m).max(1))
        .min(n)
        .max(1);

    // Counts `[total, unique]` per transcript over reads `lo..hi`.
    let count = |lo: usize, hi: usize| -> Vec<[u32; 2]> {
        let mut c = vec![[0u32; 2]; m];
        // the transcripts of one read; a read with several alignments to one transcript
        // counts once for it, and is unique if that is its only transcript
        let mut ids = Vec::<u32>::new();
        for r in lo..hi {
            let (start, end) = (bounds[r], bounds[r + 1]);
            ids.clear();
            ids.extend(alns[start..end].iter().map(|a| a.ref_id));
            if ids.len() > 1 {
                ids.sort_unstable();
                ids.dedup();
            }
            let is_unique = (ids.len() == 1) as u32;
            for &id in &ids {
                if let Some(ci) = c.get_mut(id as usize) {
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

#[cfg(test)]
mod tests {
    use super::get_aux_counts;
    use crate::em_tests::{store_of, tinfo};
    use noodles_sam::Header;

    #[test]
    fn reads_are_counted_once_per_transcript() {
        let h = Header::default();
        let mut txps = tinfo(&[1000, 1000, 1000]);
        let reads = vec![
            // two alignments to transcript 0 and none elsewhere: one read, and a unique one
            vec![(0, 1.0), (0, 0.5)],
            // transcript 1 twice and transcript 2: one ambiguous read for each
            vec![(1, 1.0), (2, 1.0), (1, 0.5)],
            vec![(2, 1.0)],
        ];
        let store = store_of(&h, &reads, &mut txps);
        let c = get_aux_counts(&store, &txps).unwrap();
        let got: Vec<(u32, u32)> = c.iter().map(|c| (c.unique_count, c.total_count)).collect();
        assert_eq!(got, [(1, 1), (0, 1), (1, 2)]);
    }
}
