use crate::TranscriptInfo;
use crate::util::oarfish_types::InMemoryAlignmentStore;
use anyhow;
use itertools::izip;

pub struct CountInfo {
    pub unique_count: u32,
    pub total_count: u32,
    #[allow(dead_code)]
    pub expected_count: f64,
}

impl CountInfo {
    pub fn new() -> Self {
        Self {
            unique_count: 0,
            total_count: 0,
            expected_count: 0.0,
        }
    }
}

pub fn get_aux_counts(
    store: &InMemoryAlignmentStore<'_>,
    txps: &[TranscriptInfo],
) -> anyhow::Result<Vec<CountInfo>> {
    let mut cinfo = Vec::with_capacity(txps.len());

    for _ in txps {
        cinfo.push(CountInfo::new())
    }

    for (alns, probs, coverage_probs) in store.iter() {
        let is_unique = alns.len() == 1;
        for (a, _p, _cp) in izip!(alns, probs, coverage_probs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            if let Some(ref mut ci) = cinfo.get_mut(target_id) {
                ci.total_count += 1;
                if is_unique {
                    ci.unique_count += 1;
                }
            }
        }
    }

    Ok(cinfo)
}
