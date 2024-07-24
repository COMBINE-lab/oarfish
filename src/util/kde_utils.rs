use crate::util::oarfish_types::{InMemoryAlignmentStore, TranscriptInfo};
use itertools::izip;
use kders::kde::{GridDimensions, KDEModel};
use tracing::info;

pub fn get_kde_model(
    txps: &[TranscriptInfo],
    store: &InMemoryAlignmentStore,
) -> anyhow::Result<KDEModel> {
    let mut max_x: f64 = 0_f64;
    let mut max_y: f64 = 0_f64;

    for (ainfs, _aprobs, _cprobs) in store.iter() {
        for ainf in ainfs {
            let txp_len = txps[ainf.ref_id as usize].lenf;
            let aln_len = ainf.alignment_span() as f64;
            max_x = max_x.max(txp_len);
            max_y = max_y.max(aln_len);
        }
    }

    let gd = GridDimensions {
        width: max_x as usize + 1,
        height: max_y as usize + 1,
    };

    info!("KDE grid maxima = ({}, {})", gd.width, gd.height);

    let kernel_bandwidth = 50_f64;
    let bin_width = 25_usize;

    let mut grid = kders::kde::KDEGrid::new(gd, bin_width, Some(kernel_bandwidth));

    for (ainfs, _aprobs, _cprobs) in store.iter() {
        let w = 1. / (ainfs.len() as f64);
        for ainf in ainfs {
            let txp_len = txps[ainf.ref_id as usize].lenf;
            let aln_len = ainf.alignment_span();
            grid.add_observation(txp_len as usize, aln_len as usize, w);
        }
    }

    let density = grid.get_kde()?;
    Ok(density)
}

#[allow(unused)]
pub fn refresh_kde_model(
    txps: &[TranscriptInfo],
    store: &InMemoryAlignmentStore,
    kde_model: &KDEModel,
    counts: &[f64],
) -> anyhow::Result<KDEModel> {
    let gd = GridDimensions {
        width: kde_model.width,
        height: kde_model.height,
    };

    info!("KDE grid maxima = ({}, {})", gd.width, gd.height);

    let kernel_bandwidth = 50_f64;
    let bin_width = 25_usize;

    let mut grid = kders::kde::KDEGrid::new(gd, bin_width, Some(kernel_bandwidth));

    for (ainfs, aprobs, cprobs) in store.iter() {
        let mut denom = 0.0_f64;
        for (a, p, _cp) in izip!(ainfs, aprobs, cprobs) {
            // Compute the probability of assignment of the
            // current read based on this alignment and the
            // target's estimated abundance.
            let target_id = a.ref_id as usize;
            let prob = *p as f64;
            let cov_prob = 1.0; //if model_coverage { *cp } else { 1.0 };
            let txp_len = txps[target_id].lenf;
            let aln_len = a.alignment_span();
            let flprob = kde_model[(txp_len as usize, aln_len as usize)];
            denom += counts[target_id] * prob * cov_prob * flprob;
        }

        // If this read can be assigned
        if denom > crate::util::constants::EM_DENOM_THRESH {
            // Loop over all possible assignment locations and proportionally
            // allocate the read according to our model and current parameter
            // estimates.
            for (a, p, _cp) in izip!(ainfs, aprobs, cprobs) {
                let target_id = a.ref_id as usize;
                let prob = *p as f64;
                let cov_prob = 1.0; //if model_coverage { *cp } else { 1.0 };
                let txp_len = txps[target_id].lenf;
                let aln_len = a.alignment_span();
                let flprob = kde_model[(txp_len as usize, aln_len as usize)];
                let w = (counts[target_id] * prob * cov_prob * flprob) / denom;
                grid.add_observation(txp_len as usize, aln_len as usize, w);
            }
        }
    }
    info!("filled grid; computing KDE");
    grid.get_kde()
}
