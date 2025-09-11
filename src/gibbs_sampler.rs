use crate::util::oarfish_types::{AlnInfo, EMInfo, TranscriptInfo};
use rand::Rng;
use rand::distr::{Distribution, Uniform};

pub fn collapsed_gibbs<
    'a,
    I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a,
    F: Fn() -> I,
>(
    em_info: &'a EMInfo,
    num_steps: usize,
    make_iter: F,
    do_log: bool,
) -> Vec<f64> {
    let eq_map = &em_info.eq_map;
    let fops = &eq_map.filter_opts;
    let tinfo: &[TranscriptInfo] = em_info.txp_info;
    let max_iter = num_steps;
    let total_weight: f64 = eq_map.num_aligned_reads() as f64;

    // initialize the estimated counts for the EM procedure
    let mut prev_est_counts: Vec<f64>;
    let mut curr_counts: Vec<u64> = vec![0u64; tinfo.len()];

    // uniform, length normalized abundance
    let avg = total_weight / (tinfo.len() as f64);
    prev_est_counts = vec![avg; tinfo.len()];

    vec![]
}
