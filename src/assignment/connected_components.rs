use crate::util::oarfish_types::{AlnInfo, EMInfo};
use aph_disjoint_set::DisjointSet;
use std::collections::HashMap;

#[derive(Default, Clone, Debug, PartialEq, PartialOrd)]
pub struct ClusterSize {
    pub ntxp: u32,
    pub nread: u32,
}

#[allow(dead_code)]
impl ClusterSize {
    pub fn new() -> Self {
        ClusterSize { ntxp: 0, nread: 0 }
    }

    pub fn with_txp_count(ntxp: u32) -> Self {
        ClusterSize { ntxp, nread: 0 }
    }

    pub fn with_read_count(nread: u32) -> Self {
        ClusterSize { ntxp: 0, nread }
    }

    pub fn with_txp_read_count(ntxp: u32, nread: u32) -> Self {
        ClusterSize { ntxp, nread }
    }

    pub fn add_read(&mut self) {
        self.nread += 1;
    }

    pub fn add_txp(&mut self) {
        self.ntxp += 1;
    }
}

pub struct TranscriptConnectedComponentLabeling {
    pub labels: Vec<u32>,
    pub cc_sizes: HashMap<u32, ClusterSize>,
}

impl TranscriptConnectedComponentLabeling {
    pub fn from_labels_and_sizes(labels: Vec<u32>, cc_sizes: HashMap<u32, ClusterSize>) -> Self {
        Self { labels, cc_sizes }
    }
}

pub fn get_connected_components<
    'a,
    I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a,
    F: Fn() -> I,
>(
    em_info: &'a EMInfo,
    make_iter: F,
) -> TranscriptConnectedComponentLabeling {
    let ntrans = em_info.txp_info.len();
    let mut txp_uf = DisjointSet::new(ntrans);

    for (alns, _probs, _coverage_probs) in make_iter() {
        if alns.len() > 1 {
            let first_txp = alns[0].ref_id as usize;
            for a in alns.iter().skip(1) {
                let target_id = a.ref_id as usize;
                txp_uf.union(target_id, first_txp);
            }
        }
    }

    let mut labels = Vec::with_capacity(ntrans);
    let mut cluster_sizes = HashMap::<u32, ClusterSize>::with_capacity(10_000);
    for tid in 0..ntrans {
        let rep = txp_uf.get_root(tid).into_inner() as u32;
        labels.push(rep);
        cluster_sizes
            .entry(rep)
            .and_modify(|cluster_size| cluster_size.add_txp())
            .or_insert(ClusterSize::with_txp_count(1));
    }

    // now add the read count to each cluster
    for (alns, _probs, _coverage_probs) in make_iter() {
        if let Some(aln) = alns.first() {
            let tid = aln.ref_id as usize;
            // get the cluster id
            let rep = txp_uf.get_root(tid).into_inner() as u32;
            cluster_sizes
                .entry(rep)
                .and_modify(|cluster_size| cluster_size.add_read());
        }
    }

    TranscriptConnectedComponentLabeling::from_labels_and_sizes(labels, cluster_sizes)
}
