use crate::variables::TranscriptInfo;
use nested_intervals::IntervalSet;
use rayon::prelude::*;
use std::sync::{RwLock, Arc};

pub fn uniform_prob(txps: &mut Vec<TranscriptInfo>, threads: usize) -> Vec<usize> {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));

    txps.par_iter_mut().enumerate().for_each(|(i, t)| {

        //fill out the number of reads and discarded reads
        {
            let mut w1 =num_reads.write().unwrap();
            (*w1)[i] = t.ranges.len();
        }

        //computing the coverage probability
        let len = t.len.get() as u32; //transcript length
        let interval_set = IntervalSet::new(&t.ranges).expect("couldn't build interval set");
        let mut interval_set = interval_set.merge_connected();
        let covered = interval_set.covered_units();
        let prob_uniform = (covered as f64) / (len as f64);

        t.coverage = prob_uniform;
    });

    let mut nr1: Vec<usize> = vec![];
    {
        nr1 = num_reads.read().unwrap().to_vec().clone();
    }
    
    nr1
}