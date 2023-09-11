use crate::variables::TranscriptInfo;
use rayon::prelude::*;
use std::sync::{RwLock, Arc};

pub fn no_coverage_prob(txps: &mut Vec<TranscriptInfo>, threads: usize) -> Vec<usize> {

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
        t.coverage = 1.0 as f64;
    });

    let mut nr1: Vec<usize> = vec![];
    {
        nr1 = num_reads.read().unwrap().to_vec().clone();
    }
    
    nr1

}