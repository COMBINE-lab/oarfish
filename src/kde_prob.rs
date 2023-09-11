use crate::variables::TranscriptInfo;
use kernel_density_estimation::prelude::*;
use nested_intervals::IntervalSet;
extern crate libc;
use rayon::prelude::*;
use std::sync::{RwLock, Arc};

#[link(name = "mlpack_kde_wrapper")] // Link against the shared library
extern "C" {
    fn mlpack_kde_wrapper(
        referenceData: *const libc::c_double, referenceRows: libc::size_t, referenceCols: libc::size_t,
        queryData: *const libc::c_double, queryRows: libc::size_t, queryCols: libc::size_t,
        bandwidth_val: libc::c_double, monte_carlo_tag: bool,
        prediction: *mut libc::c_double, predictionSize: libc::size_t,
    ) -> u32;
}


//kernel density estimation used to obtain probability of each alignment
pub fn kde_prob(txps: &mut Vec<TranscriptInfo>, rate: &str, threads: usize) -> Vec<usize> {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));  

    match rate{
        "1D" => {
            txps.par_iter_mut().enumerate().for_each(|(i, t)| {

                {
                    let mut w =num_reads.write().unwrap();
                    (*w)[i] = t.ranges.len();
                }
                
                let mut observations: Vec<_> = t.ranges
                    .iter()
                    .flat_map(|range| range.clone())
                    .collect();

                observations.sort();
                let f32_observations: Vec<f32> = observations.iter().map(|&x| x as f32).collect();
                let bandwidth = Scott;
                let kernel = Normal;
                let kde = KernelDensityEstimator::new(f32_observations, bandwidth, kernel);

                let tlen = t.len.get(); //transcript length

                let group_size: usize = 10;
                let (start_elemtns, middle_elements) = sample_middle_elements(tlen, group_size);
                let cdf_t = kde.cdf(middle_elements.as_slice());

                t.coverage_prob = vec![std::f32::NAN; tlen + 1];
                start_elemtns.iter()
                    .zip(cdf_t.iter())
                    .for_each(|(data, cdf)| {
                        t.coverage_prob[*data as usize] = *cdf;
                    });
            });
        }
        "2D" => {
                    
        }
        _ => {
            panic!("{:?} rate is not defined in the program", rate);
        }
    }

    let mut nr: Vec<usize> = vec![];
    {
        nr = num_reads.read().unwrap().to_vec().clone();
    }
    
    nr

}


pub fn kde_c_prob(txps: &mut Vec<TranscriptInfo>, rate: &str, threads: usize) -> Vec<usize> {

    rayon::ThreadPoolBuilder::new()
    .num_threads(threads)
    .build()
    .unwrap();

    let num_reads = Arc::new(RwLock::new(vec![0; txps.len()]));  

    match rate{
        "1D" => {

            txps.par_iter_mut().enumerate().for_each(|(i, t)| {
                
                {
                    let mut w =num_reads.write().unwrap();
                    (*w)[i] = t.ranges.len();
                }

                let mut observations: Vec<_> = t.ranges
                    .iter()
                    .flat_map(|range| range.clone())
                    .collect();
                observations.sort();
                
                let mut f32_observations: Vec<f64> = observations.iter().map(|&x| x as f64).collect();
                let observation_rows= 1;
                let observation_cols = f32_observations.len();

                let bandwidth = 1.0;
                let monte_carlo_tag: bool = true;

                let tlen = t.len.get(); //transcript length
                let group_size: usize = 10;
                let (start_elemtns, middle_elements) = sample_middle_elements(tlen, group_size);
                let mut middle_elements_f64: Vec<f64> = middle_elements.iter().map(|&element| element as f64).collect();
                let dataset_rows = 1;
                let dataset_cols = middle_elements_f64.len();

                let mut prediction = vec![0.0; middle_elements_f64.len()];
                let prediction_size = middle_elements_f64.len();
                
                if !f32_observations.is_empty(){
                    let result = unsafe {
                        mlpack_kde_wrapper(
                            f32_observations.as_mut_ptr(), observation_rows, observation_cols,
                            middle_elements_f64.as_mut_ptr(), dataset_rows, dataset_cols,
                            bandwidth, monte_carlo_tag, prediction.as_mut_ptr(), prediction_size
                        )
                    };
                    match result {
                        1 => {
                            panic!("CDFFailed in mlpack wrapper");
                        }
                        2 => {
                            panic!("PointerFailed in mlpack wrapper");
                        }
                        3 => {
                            panic!("PointerSizeFailed in mlpack wrapper");
                        }
                        4 => {
                            panic!("InvalidInput in mlpack wrapper");
                        }
                        _ => {

                        }
                    }
                }
                t.coverage_prob = vec![std::f32::NAN; tlen + 1];
                start_elemtns.iter()
                    .zip(prediction.iter())
                    .for_each(|(data, cdf)| {
                        t.coverage_prob[*data as usize] = *cdf as f32;
                    });
            });

        }
        "2D" => {

            let mut tlen: Vec<f64> = vec![];
            let mut tcovered: Vec<f64> = vec![];
            for t in txps.iter_mut() {
                
                let len = t.len.get() as u32; //transcript length
                tlen.push(len as f64);

                let interval_set = IntervalSet::new(&t.ranges).expect("couldn't build interval set");
                let mut interval_set = interval_set.merge_connected();
                let covered = interval_set.covered_units();
                let fraction_covered = (covered as f64) / (len as f64);
                tcovered.push(fraction_covered);
            }
            let observation_rows= 2;
            let observation_cols = tlen.len();
            let observations: Vec<f64> = tlen.into_iter().zip(tcovered.into_iter()).flat_map(|(num1, num2)| vec![num1, num2]).collect::<Vec<f64>>();

            let bandwidth = 1.0;
            let monte_carlo_tag: bool = true;

            //https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html
            //https://docs.rs/rayon/latest/rayon/iter/trait.IntoParallelRefMutIterator.html

            txps.par_iter_mut().enumerate().for_each(|(i, t)| {
                
                {
                    let mut w =num_reads.write().unwrap();
                    (*w)[i] = t.ranges.len();
                }

                let len = t.len.get() as f64; //transcript length
                let rlen: Vec<f64> = t.ranges.iter().map(|range| (range.end - range.start) as f64).collect();
                let rcovered: Vec<f64> = t.ranges.iter().map(|range| (range.end - range.start) as f64 / len).collect();
                let dataset_rows = 2;
                let dataset_cols = rlen.len();
                let dataset: Vec<f64> = rlen.into_iter().zip(rcovered.into_iter()).flat_map(|(num1, num2)| vec![num1, num2]).collect::<Vec<f64>>();

                let prediction_size = if dataset_cols != 0 {dataset_cols} else {1};
                let mut prediction = vec![0.0; prediction_size];
                

                if !dataset.is_empty(){
                    unsafe {
                        mlpack_kde_wrapper(
                            observations.as_ptr(), observation_rows, observation_cols,
                            dataset.as_ptr(), dataset_rows, dataset_cols,
                            bandwidth, monte_carlo_tag, prediction.as_mut_ptr(), prediction_size
                        );
                    }
                }
                t.coverage_prob = prediction.iter().map(|&x| x as f32).collect();
            });


        }
        _ => {
            panic!("{:?} rate is not defined in the program", rate);
        }
    }
        
    let mut nr: Vec<usize> = vec![];
    {
        nr = num_reads.read().unwrap().to_vec().clone();
    }
    
    nr

}

//sampling function for KDE
pub fn sample_middle_elements(tlen: usize, group_size: usize) -> (Vec<f32>, Vec<f32>) {

    let num_groups = (tlen as f32 / group_size as f32).ceil();
    let mut start_elemtns: Vec<f32> = Vec::with_capacity(num_groups as usize);
    let mut middle_elements: Vec<f32> = Vec::with_capacity(num_groups as usize);

    for i in 0..num_groups as usize {
        let start_index = i * group_size;
        let middle_index = start_index as f32 + (group_size as f32 / 2.0);
        start_elemtns.push(start_index as f32);
        middle_elements.push(middle_index);
    }

    (start_elemtns, middle_elements)
}
