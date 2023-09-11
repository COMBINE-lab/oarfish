extern crate libc;

#[link(name = "mlpack_kde_wrapper")] // Link against the shared library
extern "C" {
    fn mlpack_kde_wrapper(
        referenceData: *mut libc::c_double, referenceRows: libc::size_t, referenceCols: libc::size_t,
        queryData: *mut libc::c_double, queryRows: libc::size_t, queryCols: libc::size_t,
        bandwidth_val: libc::c_double, monte_carlo_tag: bool,
        prediction: *mut libc::c_double, predictionSize: libc::size_t,
    );
}

fn main() {
    // Define the data as arrays
    let mut ref_data = vec![3.0, 8.0, 9.0, 20.0];
    let mut q_data = vec![8.0, 9.0, 0.01];

    let ref_rows = 1;
    let ref_cols = ref_data.len();
    let q_rows = 1;
    let q_cols = q_data.len();

    let bandwidth = 1.0;
    let monte_carlo_tag = true;

    let mut prediction: Vec<f64> = vec![0.0; 3];



    // Call the mlpack_kde_wrapper function
    unsafe{
    mlpack_kde_wrapper(
        ref_data.as_mut_ptr(), ref_rows, ref_cols,
        q_data.as_mut_ptr(), q_rows, q_cols,
        bandwidth,
        monte_carlo_tag,
        prediction.as_mut_ptr(), prediction.len()
    );
}

    println!("Prediction: {:?}", prediction);
    println!("it is done");
}