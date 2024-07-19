//use pyo3::prelude::*;
//use pyo3::types::PyDict;
//use numpy::{PyArray1, PyArray2};
//use ndarray::{Array1, Array2, ArrayView2};
//use numpy::IntoPyArray;
//use pyo3::prepare_freethreaded_python;
use crate::util::oarfish_types::InMemoryAlignmentStore;
use kders::kde::GridDimensions;

pub fn kde_computation(
    data: &Vec<f64>,
    weights: &Vec<f64>,
    store: &mut InMemoryAlignmentStore,
) -> anyhow::Result<()> {
    
    let n = data.len() / 2;
    let (mut max_x, mut max_y) = (0_f64, 0_f64);
    for i in 0..n {
        max_x = data[2 * i].max(max_x);
        max_y = data[(2*i) + 1].max(max_y);
    }

    let gd = kders::kde::GridDimensions {
        width: max_x as usize + 1,
        height: max_y as usize + 1,
    };

    let kernel_bandwidth = 1_f64;
    let bin_width = 10_usize;

    let mut grid = kders::kde::KDEGrid::new(gd, bin_width, Some(kernel_bandwidth));
    for (chunk, w) in data.chunks(2).zip(weights.iter()) {
        grid.add_observation(chunk[0] as usize, chunk[1] as usize, *w);
    }

    let density = grid.get_nearest_neighbor_kde()?;

    let mut lookups = Vec::<f64>::with_capacity(n);
    for chunk in data.chunks(2) {
        lookups.push(density[(chunk[0] as usize, chunk[1] as usize)]);
    }
        
    store.kde_prob = lookups;

    Ok(())
    
}



//old version using python library
//pub fn kde_computation(
//    data: &Vec<f64>,
//    weight: &Vec<f64>,
//    store: &mut InMemoryAlignmentStore,
//) -> PyResult<()> {
//    // Initialize the Python interpreter
//    prepare_freethreaded_python();
//    
//    Python::with_gil(|py| {
//        // Import the Python module and function
//        let module = PyModule::from_code(
//            py,
//            include_str!("/mnt/scratch4/zahra/oarfish_write_change/oarfish/src/util/2D_kde_function.py"),
//            "2D_kde_function",
//            "/mnt/scratch4/zahra/oarfish_write_change/oarfish/src/util/2D_kde_function.py",
//        )?;
//        let function = module.getattr("calculate_kde")?;
//        println!("data length is: {:?}", data.len());
//        println!("weight length is: {:?}", weight.len());
//
//        // Example 2D data
//        let data = Array2::from_shape_vec(((data.len() / 2), 2), data.to_vec()).unwrap();
//
//        // Convert data to PyArray2
//        let data_py: &PyArray2<f64> = data.view().to_owned().into_pyarray(py);
//
//        // Dummy weights for demonstration purposes
//        let weights = Array2::from_shape_vec((weight.len(), 1), weight.to_vec()).unwrap();
//
//        // Convert weights to PyArray2
//        let weights_py: &PyArray2<f64> = weights.view().to_owned().into_pyarray(py);
//
//        // Call the Python function
//        println!("it is before python function");
//        let result: Py<PyArray1<f64>> = function.call1((data_py, weights_py))?.extract()?;
//        //println!("result length: {:?}", result.as_ref().len());
//        println!("it is after python function");
//        println!("result is: {:?}", result);
//
//        // Convert the result to Rust ndarray
//        let result_array: &PyArray1<f64> = result.extract(py)?;
//        let result_slice = unsafe {
//            result_array.as_slice().expect("Failed to convert PyArray to slice")
//        };
//        let result_array_1d: Array1<f64> = Array1::from_shape_vec((result_slice.len(),), result_slice.to_vec()).expect("Failed to create Array1 from slice");
//        println!("data length: {:?}", (data.len() / 2));
//        println!("data length: {:?}", (weight.len()));
//        println!("prob length: {:?}", (result_array_1d.len()));
//        //println!("{:?}", result_array); // Print the result
//        store.kde_prob = result_array_1d.to_vec();
//
//        Ok(())
//    })
//}