#include "mlpack_kde_wrapper.h"
#include <iostream>
#include <string>

void mlpack_kde_wrapper(
    double* referenceData, size_t referenceRows, size_t referenceCols,
    double* queryData, size_t queryRows, size_t queryCols,
    double bandwidth_val, bool monte_carlo_tag,
    double* prediction, size_t predictionSize
) {
    // Convert inputs to appropriate MLPack data types
    arma::mat referenceVec(referenceData, referenceRows, referenceCols, false, false);
    arma::mat queryVec(queryData, queryRows, queryCols, false, false);
    //arma::vec predictionVec(predictionSize, arma::fill::zeros);
    mlpack::KDE<> kde;  // Create a KDE object
    kde.Kernel() = mlpack::GaussianKernel();  // Use Gaussian kernel
    kde.Kernel().Bandwidth(bandwidth_val);  // Set bandwidth
    kde.RelativeError(0.05);
    kde.AbsoluteError(0.0);

    kde.MonteCarlo() = monte_carlo_tag; // Set Monte Carlo
    if (monte_carlo_tag) {
          // Set other Monte Carlo parameters here
          kde.MCProb(0.95);
          kde.MCInitialSampleSize() = 200;
          kde.MCEntryCoef(3.0);
          kde.MCBreakCoef(0.4);
    }

    // Train the KDE model
    kde.Train(referenceVec);

    // Perform evaluation
    arma::vec estimations;
    kde.Evaluate(queryVec, estimations);

    // Copy estimations to prediction
    std::memcpy(prediction, estimations.memptr(), sizeof(double) * predictionSize);
}

#include <vector>

int main(){
        // Define the data as arrays
    double ref_data[] = {3, 8, 9, 20};
    double q_data[] = {8, 9};

    // Create vectors from the arrays
    std::vector<double> ref(ref_data, ref_data + sizeof(ref_data) / sizeof(double));
    std::vector<double> q(q_data, q_data + sizeof(q_data) / sizeof(double));

    size_t ref_c = ref.size();
    size_t ref_r = 1;
    size_t q_c = q.size();
    size_t q_r = 1;

    double bandwidth = 1.0;
    bool monte_carlo_tag = true;

    // Create a vector to hold the prediction
    std::vector<double> prediction(2);

    // Call the mlpack_kde_wrapper function
    mlpack_kde_wrapper(
        ref.data(), ref_r, ref_c,
        q.data(), q_r, q_c,
        bandwidth, monte_carlo_tag,
        prediction.data(), prediction.size()
    );

    std::cout << "Prediction: ";
    for (double value : prediction) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::cout << "it is done" << std::endl;
    
    return 0;
}


fn main() {
    // Define the data as arrays
    let mut ref_data = vec![3.0, 8.0, 9.0, 20.0];
    let mut q_data = vec![8.0, 9.0];

    let ref_rows = 1;
    let ref_cols = ref_data.len();
    let q_rows = 1;
    let q_cols = q_data.len();

    let bandwidth = 1.0;
    let monte_carlo_tag = true;

    let mut prediction: Vec<f64> = vec![0.0; 2];



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
