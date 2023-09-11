#include "mlpack_kde_wrapper.h"
#include <iostream>
#include <string>

int mlpack_kde_wrapper(
    double* referenceData, size_t referenceRows, size_t referenceCols,
    double* queryData, size_t queryRows, size_t queryCols,
    double bandwidth_val, bool monte_carlo_tag,
    double* prediction, size_t predictionSize
) {
    try {
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

        // cumulative summation function
        arma::vec cdf_estimation = arma::cumsum(estimations);

        // Normalizing cumulative distribution function (CDF)
        cdf_estimation /= cdf_estimation.back();

        if (cdf_estimation.is_empty() || !cdf_estimation.is_finite()) {
            return static_cast<int>(ErrorCode::CDFFailed);
        }

        // Copy cdf_estimation to prediction
        if (prediction == nullptr || cdf_estimation.memptr() == nullptr) {
            return static_cast<int>(ErrorCode::PointerFailed);
        }

        if (predictionSize != cdf_estimation.n_elem) {
            return static_cast<int>(ErrorCode::PointerSizeFailed);
        }

        std::memcpy(prediction, cdf_estimation.memptr(), sizeof(double) * predictionSize);

        return static_cast<int>(ErrorCode::Success);

    } catch (const std::exception& e) {
        return static_cast<int>(ErrorCode::InvalidInput); 
    }
}

