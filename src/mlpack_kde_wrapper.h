#ifndef MLPACK_KDE_WRAPPER_H
#define MLPACK_KDE_WRAPPER_H

#include <mlpack.hpp>

extern "C" {
    enum class ErrorCode {
        Success = 0,
        CDFFailed = 1,
        PointerFailed = 2,
        PointerSizeFailed = 3,
        InvalidInput = 4,
    };

    int mlpack_kde_wrapper(double* referenceData, size_t referenceRows, size_t referenceCols,
                            double* queryData, size_t queryRows, size_t queryCols,
                            double bandwidth_val, bool monte_carlo_tag,
                            double* prediction, size_t predictionSize);
}

#endif // MLPACK_KDE_WRAPPER_H
