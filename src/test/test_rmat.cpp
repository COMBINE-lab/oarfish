#include "mlpack_kde_wrapper.h"
#include <iostream>
#include <string>

int main(){

    std::vector<double> a = {1, 2, 3, 4, 5, 6};

    arma::mat queryVec(a.data(), 2, 3, false, false);
    

    queryVec.print();


    return 0;
}