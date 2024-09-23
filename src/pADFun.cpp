#include "utils/pADFun.h"



RCPP_MODULE(cppad_module) {
    Rcpp::class_<pADFun>("ADFun")
        .constructor()
        .property("size_order", &pADFun::size_order);
}

