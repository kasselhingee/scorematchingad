#include "utils/pADFun.h"

using namespace Rcpp;

RCPP_MODULE(cppad_module) {
    Rcpp::class_<pADFun>("ADFun")
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .property("size_order", &pADFun::size_order);
}

