#include "utils/pADFun.h"

using namespace Rcpp;

RCPP_MODULE(cppad_module) {
    Rcpp::class_<pADFun>("ADFun")
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .property("size_order", &pADFun::size_order, "Number of Taylor coefficient orders, per variable,direction, currently calculated and stored")
        .property("Domain", &pADFun::size_order, "Dimension of domain space (i.e. length of independent vector).");
}

