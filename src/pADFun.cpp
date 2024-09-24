#include "utils/pADFun.h"

using namespace Rcpp;

pADFun * XPtr2ppADFun(Rcpp::XPtr<CppAD::ADFun<double>> ptape){
  pADFun * out = new pADFun(ptape);
  return(out);
}

RCPP_MODULE(cppad_module) {
    Rcpp::class_<pADFun>("ADFun")
        .factory<Rcpp::XPtr<CppAD::ADFun<double>>>(XPtr2ppADFun)
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .property("size_order", &pADFun::size_order);
}

