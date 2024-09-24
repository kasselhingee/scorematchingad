#include "utils/pADFun.h"

//small function to run that will force Module to initialise as
//I'm seeing strange loading when using both the manifolds module
//and the cppad_module
void empty(){}


RCPP_MODULE(cppad_module) {
    function("empty", &empty);

    Rcpp::class_<pADFun>("ADFun")
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .property("size_order", &pADFun::size_order, "Number of Taylor coefficient orders, per variable,direction, currently calculated and stored")
        .property("Domain", &pADFun::Domain, "Dimension of domain space (i.e. length of independent vector).");
}

