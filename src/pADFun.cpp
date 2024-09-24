#include "utils/pADFun.h"

//small function to run that will force Module to initialise as
//I'm seeing strange loading when using both the manifolds module
//and the cppad_module
void empty(){}


RCPP_MODULE(cppad_module) {
    Rcpp::function("empty", &empty);

    Rcpp::class_<pADFun>("ADFun")
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .property("size_order", &pADFun::size_order, "Number of Taylor coefficient orders, per variable,direction, currently calculated and stored")
        .property("Domain", &pADFun::Domain, "Dimension of domain space (i.e. length of independent variables vector).")
        .property("Range", &pADFun::Range, "Dimension of range space.")
        .property("size_dyn_ind", &pADFun::size_dyn_ind, "Number of (independent) dynamic parameters.")
        
        .method("Parameter", &pADFun::Parameter, "Returns true if the ith component of the range space corresponds to a 'Parameter' and is thus constant.");
//        .method("Forward", &pADFun::Forward, "Forward mode evaluation.");


}

