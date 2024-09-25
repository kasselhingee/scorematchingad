#include "utils/pADFun.h"
#include <RcppEigen.h> //need RcppEigen fully included here for Rcpp to find vecd wrap and as for items below

//small function to run that will force Module to initialise as
//I'm seeing strange loading when using both the manifolds module
//and the cppad_module
void empty(){}


RCPP_MODULE(cppad_module) {
    Rcpp::function("empty", &empty);

    Rcpp::class_<pADFun>("ADFun")
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .property("size_order", &pADFun::size_order, "Number of Taylor coefficient orders, per variable,direction, currently calculated and stored")
        .property("domain", &pADFun::Domain, "Dimension of domain space (i.e. length of independent variables vector).")
        .property("range", &pADFun::Range, "Dimension of range space.")
        .property("size_dyn_ind", &pADFun::size_dyn_ind, "Number of (independent) dynamic parameters.")

        .method("new_dynamic", &pADFun::new_dynamic, "Specify new values of the dynamic parameters.")
        .method("forward", &pADFun::Forward, "Forward mode evaluation.")
        .method("Jacobian", &pADFun::Jacobian, "Evaluate Jacobian.")
        .method("Hessiani", &pADFun::Hessiani, "Evaluate Hessian for ith element of range i = 0, 1, ...")
        .method("Hessian0", &pADFun::Hessian0, "Evaluate Hessian for first element of range.")
        .method("Hessianw", &pADFun::Hessianw, "Evaluate Hessian for weighted sum of range.")
        
        .method("eval", &pADFun::eval, "Evaluation with new dynamic.")
        .method("Jac", &pADFun::Jac, "Jacobian with new dynamic.")
        .method("Hes", &pADFun::Hes, "Hessian with new dynamic.")
        
        .method("parameter", &pADFun::Parameter, "Returns true if the ith component of the range space corresponds to a 'Parameter' and is thus constant.");


}

