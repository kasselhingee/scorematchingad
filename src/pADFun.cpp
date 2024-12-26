#include <scorematchingad_forward.h>
#include "utils/pADFun.h"
#include <RcppEigen.h> //need RcppEigen fully included here for Rcpp to find vecd wrap and as for items below

// LOOK TO R/ADFun.R FOR HELP ON THESE OBJECTS

RCPP_MODULE(ADFun) {
    Rcpp::class_<pADFun>("ADFun")
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>>()
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>, veca1, veca1>()
        .constructor<Rcpp::XPtr<CppAD::ADFun<double>>, veca1, veca1, std::string>() //don't pass by reference - Rcpp generated errors as compile time
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
        .method("set_check_for_nan", &pADFun::set_check_for_nan, "Set the tape to check for nan values when C++ debugging enabled.")
        .method("get_check_for_nan", &pADFun::get_check_for_nan, "Return whethe the tape will check for nan values when C++ debugging enabled.")
        
        .method("eval", &pADFun::eval, "Evaluation with new dynamic.")
        .method("Jac", &pADFun::Jac, "Jacobian with new dynamic.")
        .method("Hes", &pADFun::Hes, "Hessian with new dynamic.")

        .field("name", &pADFun::name, "An optional name for the tape.")
        .field_readonly("xtape", &pADFun::xtape, "The values of the independent variables used for taping.")
        .field_readonly("dyntape", &pADFun::dyntape, "The values of the dynamic variables used for taping.")
        
        .method("parameter", &pADFun::Parameter, "Returns true if the ith component of the range space corresponds to a 'Parameter' and is thus constant.");


}

