#ifndef exposemanifold 
#define exposemanifold

#include "scorecompdir_types.h"
#include "mycpp/mantrans.hpp"
#include <Rcpp.h>


////////////// Create Pointers to Manifold Objects ///////////////
//in R store a pointer to the ADFun object
//' @noRd
//' @title Generate manifold with transformation object
//' @param manifoldname The name of the manifold to transform to. Either 'sphere' or 'simplex'
//' @return An RCpp::XPtr object pointing to the C++ manifold object
//' @details
//'  + "sphere" for square-root transformation from the simplex to the positive orthant of the sphere
//'  + "simplex" for the simplex without any transformation.
//'  + "Ralr" for the additive log-ratio transformation from the simplex to Euclidean space, using the final component of vectors in the denominator of the ratio.
//'  + "Snative" for the sphere without any transformation
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< manifold<a1type> > pmanifold(std::string manifoldname);

//' @noRd
//' @title Test a manifold object
//' @description A lightweight test of a manifold object.
//' Its main benefit is to force compilation of templated functions for the manifold,
//' and to print results to standard output.
//' @param pman An Rcpp::XPtr to a manifold object. Created by `pmanifold()`
//' @return An integer. 0 if the testable parts pass.
// [[Rcpp::export]]
int testmanifold(Rcpp::XPtr< manifold<a1type> > pman, veca1 u_ad);


#endif

