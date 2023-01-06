#ifndef util
#define util

#include "scorecompdir_types.h"
#include "mycpp/mantrans.hpp"
#include "mycpp/approx.hpp"
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

//' @noRd
//' @title Apply to `toM` function of a manifold object
//' @description Apply the `toM` function of a manifold object.
//' @param pman An Rcpp::XPtr to a manifold object. Created by `pmanifold()`.
//' @param u A vector to be transformed to the manifold via `toM`.
//' @return A vector on the manifold.
// [[Rcpp::export]]
veca1 ptoM(Rcpp::XPtr< manifold<a1type> > pman, veca1 u_ad);


//' @title The value of a recorded function approximated by Taylor expansion
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with independent values that are the points to be differentiated with
//' @param u A vector in the domain of the taped function.
//' @param centre A vector in the domain of the taped function to approximate the value at `u` from.
//' @param dynparam a vector of the dynamic parameters
//' @param order The order of Taylor expansion to use.
//' @description Approximates the value of a `CppAD` tape at `u` using a Taylor approximation at `centre`. The dynamic parameters of the tape are set by `dynparam`.
//' @return The approximate value of pfun
//' @export
// [[Rcpp::export]]
vecd pTaylorApprox(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                     vecd u, vecd centre,
                     vecd dynparam, size_t order);

// A function that creates a tape of fromM
CppAD::ADFun<double> tapefromM(veca1 z,
                               manifold<a1type> *pman);

// And its version for R
// for calculating the determinant of the transform to a manifold
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > ptapefromM(veca1 z,
                               Rcpp::XPtr<manifold<a1type> > pman);



#endif
