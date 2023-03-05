#ifndef util
#define util

#include "scorecompdir_types.h"
#include "mycpp/manifolds.hpp"
#include "mycpp/approx.hpp"
#include <Rcpp.h>

//' @noRd
//' @title Apply to `toM` function of a manifold object
//' @description Apply the `toM` function of a manifold object.
//' @param tran An Rcpp_transform_ad object.
//' @param u A vector to be transformed to the manifold via `toM`.
//' @return A vector on the manifold.
// [[Rcpp::export]]
veca1 ptoM(transform_a1type & tran, veca1 u_ad);


//' @describeIn evaltape_internal The value of a recorded function approximated by Taylor expansion.
//' Returns the approximate value of `pfun` at `x`.
//' @details
//' # pTaylorApprox 
//' Approximates the value of a `CppAD` tape at `x` using a Taylor approximation at `centre`. The dynamic parameters of the tape are set by `dynparam`.
//' @param centre For pTaylorApprox. A vector in the domain of the taped function to approximate the value at `x` from.
//' @param order For pTaylorApprox. The order of Taylor expansion to use.
//' @export
// [[Rcpp::export]]
vecd pTaylorApprox(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                     vecd x, vecd centre,
                     vecd dynparam, size_t order);

// A function that creates a tape of fromM
CppAD::ADFun<double> tapefromM(veca1 z,
                               transform<a1type> & tran);

// And its version for R
// for calculating the determinant of the transform to a manifold
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > ptapefromM(veca1 z,
                               transform<a1type> & tran);



#endif
