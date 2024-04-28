#ifndef tapell_h
#define tapell_h

# include <RcppEigen.h>
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "scorematchingad.h"
# include "likelihoods/likelihoods.hpp"
# include "utils/PrintFor.hpp"
# include "tapellcpp.h"

//' @noRd
//' @title Tape of a log-likelihood calculation 2
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > ptapell2(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     Rcpp::XPtr<llPtr> llfXPtr, //the log likelihood function
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     );

//' @noRd
//' @title Get an XPtr to a named log-likelihood function in source code of package
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to a `llPtr` object of the log-likelihood function. Since `llPtr` is itself a pointer object, we have an XPtr pointing to a pointer that points to a function.
// [[Rcpp::export]]
Rcpp::XPtr<llPtr> getllptr(std::string llname);

//' @name evalll
//' @title Evaluate a log-likelihood function
//' @param ll A compiled log-likelihood function created by [`customll()`].
// ( ll is an XPtr to a llPtr object that points to a log-likelihood function )
//' @param u A vector of measurements for an individual
//' @param theta A vector of parameters
//' @return The value of the log-likelihood at `u` with parameters `theta`.
//' @export
// [[Rcpp::export]]
a1type evalll(Rcpp::XPtr<llPtr> ll, const veca1& u, const veca1& theta);

#endif
