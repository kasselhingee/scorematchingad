#ifndef tapell_h
#define tapell_h

# include <RcppEigen.h>
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "scorematchingad.h"
# include "likelihoods/likelihoods.hpp"
# include "utils/PrintFor.hpp"

// declare a function that tapes a log likelihood
CppAD::ADFun<double> tapellcpp(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                               transform<a1type> & tran, //it seems points must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)
                               Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
                               bool verbose
                               );

//' @noRd
//' @title Tape of a log-likelihood calculation
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > ptapell(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     std::string llname,
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     );

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

#endif
