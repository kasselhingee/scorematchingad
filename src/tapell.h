#ifndef tapell_h
#define tapell_h

# include <RcppEigen.h>
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "scorematchingad_forward.h"
# include "likelihoods/likelihoods.hpp"
# include "utils/PrintFor.hpp"
# include "utils/pADFun.h"

// declare a function that tapes a log likelihood
CppAD::ADFun<double> tapellcpp(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                               transform<a1type> & tran, //it seems points must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)
                               Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
                               bool verbose
                               );


//' @noRd
//' @title Tape of a log-density calculation 2
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
// [[Rcpp::export]]
pADFun ptapell2(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     Rcpp::XPtr<llPtr> llfXPtr, //the log likelihood function
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     );

//' @noRd
//' @title Get an XPtr to a named log-density function in source code of package
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to a `llPtr` object of the log-density function. Since `llPtr` is itself a pointer object, we have an XPtr pointing to a pointer that points to a function.
// [[Rcpp::export]]
Rcpp::XPtr<llPtr> getllptr(std::string llname);

//' @rdname tape_uld
//' @name tape_uld
//' @param name Name of an inbuilt function. See details.
//' @details
//' For `tape_uld_inbuilt()`, currently available unnormalised log-density functions are:
//'
//' ```{r, results = "asis", echo = FALSE}
//' cat(paste(" +", llnames), sep = "\n")
//' ```
//' @export
// [[Rcpp::export]]
pADFun tape_uld_inbuilt(std::string name, veca1 x, veca1 theta);

#endif
