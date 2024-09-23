#ifndef tapelltape_h
#define tapelltape_h

#include "scorematchingad_forward.h" //includes RcppEigenForward.h
#include <Rcpp.h>
#include "utils/pADFun.h"

// define a function that converts a plain taped likelihood to:
// + have a different metric
// + fix some of the parameters
// On 19 Sep 2024: this will start as a bit of a mess because the tape comes with fixed dimension that must be matched by z and theta
// It is a duplicate nearly of tapellcpp()
CppAD::ADFun<double> tapelltape(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               CppAD::ADFun<double> & llf, //the log likelihood function as a tape
                               transform<a1type> & tran, //it seems pointer or references must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)
                               Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
                               bool verbose
                               );

//same as above but for pointers
// [[Rcpp::export]]
pADFun ptapelltape(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     Rcpp::XPtr< CppAD::ADFun<double> > pllf, //the log likelihood function
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     );

#endif
