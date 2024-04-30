#ifndef tapellcpp_h
#define tapellcpp_h

# include <RcppEigen.h>
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "scorematchingad.h"
# include "likelihoods/likelihoods.hpp"
# include "utils/PrintFor.hpp"

// define a function that tapes a log likelihood
CppAD::ADFun<double> tapellcpp(veca1 z, //data measurement in domain of llf 
                            veca1 theta, //theta parameter including fixed values
                               a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                               transform<a1type> & tran, //it seems pointer or references must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)
                               Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
                               bool verbose
                               ){
  if (theta.size() != fixedtheta.size()){
    Rcpp::stop("theta and fixedtheta must have the same length");
  }

  //separate fixed and variable theta
  veca1 thetavar(theta.size() - fixedtheta.sum());
  veca1 thetafxd(fixedtheta.sum());
  size_t idx_var(0);
  size_t idx_fxd(0);
  for (long int i=0;i<theta.size();i++){
    if (fixedtheta[i]){
      thetafxd[idx_fxd] = theta[i];
      idx_fxd += 1;
    } else {
      thetavar[idx_var] = theta[i];
      idx_var += 1;
    }
  }

  if (verbose){
    Rcpp::Rcout << "Fixing according to pattern: " << std::endl;
    for (long int i=0;i<fixedtheta.size();i++){
      Rcpp::Rcout << " " << fixedtheta[i];
    }
    Rcpp::Rcout << std::endl;

    Rcpp::Rcout << "Fixed theta is:";
    if (thetafxd.size() == 0){
      Rcpp::Rcout << " none" << std::endl;
    } else {
      for (long int i=0;i<thetafxd.size();i++){
        Rcpp::Rcout << " " << thetafxd[i];
      }
      Rcpp::Rcout << std::endl;
    }
  }

  //tape relationship between x and log-likelihood
  CppAD::Independent(z, thetavar);  //for this tape, theta must be altered using new_dynamic
  if (verbose){
    Rcpp::Rcout << "x is: " << z.transpose() << std::endl;
    PrintForVec("\n x is: ", z);
    Rcpp::Rcout << "dynamic theta elements are: " << thetavar.transpose() << std::endl;
    PrintForVec("\n dynamic theta elements are: ", thetavar);
  }

  //combine fixed and variable theta
  veca1 thetarecom(theta.size());
  idx_var = 0;
  for (long int i=0;i<theta.size();i++){
    if (fixedtheta[i]){
      thetarecom[i] = theta[i];
    } else {
      thetarecom[i] = thetavar[idx_var];
      idx_var += 1;
    }
  }
  if (verbose){
    Rcpp::Rcout << "full theta is: " << thetarecom.transpose() << std::endl;
    PrintForVec("\n full theta is: ", thetarecom);
  }

  // range space vector
  veca1 y(1); // vector of ranges space variables
  y.setZero();
  y[0] += llf(z, thetarecom);

  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  if (verbose){
    Rcpp::Rcout << "tape has " << tape.size_dyn_ind() << " independent dynamic parameters" << std::endl;
    Rcpp::Rcout << "tape requires vectors of length " << tape.Domain() << std::endl;
    Rcpp::Rcout << "tape returns vectors of length " << tape.Range() << std::endl;
  }
  return(tape);
}


#endif
