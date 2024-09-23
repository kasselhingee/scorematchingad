#include "tapelltape.h"
#include "utils/PrintFor.hpp"

CppAD::ADFun<double> tapelltape(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               CppAD::ADFun<double> & llf, //the log likelihood function as a tape
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

  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> llfhigher;
  llfhigher = llf.base2ad();

  //tape relationship between x and log-likelihood
  CppAD::Independent(z, thetavar);  //for this tape, theta must be altered using new_dynamic
  if (verbose){
    Rcpp::Rcout << "thetavar is: " << thetavar.transpose() << std::endl;
    PrintForVec("\n thetavar is: ", thetavar);
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
    Rcpp::Rcout << "thetarecom is: " << thetarecom.transpose() << std::endl;
    PrintForVec("\n thetarecom is: ", thetarecom);
  }

  // range space vector
  veca1 y(1); // vector of ranges space variables
  veca1 u(0); //0 here because size dictated by fromM
  u = tran.fromM(z);
  y.setZero();
  llfhigher.new_dynamic(thetarecom);
  y[0] += llfhigher.Forward(0, u)[0];

  //get log determinant of fromM
  y[0] += tran.logdetJfromM(z);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  if (verbose){
    Rcpp::Rcout << "tape has " << tape.size_dyn_ind() << " independent dynamic parameters" << std::endl;
    Rcpp::Rcout << "tape requires vectors of length " << tape.Domain() << std::endl;
    Rcpp::Rcout << "tape returns vectors of length " << tape.Range() << std::endl;
  }
  return(tape);
}


CppAD::ADFun<double> ptapelltape(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     Rcpp::XPtr< CppAD::ADFun<double> > pllf, //the log likelihood function
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     ){
  CppAD::ADFun<double> out;
  out = tapelltape(z_ad,
                theta_ad,
                *pllf,
                tran,
                fixedtheta,
                verbose);
  return(out);
}

