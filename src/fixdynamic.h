#ifndef FIXDYNAMIC_H
#define FIXDYNAMIC_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>


// [[Rcpp::export]]
pADFun fixdynamic(pADFun & uld, // the unnormalised log density tape
                  veca1 theta, //new theta to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta); //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)

// Rembed the unnormalised log-density using a different metric, which means the embedding in Euclidean space changes (so the 'x' changes too)
// [[Rcpp::export]]
pADFun reembed(pADFun & uld, 
               transform<a1type> & tran); //it seems pointer or references must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)

#endif
