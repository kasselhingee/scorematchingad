#ifndef FIXDYNAMIC_H
#define FIXDYNAMIC_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>


// [[Rcpp::export]]
pADFun fixdynamic(pADFun & uld, // the unnormalised log density tape
                  veca1 theta, //new theta to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta); //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)


#endif
