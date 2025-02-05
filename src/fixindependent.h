#ifndef FIXINDEPENDENT_H
#define FIXINDEPENDENT_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>

//' @title Fix Independent Arguments of a Tape
//' @family tape builders
//' @param pfun An `Rcpp_ADFun` object.
//' @param x A numerical vector specifying the value of all independent arguments of `pfun`. Some of these will be fixed according to `isfixed`, the remainder will remain as independent arguments.
//' @param isfixed A boolean vector same length as `x`. `TRUE` values are fixed at the value of `x`, `FALSE` values are left as independent arguments.
//' @description Retapes an existing `CppAD` tape but with some of original independent arguments fixed to specified values.
//' For creating this tape, the values of `pfun$dyntape` are used.
//' @return An `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun fixindependent(pADFun & pfun, // the unnormalised log density tape
                  veca1 x, //new x to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> isfixed); //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)


#endif
