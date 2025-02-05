#ifndef FIXDYNAMIC_H
#define FIXDYNAMIC_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>

//' @title Fix Dynamic Parameters of a Tape
//' @family tape builders
//' @param pfun An `Rcpp_ADFun` object.
//' @param theta A numerical vector specifying the value of all dynamic parameters of `pfun`. Some of these will be fixed according to `isfixed`, the remainder will remain dynamic.
//' @param isfixed A boolean vector same length as `theta`. `TRUE` values are fixed at the value of `theta`, `FALSE` values are left dynamic.
//' @description Retapes an existing `CppAD` tape but with some of original dynamic parameters fixed to specified values.
//' For creating this tape, the values of `pfun$xtape` is used.
//' @return An `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun fixdynamic(pADFun & pfun, // the unnormalised log density tape
                  veca1 theta, //new theta to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> isfixed); //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)

// Rembed the unnormalised log-density using a different metric, which means the embedding in Euclidean space changes (so the 'x' changes too)
// [[Rcpp::export]]
pADFun reembed(pADFun & uld, 
               transform<a1type> & tran); //it seems pointer or references must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)

#endif
