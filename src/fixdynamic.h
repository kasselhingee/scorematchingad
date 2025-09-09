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

//' @title Build unnormalised log-density in new embedding
//' @param pfun An `Rcpp_ADFun` object.
//' @param tran A transform object.
//' @description Build a tape of the unnormalised log-density on a different embedding of the manifold.
//' @details
//' When the embedding of a manifold changes according to `tran`, then the Riemannian metric on the manifold changes and so the uniform measure of the manifold also changes.
//' This change is accounted for using the `logdetJfromM()` property of transform objects.
//'
//' The returned tape has an independent variable that is on the newly embedded manifold, and the value used for taping was `tran$toM(pfun$xtape)`.
//' @export
// [[Rcpp::export]]
pADFun reembed(pADFun & uld, 
               transform<a1type> & tran); //it seems pointer or references must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)

#endif
