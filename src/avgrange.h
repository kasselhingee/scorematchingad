#ifndef AVGRANGE_H
#define AVGRANGE_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>

//' @title Average Across Range of a Tape
//' @family tape builders
//' @param pfun An `Rcpp_ADFun` object.
//' @description Creates a `CppAD` tape that is the average of the returned values of `pfun`.
//' For creating this tape, the values of `pfun$dyntape` and `pfun$xtape` are used.
//' @return An `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun avgrange(pADFun & pfun);

#endif

