#ifndef KEEPRANGE_H
#define KEEPRANGE_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>

//' @title Reduce Range Dimension of a Tape
//' @family tape builders
//' @param pfun An `Rcpp_ADFun` object.
//' @param keep Integers (lowest of 1, highest of `pfun$range`) specifying which elements of the range to keep. To keep all pass `keep = seq(1, pfun$range)`.
//' @description Retapes an existing `CppAD` tape omitting some of the returned elements.
//' For creating this tape, the values of `pfun$dyntape` and `pfun$xtape` are used.
//' @return An `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun keeprange(pADFun & pfun,
                 Eigen::Matrix<int, Eigen::Dynamic, 1> keep);


#endif
