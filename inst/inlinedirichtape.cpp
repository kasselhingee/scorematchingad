#include <scorematchingad_forward.h>
#include <Rcpp.h>
#include <scorematchingad.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(scorematchingad)]]


// [[Rcpp::export]]
a1type dirich(const veca1 &u, const veca1 &beta) {
  size_t d  = u.size();
  a1type y(0.);  // initialize summation at 0
  for(size_t i = 0; i < d; i++)
  {   y   += beta[i] * CppAD::log(u[i]);
  }
  return y;
}

// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > tapedirich(veca1 & u, veca1 & beta){
  CppAD::ADFun<double> * ptape = new CppAD::ADFun<double>; 
  CppAD::Independent(u, beta);
  veca1 y(1);
  y(0) = dirich(u, beta);
  ptape -> Dependent(u, y);
  Rcpp::XPtr< CppAD::ADFun<double> > pout(ptape, true);
  return(pout);
}
