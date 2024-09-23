#ifndef PADFUN_H
#define PADFUN_H

#include <scorematchingad_forward.h>
#include <Rcpp.h>

class pADFun {
private:
  Rcpp::XPtr < CppAD::ADFun<double> > ptr;

  static Rcpp::XPtr < CppAD::ADFun<double> > movetoXPtr(CppAD::ADFun<double> & tape);

public:
  // Default constructor
  pADFun();

  // Constructor taking tape
  pADFun(CppAD::ADFun<double> & tape);

  // Constructor taking pADFundouble
  pADFun(const Rcpp::XPtr<CppAD::ADFun<double>> & p);

  size_t size_order() const;

};

RCPP_EXPOSED_CLASS_NODECL(pADFun)

#endif
