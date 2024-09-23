#ifndef PADFUN_H
#define PADFUN_H

//both Rcpp::wrap and Rcpp::as use copy operations, which ADFun explicitly does not have (the default is 'deleted'). [UNLESS IT WAS JUST ME ACCIDENTALLY USING COPY]?
//So exposing pointers of ADFun objects to make interrogation possible
//Hopefully the wrapping for Rcpp::XPtr won't override anything I try here


#include <scorematchingad_forward.h>
#include <Rcpp.h>

class pADFun {
private:
  Rcpp::XPtr < CppAD::ADFun<double> > ptr;

//MOVES the tape to a spot with memory management. The previous version of the tape will NOT be available because it this is a MOVE. Move needed because copy operation of ADFun not allowed
Rcpp::XPtr < CppAD::ADFun<double> > movetoXPtr(CppAD::ADFun<double> & tape){
  CppAD::ADFun<double> * out = new CppAD::ADFun<double>(std::move(tape));//reserve memory for a tape
  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

public:
// Default constructor
pADFun() : ptr(Rcpp::XPtr< CppAD::ADFun<double> >(nullptr, false)) {}

// Constructor taking tape
pADFun(CppAD::ADFun<double> & tape) : ptr(movetoXPtr(tape)) {}

// Constructor taking pADFundouble
pADFun(const Rcpp::XPtr<CppAD::ADFun<double>> & p) : ptr(p) {}

size_t size_order() const {
   return ptr->size_order();
}

};

RCPP_EXPOSED_CLASS_NODECL(pADFun)

#endif
