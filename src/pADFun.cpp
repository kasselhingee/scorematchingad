#include "utils/pADFun.h"

//both Rcpp::wrap and Rcpp::as use copy operations, which ADFun explicitly does not have (the default is 'deleted').
//So exposing pointers of ADFun objects to make interrogation possible
//Hopefully the wrapping for Rcpp::XPtr won't override anything I try here

//MOVES the tape to a spot with memory management. The previous version of the tape will NOT be available because it this is a MOVE. Move needed because copy operation of ADFun not allowed
Rcpp::XPtr < CppAD::ADFun<double> > pADFun::movetoXPtr(CppAD::ADFun<double> & tape){
  CppAD::ADFun<double> * out = new CppAD::ADFun<double>(std::move(tape));//reserve memory for a tape
  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

// Default constructor
pADFun::pADFun() : ptr(Rcpp::XPtr< CppAD::ADFun<double> >(nullptr, false)) {}

// Constructor taking tape
pADFun::pADFun(CppAD::ADFun<double> & tape) : ptr(movetoXPtr(tape)) {}

// Constructor taking pADFundouble
pADFun::pADFun(const Rcpp::XPtr<CppAD::ADFun<double>> & p) : ptr(p) {}

size_t pADFun::size_order() const {
   return ptr->size_order();
}


RCPP_MODULE(cppad_module) {
    Rcpp::class_<pADFun>("ADFun")
        .constructor()
        .property("size_order", &pADFun::size_order);
}

