#ifndef PADFUN_H
#define PADFUN_H

//both Rcpp::wrap and Rcpp::as use copy operations, which ADFun explicitly does not have (the default is 'deleted'). [UNLESS IT WAS JUST ME ACCIDENTALLY USING COPY]?
//So exposing pointers of ADFun objects to make interrogation possible
//Hopefully the wrapping for Rcpp::XPtr won't override anything I try here


#include <scorematchingad_forward.h>
#include <Rcpp.h>

class pADFun {
private:
  Rcpp::XPtr < CppAD::ADFun<double> > ptr; //Using Rcpp::XPtr here for automatic management of memory. I could using something Cpp specific for probably faster results

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

// Constructor taking XPtr - and copying it
pADFun(Rcpp::XPtr<CppAD::ADFun<double>> p) : ptr(p) {}

size_t size_order() const { return ptr->size_order(); }
size_t Domain() const { return ptr->Domain(); }
size_t Range() const { return ptr->Domain(); }
size_t size_var() const { return ptr->size_var(); }
size_t size_par() const { return ptr->size_par(); }
size_t size_op() const { return ptr->size_op(); }
size_t size_op_arg() const { return ptr->size_op_arg(); }
size_t size_dyn_ind() const { return ptr->size_dyn_ind(); }
bool Parameter(size_t i) const { return ptr->Parameter(i); }

void new_dynamic(const vecd & dynamic) {ptr->new_dynamic(dynamic);}
vecd Forward(size_t q, const vecd & xq) {return ptr->Forward(q, xq);}
vecd Jacobian(const vecd & x) {return ptr->Jacobian(x);}
vecd Hessiani(const vecd & x, size_t i) {return ptr->Hessian(x, i);}
vecd Hessian0(const vecd & x) {return ptr->Hessian(x, 0);}
vecd Hessianw(const vecd & x, const vecd & w) {return ptr->Hessian(x, w);}

//methods to simplify use
vecd eval(const vecd & x, const vecd & dyn) {
  ptr->new_dynamic(dyn);
  return ptr->Forward(0, x);}
vecd Jac(const vecd & x, const vecd & dyn) {
  ptr->new_dynamic(dyn);
  return ptr->Jacobian(x);}
vecd Hes(const vecd & x, const vecd & dyn) {
  ptr->new_dynamic(dyn);
  return ptr->Hessian(x, 0);}



};

RCPP_EXPOSED_CLASS(pADFun)

#endif
