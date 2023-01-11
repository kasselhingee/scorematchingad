# ifndef MYCPP_WRAPAS
# define MYCPP_WRAPAS

# include "../scorecompdir_types.h"
# include <RcppEigen.h> //for Rcpp::as to Eigen entities
# include <Rcpp.h>
# include <cppad/cppad.hpp>

// definitions for Rcpp::wrap and Rcpp::as for various data types
namespace Rcpp {
  // from an SEXP to an eigen vector if a1type
  template <> veca1 as( SEXP invec) {
    vecd inevec;
    inevec = Rcpp::as<vecd>(invec);
    veca1 out(inevec.size());
    for (long int i = 0; i<out.size(); i++){
      out[i] = inevec[i];
    }
    return(out);
  }

  // eigen vector of a1type to SEXP
  template <> SEXP wrap(const veca1 &invec){
    vecd out(invec.size());
    for (long int i=0; i<out.size(); i++){
      out[i] = CppAD::Value(invec[i]);
    }
    return(Rcpp::wrap(out)); //returns SEXP
  } 



  // from an SEXP to an eigen matrix of a1type
  //template <> mata1 as( SEXP inmat) {
  //  Rcpp::NumericMatrix insmat(inmat);
  //  mata1 out(insmat.size());
  //  for (long int i = 0; i<out.size(); i++){
  //    out[i] = insmat[i];
  //  }
  //  return(out);
  //}

  // eigen matrix of a1type to SEXP
  //template <> SEXP wrap(const mata1 &inmat){
  //  svecd out(inmat.size());
  //  for (long int i=0; i<out.size(); i++){
  //    out[i] = CppAD::Value(inmat[i]);
  //  }
  //  return(Rcpp::wrap(out)); //returns SEXP
  //} 

}

# endif
