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
  template <> mata1 as( SEXP inmat) {
    matd inemat;
    inemat = Rcpp::as<matd>(inmat);
    mata1 out(inemat.rows(), inemat.cols());
    for (long int i = 0; i<out.cols(); i++){
      for (long int j = 0; j<out.rows(); j++){
        out(j, i) = inemat(j, i);
      }
    }
    return(out);
  }

  // eigen matrix of a1type to SEXP
  template <> SEXP wrap(const mata1 &inmat){
    matd out(inmat.rows(), inmat.cols());
    for (long int i = 0; i<out.cols(); i++){
      for (long int j = 0; j<out.rows(); j++){
        out(j, i) = CppAD::Value(inmat(j, i));
      }
    }
    return(Rcpp::wrap(out)); //returns SEXP
  } 

}

# endif
