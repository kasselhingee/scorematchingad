# ifndef MYCPP_WRAPAS
# define MYCPP_WRAPAS

// definitions for Rcpp::wrap and Rcpp::as for various data types
namespace Rcpp {
  // from an SEXP to an eigen vector if a1type
  template <> veca1 as( SEXP invec) {
    Rcpp::NumericVector insvec(invec);
    veca1 out(insvec.size());
    for (long int i = 0; i<out.size(); i++){
      out[i] = insvec[i];
    }
    return(out);
  }

  // eigen vector of a1type to SEXP
  template <> SEXP wrap(const veca1 &invec){
    svecd out(invec.size());
    for (long int i=0; i<out.size(); i++){
      out[i] = CppAD::Value(invec[i]);
    }
    return(Rcpp::wrap(out)); //returns SEXP
  } 
}

# endif
