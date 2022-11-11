// definitions for Rcpp::wrap and Rcpp::as for various data types
namespace Rcpp {
  // from an SEXP to an eigen vector
  template <> vecd as( SEXP invec) {
    Rcpp::NumericVector insvec(invec);
    vecd out(insvec.size());
    for (long int i = 0; i<out.size(); i++){
      out[i] = insvec[i];
    }
    return(out);
  }

  // eigen vector using the std::vector (called svecd for short in my code)
  template <> SEXP wrap(const vecd &invec){
    svecd out(invec.size());
    for (long int i = 0; i<invec.size(); i++){
      out[i] = invec[i];
    }
    return(Rcpp::wrap(out)); //returns SEXP
  } 
}
