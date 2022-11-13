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
}
