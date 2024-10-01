// A function that calls the CppAD abort_recording() function from within R

#include <scorematchingad_forward.h>

// [[Rcpp::export]]
void abort_recording(){
  a1type::abort_recording();
  a2type::abort_recording();
}

