# ifndef APPROX_H
# define APPROX_H

#include <scorematchingad_forward.h>
#include <utils/pADFun.h>

//' @noRd
//' @describeIn evaltape_internal The value of a recorded function approximated by Taylor expansion.
//' Returns the approximate value of `pfun` at `x`.
//' @details
//' # taylorApprox_currentdynparam evaluates the tape without updating the dynamic parameter value 
//' Approximates the value of a `CppAD` tape at `x` using a Taylor approximation at `centre`. The dynamic parameters of the tape are set by `dynparam`.
//' @param centre For pTaylorApprox. A vector in the domain of the taped function to approximate the value at `x` from.
//' @param order For pTaylorApprox. The order of Taylor expansion to use.
// [[Rcpp::export]]
vecd taylorApprox_currentdynparam(pADFun & pfun,  //a tape with independent values that are points on the manifold (not the parameters)
		  vecd x,
                  vecd centre,
		  const size_t order);


//' @noRd
//' @describeIn evaltape_internal The value of a recorded function approximated by Taylor expansion.
//' Returns the approximate value of `pfun` at `x`.
// [[Rcpp::export]]
vecd taylorApprox(pADFun & pfun,  //a tape with independent values that are points on the manifold (not the parameters)
		  vecd x,
                  vecd centre,
                  vecd dynparam,
		  const size_t order);

# endif


