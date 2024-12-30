# ifndef PADFUN_RETAPE
# define PADFUN_RETAPE

// things for manipulating and evaluating `Rcpp_ADFun` objects from R

//for content that is Rcpp specific
#include "scorematchingad_forward.h"
#include <Rcpp.h>
#include <utils/pADFun.h>

//' @title Tape the Jacobian of CppAD Tape
//' @family tape builders
//' @param pfun An `Rcpp_ADFun` object.
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters.
//' @description Creates a tape of the Jacobian of function taped by `CppAD`.
//' When the function returns a real value (as is the case for densities and the score matching objective) the Jacobian is equivalent to the gradient.
//' The `x` vector is used as the value to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say `$eval()`, the resultant vector contains the Jacobian in long format (see <https://cppad.readthedocs.io/latest/Jacobian.html>).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{m}-dimensional space, then
//' the first \eqn{n} elements of vector is the gradient of the first component of function output.
//' The next \eqn{n} elements of the vector is the gradient of the second component of the function output.
//' The Jacobian as a matrix, could then be obtained by [`as.matrix()`] with `byrow = TRUE` and `ncol = n`.
//' @return A `Rcpp::XPtr` to a `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun  tapeJacobian(pADFun & pfun);

//' @title Tape the Hessian of a CppAD Tape
//' @family tape builders
//' @inheritParams pTapeJacobian
//' @description Creates a tape of the Hessian of a function taped by CppAD.
//' The taped function represented by `pfun` must be scalar-valued (i.e. a vector of length 1).
//' The `x` vector and `dynparam` are used as the values to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Hessian in long format (see <https://cppad.readthedocs.io/latest/Hessian.html>).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{1}-dimensional space, then
//' the first \eqn{n} elements of the vector is the gradient of the partial derivative with respect to the first dimension of the function's domain.
//' The next \eqn{n} elements of the vector is the gradient of the partial derivative of the second dimension of the function's domain.
//' The Hessian as a matrix, can be obtained by using [`as.matrix()`] with `ncol = n`.
//' @return A `Rcpp::XPtr` to a `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun  tapeHessian(pADFun & pfun);


//' @title Tape the Gradient Offset of a Quadratic CppAD Tape
//' @family tape builders
//' @inheritParams pTapeJacobian
//' @return A `Rcpp::XPtr` to a `Rcpp_ADFun` object. The independent argument to the function are the dynamic parameters of `pfun`.
//' @details A quadratic function can be written as
//' \deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
//' The gradient of \eqn{f(x; \theta)} with respect to \eqn{x} is
//' \deqn{\Delta f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T)x + b(\theta).}
//' The Hessian is 
//' \deqn{H f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T),}
//' which does not depend on \eqn{x},
//' so the gradient of the function can be rewritten as
//' \deqn{\Delta f(x;\theta) = H f(x; \theta) x + b(\theta)^T.}
//' The tape calculates \eqn{b(\theta)} as
//'  \deqn{b(\theta) = \Delta f(x;\theta) - H f(x; \theta) x,}
//' which does not depend on \eqn{x}.
//' @export
// [[Rcpp::export]]
pADFun tapeGradOffset(pADFun & pfun);


//' @title Tape the log of Jacobian determinant of a CppAD Tape
//' @family tape builders
//' @param pfun An `Rcpp_ADFun` object.
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters (vector of size `0` if ther are no dynamic parameters).
//' @description Creates a tape of the log of the Jacobian determinant of a function taped by CppAD.
//' The `x` vector is used as the value to conduct the taping.
//' @return An `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun  tapeLogJacDet(pADFun & pfun);

//' @title Switch Dynamic and Independent Values of a Tape
//' @family tape builders
//' @description Convert an `Rcpp_ADFun` object so that the independent values become dynamic parameters
//' and the dynamic parameters become independent values
//' @param pfun An `Rcpp_ADFun` object. 
//' @param newvalue The independent value (in the sense after the switch has occurred) at which to tape the `Rcpp_ADFun`
//' @param newdynparam The value of the dynamic parameters (after the switch) at which to tape the `Rcpp_ADFun`
//' @return An `Rcpp_ADFun` object.
//' @export
// [[Rcpp::export]]
pADFun tapeSwap(pADFun & pfun);

# endif
