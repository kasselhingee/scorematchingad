# ifndef CPPAD_INTERFACE
# define CPPAD_INTERFACE

//for content that is Rcpp specific
#include "scorecompdir_types.h"
#include "mycpp/wrapas.hpp"  //needed because converting veca1 from R

//' @title Switch Dynamic and Independent Values of a Tape
//' @description Convert an ADFun so that the independent values become dynamic parameters
//' and the dynamic parameters become independent values
//' @param newvalue The independent value (in the sense after the switch has occurred) at which to tape the ADFun
//' @param newdynparam The value of the dynamic parameters (after the switch) at which to tape the ADFun
//' @return A pointer to an ADFun
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > swapDynamic(Rcpp::XPtr< CppAD::ADFun<double> > pfun, veca1 newvalue, veca1 newdynparam);

//' @title Evaluate the Jacobian of a tape
//' @param pfun Rcpp::XPtr to an ADFun
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters. If `pfun` has no dynamic parameters then set `dynparam = vector(mode = "numeric")`.
//' @return The Jacobian of pfun
//' @export
// [[Rcpp::export]]
vecd pJacobian(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta);

//' @title Evaluate a CppAD tape
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters.
//' @return The value of `pfun` evaluated at `x` with parameters `dynparam`.
//' @export
// [[Rcpp::export]]
vecd pForward0(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam);

//' @noRd
//' @title OBSOLETE: The Hessian of recorded function. Used only in smobj.R
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The Hessian of pfun
// [[Rcpp::export]]
vecd pHessian(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta);

//' @title Tape the Jacobian of CppAD Tape
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters
//' @description Creates a tape of the Jacobian of function taped by CppAD.
//' When the function returns a real value (as is the case for densities and the score matching objective) the Jacobian is equivalent to the gradient.
//' The `x` vector is used as the value to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Jacobian in long format (see [https://coin-or.github.io/CppAD/doc/jacobian.htm]).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{m}-dimensional space, then
//' the first \eqn{n} elements of vector is the gradient of the first component of function output.
//' The next \eqn{n} elements of the vector is the gradient of the second component of the function output.
//' The Jacobian as a matrix, could then be obtained by [`as.matrix()`] with `byrow = TRUE` and `ncol = n`.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  pTapeJacobian(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);

//' @title Tape the Hessian of a CppAD Tape
//' @inheritParams pTapeJacobian
//' @description Creates a tape of the Hessian of a function taped by CppAD.
//' The taped function represented by `pfun` must be scalar-valued (i.e. a vector of length 1).
//' The `x` vector and `dynparam` are used as the values to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Hessian in long format (see [https://coin-or.github.io/CppAD/doc/hessian.htm]).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{1}-dimensional space, then
//' the first \eqn{n} elements of the vector is the gradient of the partial derivative with respect to the first dimension of the function's domain.
//' The next \eqn{n} elements of the vector is the gradient of the partial derivative of the second dimension of the function's domain.
//' The Hessian as a matrix, can be obtained by using [`as.matrix()`] with `ncol = n`.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  pTapeHessian(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);

//' @title Indicate Constant Components of Range
//' @description Use `CppAD`'s `Parameter()` function for `ADFun` objects to see if the returned values of a tape are constant with respect to the independent values.
//' @param pfun A CppAD tape.
//' @param dynparam A set of dynamic parameters for `pfun`.
//' @return A vector logical values. `TRUE` indicates that element of the tape result is constant.
//' @details The `CppAD` function `Parameter(i)` [https://coin-or.github.io/CppAD/doc/fun_property.htm] returns `TRUE` when the `i`th component of the range does not depend on the independent value
//' (the `i`th component may still depend on the value of the dynamic parameters (see 'Dynamic' in [https://coin-or.github.io/CppAD/doc/glossary.htm#Parameter]) ).
//' @export
// [[Rcpp::export]]
std::vector<bool> pParameter(Rcpp::XPtr< CppAD::ADFun<double> > pfun);
// According to the help, applying Variable(u) to each return value would be false if u depends on the dynamic parameters and does not depend on the independent variable vector.

//' @title Tape the Gradient Offset of a Quadratic CppAD Tape
//' @inheritParams pTapeJacobian
//' @description A quadratic function can be written as
//' \deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
//' The function `pTapeGradOffset` creates a tape of \eqn{b(\theta)} where \eqn{\theta} is the independent variable.
//' @param pfun A quadratic CppAD Tape. Test for quadratic form using [`testquadratictape()`].
//' @details
//' The gradient of \eqn{f(x; \theta)} with respect to \eqn{\theta} is
//' \deqn{\Delta f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T)\theta + b(\theta),}
//' and the Hessian is 
//' \deqn{H f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T),}
//' which does not depend on \eqn{x}.
//' The gradient of the function can be rewritten as
//' \deqn{\Delta f(x;\theta) = H f(x; \theta) x + b(\theta)^T x.}
//' The tape calculates \eqn{b(\theta)} as
//'  \deqn{b(\theta) = \Delta f(x;\theta) - H f(x; \theta) x},
//' which does not depend on \eqn{x}.
//' In `pTapeGradOffset()` the `x` provided as an argument is used as the template for calculating \eqn{b(\theta)}.
//' The `x` vector and `dynparam` are used as the values to conduct the taping.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object. The independent argument to the function are the dynamic parameters of `pfun`.
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  pTapeGradOffset(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);


# endif
