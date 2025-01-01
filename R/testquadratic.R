#' @title Test Whether a CppAD Tape is a Quadratic Function
#' @family tape evaluators
#' @description
#' Uses the [`CppAD` parameter property](https://cppad.readthedocs.io/latest/fun_property.html#parameter) and derivatives (via [`tape_Jacobian()`]) to test whether
#' the tape is quadratic.
#' @param tape An [`Rcpp_ADFun`] object.
#' @title Test Whether a CppAD Tape is a Quadratic Function
#' @family tape evaluators
#' @description
#' Uses the [`CppAD` parameter property](https://cppad.readthedocs.io/latest/fun_property.html#parameter) and derivatives (via [`tape_Jacobian()`]) to test whether
#' the tape is quadratic.
#' @param tape An `ADFun` object.
#' @param xmat The third-order derivatives at independent variable values of the rows of `xmat` and dynamic parameters from the rows of `dynmat` are tested.
#' @param dynmat  The third-order derivatives at independent variable values of the rows of `xmat` and dynamic parameters from the rows of `dynmat` are tested.
#' @param verbose If TRUE information about the failed tests is printed.
#' @details Uses the `xtape` and `dyntape` values stored in `tape` to create new tapes.
#' A tape of the Hessian is obtained by applying [`tape_Jacobian()`] twice, and then uses the [`CppAD` parameter property](https://cppad.readthedocs.io/latest/fun_property.html#parameter) to test whether the Hessian is constant. A function of quadratic form should have constant Hessian.
#'
#' If `xmat` and `dynmat` are non-`NULL` then `testquadratic()` also checks the Jacobian of the Hessian at `xmat` and `dynmat` values. For quadratic form functions the Jacobian of the Hessian should be zero.
#' @return `TRUE` or `FALSE`
#' @examples
#' tapes <- tape_smd(
#'    "sim", "sqrt", "sph",
#'    ll = "ppi",
#'    ytape = c(0.2, 0.3, 0.5),
#'    usertheta = ppi_paramvec(p = 3), 
#'    bdryw = "minsq",
#'    acut = 0.1,
#'    verbose = FALSE)
#'
#'  testquadratic(tapes$smdtape)
#' @export
testquadratic <- function(tape, xmat = matrix(tape$xtape, nrow = 1), dynmat = matrix(tape$dyntape, nrow = 1), verbose = FALSE){
  stopifnot(inherits(tape, "Rcpp_ADFun"))
  tapeJ <- tape_Jacobian(tape)
  tapeH <- tape_Jacobian(tapeJ)

  #parameter() test
  isparameter <- vapply(1:tapeH$range, function(i){tapeH$parameter(i-1)}, FUN.VALUE = TRUE)
  result_parameter <- all(isparameter)
  if (verbose && !result_parameter){
    message(sprintf("The Hessian was non-constant according to parameter() for elements %s.",
                    paste(which(!isparameter), collapse = ", ")))
  }

  #try Jacobian
  stopifnot(isTRUE(nrow(xmat) == nrow(dynmat)))
  thirdderivs <- lapply(1:nrow(xmat), function(i){
    tapeH$Jac(xmat[i, ], dynmat[i, ])
  })
  isallzero <- unlist(lapply(thirdderivs, function(vec){all(vec == 0)}))
  result_thirdderiv <- all(isallzero)
  if (verbose && !result_thirdderiv){
    message(sprintf("The Jacobian of the Hessian was non-zero for row %s of xmat and dynmat",
                    paste(which(!isallzero), collapse = ", ")))
  }

  finalresult <- result_thirdderiv && result_parameter
  return(finalresult)
}



