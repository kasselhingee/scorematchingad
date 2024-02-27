#' @title Test Whether a CppAD Tape is a Quadratic Function
#' @family tape evaluators
#' @description
#' Uses [`pParameter()`] and derivatives [`tapeJacobian()`], [`pJacobian()`] to test whether
#' the tape is quadratic.
#' @param tape An `ADFun` object.
#' @param xmat If passed, the third-order derivatives at values of the rows of `xmat` are tested.
#' @param dynparammat The dynamic parameters for the tape. If passed, the rows of `dynparammat` are passed to the tape as `dynparam`.
#' @param verbose If TRUE information about the failed tests is passed.
#' @details Uses the `xtape` and `dyntape` values stored in `tape` to create new tapes.
#' A tape of Hessian is obtained by applying [`tapeJacobian()`] twice. Using [`tapeHessian()`] directly did not show the correct constant parameters via [`pParameter()`] in tests.
#'
#' Two tests are conducted on the tape of the Hessian.
#' For a function of quadratic form, [`pParameter()`] should return a vector of `TRUE` values.
#' The other test evaluates the Jacobian at each row of `xmat` and `dynparammat` (if `xmat` and `dynparammat` are non-NULL). The result should be zero when `tape` represents a function of quadratic form.
#' If the results of the tests differ the returned value is `FALSE` and a message is printed indicating which test failed.
#' @return `TRUE` or `FALSE`
#' @examples
#' tapes <- buildsmotape(
#'    "sim", "sqrt", "sph",
#'    llname = "ppi",
#'    ytape = c(0.2, 0.3, 0.5),
#'    usertheta = ppi_paramvec(p = 3), 
#'    bdryw = "minsq",
#'    acut = 0.1,
#'    verbose = FALSE)
#'
#'  testquadratictape(tapes$smotape)
#' @export
testquadratictape <- function(tape, xmat = NULL, dynparammat = NULL, verbose = FALSE){
  stopifnot(inherits(tape, "ADFun"))
  tapeJ <- tapeJacobian(tape)
  tapeH <- tapeJacobian(tapeJ)

  #pParameter() test
  isparameter <- pParameter(tapeH$ptr)
  result_pParameter <- all(isparameter)
  if (verbose && !result_pParameter){
    message(sprintf("The Hessian was non-constant according to pParameter() for elements %s.",
                    paste(which(!isparameter), collapse = ", ")))
  }

  if (is.null(xmat) && is.null(dynparammat)){return(result_pParameter)}

  #try Jacobian
  stopifnot(isTRUE(nrow(xmat) == nrow(dynparammat)))
  thirdderivs <- lapply(1:nrow(xmat), function(i){
    pJacobian(tapeH$ptr, xmat[i, ], dynparammat[i, ])
  })
  isallzero <- unlist(lapply(thirdderivs, function(vec){all(vec == 0)}))
  result_thirdderiv <- all(isallzero)
  if (verbose && !result_thirdderiv){
    message(sprintf("The Jacobian of the Hessian was non-zero for row %s of xmat and dynparammat",
                    paste(which(!isallzero), collapse = ", ")))
  }

  finalresult <- result_thirdderiv && result_pParameter
  return(finalresult)
}


