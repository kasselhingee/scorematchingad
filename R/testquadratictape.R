#' @title Test Whether a CppAD Tape is a Quadratic Function
#' @description
#' Uses [`pParameter()`] and derivatives [`pTapeJacobian()`], [`pJacobian()`] to test whether
#' the tape is quadratic.
#' @param tape An `Rcpp::XPtr` to a CppAD tape.
#' @param xmat If passed, the third-order derivatives at values of the rows of `xmat` are tested.
#' @param dynparammat The dynamic parameters for the tape. If passed, the rows of `dynparammat` are passed to the tape as `dynparam`.
#' @param verbose If TRUE information about the failed tests is passed.
#' @details Uses the `xtape` and `dyntape` attributes of `tape` to create new tapes.
#' A tape of Hessian is obtained by applying [`pTapeJacobian()`] twice. Using [`pTapeHessian()`] directly did not show constant parameters via [`pParameter()`] in tests.
#'
#' Two tests are conducted on the tape of the Hessian.
#' For a function of quadratic form, [`pParameter()`] should return a vector of `TRUE` values.
#' The other test evaluates the Jacobian at each row of `xmat` and `dynparammat` (if `xmat` and `dynparammat` are non-NULL). The result should be zero when `tape` represents a function of quadratic form.
#' If the results of the tests differ the returned value is `FALSE` and a message is printed indicating which test failed.
#' @return `TRUE` or `FALSE`
#' @examples
#'  sqrtman <- pmanifold("sphere")
#'  ppitape <- tapell(llname = "ppi",
#'                    xtape = c(0.2, 0.3, 0.5),
#'                    usertheta = ppi_paramvec(p = 3), 
#'                    pmanifoldtransform = sqrtman)
#'  ppismotape <- tapesmo(lltape = ppitape,
#'                        pmanifoldtransform = sqrtman,
#'                        divweight = "minsq",
#'                        acut = 0.1,
#'                        verbose = FALSE)
#'
#'  testquadratictape(ppismotape)
#' @export
testquadratictape <- function(tape, xmat = NULL, dynparammat = NULL, verbose = FALSE){
  tapeJ <- pTapeJacobian(tape, attr(tape, "xtape"), attr(tape, "dyntape"))
  tapeH <- pTapeJacobian(tapeJ, attr(tape, "xtape"), attr(tape, "dyntape"))

  #pParameter() test
  isparameter <- pParameter(tapeH)
  result_pParameter <- all(isparameter)
  if (verbose && !result_pParameter){
    message(sprintf("The Hessian was non-constant according to pParameter() for elements %s.",
                    paste(which(!isparameter), collapse = ", ")))
  }

  if (is.null(xmat) && is.null(dynparammat)){return(result_pParameter)}

  #try Jacobian
  stopifnot(isTRUE(nrow(xmat) == nrow(dynparammat)))
  thirdderivs <- lapply(1:nrow(xmat), function(i){
    pJacobian(tapeH, xmat[i, ], dynparammat[i, ])
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


