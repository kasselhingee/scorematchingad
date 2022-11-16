#' @title Test Whether a CppAD Tape is a Quadratic Function
#' @description
#' Uses [`pParameter()`] and derivatives [`pTapeJacobian()`], [`pJacobian()`] to test whether
#' the tape is quadratic.
#' Also returns some information about the Hessian.
#' @param xmat If passed, the third-order derivatives at values of the rows of `xmat` are tested.
#' @param dynparammat The dynamic parameters for the tape. If passed, the rows of `dynparammat` are passed to the tape as `dynparam`.
#' @details Uses the `xtape` and `dyntape` attributes of `tape` to create new tapes.
#' A tape of Hessian is obtained by applying [`pTapeJacobian()`] twice. Using [`pTapeHessian()`] directly did not show constant parameters via [`pParameter()`] in tests.
#' [`pParameter()`] is evaluated with dynamic parameters given by the `dyntape` attribute of `tape`.
testquadratictape <- function(tape, xmat = NULL, dynparammat = NULL){
  tapeJ <- pTapeJacobian(tape, attr(tape, "xtape"), attr(tape, "dyntape"))
  tapeH <- pTapeJacobian(tapeJ, attr(tape, "xtape"), attr(tape, "dyntape"))

  #pParameter() test
  isparameter <- pParameter(tapeH, attr(tape, "dyntape"))
  result_pParameter <- all(isparameter)

  #try Jacobian
  result_thirdderiv <- NA
  if (!is.null(xmat) || !is.null(dynparammat)){
    stopifnot(nrow(xmat) == nrow(dynparammat))
    thirdderivs <- lapply(1:nrow(xmat), function(i){
      pJacobian(tapeH, xmat[i, ], dynparammat[i, ])
    })
    isallzero <- unlist(lapply(thirdderivs, function(vec){all(vec == 0)}))
    result_thirdderiv <- all(isallzero)
  }

  #properties of Hessian. If Hessian is constant it should matter the values used
  propsHessian <- list()
  if (isTRUE(result_pParameter) | isTRUE(result_thirdderiv)){
    if (is.null(xmat)){
      hess <- pForward0(tapeH, attr(tape, "xtape"), attr(tape, "dyntape"))
    } else {
      hess <- pForward0(tapeH, xmat[1, ], dynparammat[1, ])
    }
    hessmat <- matrix(hess, byrow = TRUE, ncol = sqrt(length(hess)))
    propsHessian$symmetric <- isSymmetric(hessmat)
    propsHessian$determinant <- determinant(hessmat, logarithm = FALSE)
    propsHessian$evalues <- eigen(hessmat)$values
    propsHessian$posdefinite <- all(propsHessian$evalues > 0) #a symmetric matrix is positive definite iff all its eigen values are positive
  } 

}
