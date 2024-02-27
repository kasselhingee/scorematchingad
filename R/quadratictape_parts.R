#' @title Evaluate the Hessian and gradient offset of a Taped Quadratic Function
#' @family tape evaluators
#' @description
#' When the score matching discrepancy function is quadratic then the gradient of the score matching discrepancy function can be written using the Hessian and an offset term. This can be useful for solving for the situation when the gradient is zero.
#' The Hessian and offset term are computed using `CppAD` tapes.
#' Taylor approximation can be used for locations at removable singularities.
#' @details
#' A quadratic function can be written
#' \deqn{f(x; t) = \frac{1}{2} x^T W(t) x + b(t)^T x + c,}
#' where \eqn{t} is considered a vector that is constant with respect to the differentiation.
#' The Hessian of the function is with respect to \eqn{x} is
#' \deqn{H f(x; t) = \frac{1}{2}(W(t) + W(t)^T).}
#' The gradient of the function with respect to \eqn{x} can then be written
#' \deqn{\Delta f(x;t) = H f(x; t) x + b(t)^T x,}
#' where the Hessian and offset \eqn{b(t)} depend only on \eqn{t}.
#'
#' The functions here evaluate the Hessian and offset \eqn{b(t)} for many values of \eqn{t}.
#' Tapes of the Hessian and gradient offset are created using [`tapeHessian()`] and [`tapeGradOffset()`] respectively.
#' These tapes are then evaluated for every row of `tmat`.
#' When the `tcentres` row is not `NA`, then approximate results are calculated using [`pTaylorApprox()`]
#' @param tape A tape of a quadratic function (such as score matching discrepancy function)
#' @param tmat A matrix of `t` vectors. Each row corresponds to a vector.
#' @param tcentres A matrix of Taylor approximation centres for rows of `tmat` that require approximation. `NA` for rows that do not require approximation.
#' @param approxorder The order of the Taylor approximation to use.
#' @return A list of
#'  + `offset` Array of offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#'  + `Hessian` Array of offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#' @export
quadratictape_parts <- function(tape, tmat, tcentres = NA * tmat, approxorder = 10){
  stopifnot(inherits(tape, "ADFun"))
  stopifnot(nrow(tmat) == nrow(tcentres))
  stopifnot(testquadratictape(tape))
  toapprox <- !is.na(tcentres[, 1])

  Hesstape <- tapeHessian(tape)
  OffsetTape <- tapeGradOffset(tape)
  Hesstape_switched <- tapeSwap(Hesstape) #Hesstape is wrt to x (which in smo world is actually the model parameter set), but we want it to be wrt to the dynamic parameter like OffsetTape is

  #exact and approximat evalution of Hess
  fakeparametermat <- matrix(0, 
                            ncol = length(tape$xtape),
                            nrow = nrow(tmat))   #fake because the results don't depend on it
  emptyparametermat <- matrix(vector(mode = "double"), ncol = 0, nrow = nrow(tmat)) 
  # for OffsetTape which has no dynamic parameters

  Hesss <- evaltape(Hesstape_switched, xmat = tmat, pmat = fakeparametermat, xcentres = tcentres, approxorder = approxorder)
  offsets <- evaltape(OffsetTape, xmat = tmat, pmat = emptyparametermat, xcentres = tcentres, approxorder = approxorder)

  return(list(
    offset = offsets,
    Hessian = Hesss
  ))
}

