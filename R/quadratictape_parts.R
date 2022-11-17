#' @title Evaluate the Hessian and gradient offset of a Taped Quadratic Function
#' @description
#' When the score matching objective function is quadratic then the gradient of the score matching objective function can be written using the Hessian and an offset term. This can be useful for solving for the situation when the gradient is zero.
#' The Hessian and offset term are computed using `CppAD` tapes.
#' The function `quadratictape_parts_approx()` evaluates the Hessian and offset using Taylor approximation.
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
#'
#' Tapes of the Hessian and gradient offset are created using [`pTapeHessian()`] and [`pTapeGradOffset()`] respectively.
#' These tapes are then evaluated for every row of `tmat`.
#' @param tape A tape of a quadratic function (such as score matching objective function)
#' @param tmat A matrix of `t` vectors. Each row corresponds to a vector.
#' @return A list of
#'  + `offset` Array of computed offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#'  + `Hessian` Array of computed offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#' @export
quadratictape_parts <- function(tape, tmat){

  stopifnot(testquadratictape(tape))

  Hesstape <- pTapeHessian(tape, attr(tape, "xtape"), attr(tape, "dyntape"))
  OffsetTape <- pTapeGradOffset(tape, attr(tape, "xtape"), attr(tape, "dyntape"))

  #evaluate Hessians
  Hesss <- apply(tmat, MARGIN = 1, function(x){pForward0(Hesstape, 0*attr(tape, "xtape"), x)}, simplify = FALSE)
  Hesss <- do.call(rbind, Hesss)
 
  #evaluate Offsets 
  offsets <- apply(tmat, MARGIN = 1, function(x){pForward0(OffsetTape, x, vector(mode = "double"))}, simplify = FALSE)
  offsets <- do.call(rbind, offsets)
  return(list(
    offset = offsets,
    Hessian = Hesss
  ))
}
