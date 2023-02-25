#' @name moretapebuilders
#' @title Build New `CppAD` Tapes from Existing Tapes
#' @description Build new `CppAD` tapes from existing tapes, including differentiation, swapping independent and dynamic variables, and Jacobian determinants.
#' @family tape builders
#' @param tape An [`ADFun`] object.
#' @return An [`ADFun`] object.
#' @details The information in the fields `xtape` and `dyntape` of `tape` are used for the taping.
#' @seealso [`ADFun`]
NULL

#' @describeIn moretapebuilders Tape the Jacobian of a tape. The resulting tape returns the Jacobian as a vector.
#' @details
#' ## tapeJacobian 
#' The returned vector is ordered with the range elements iterating fastest, then the domain elements. See <https://cppad.readthedocs.io/en/latest/Jacobian.html>.
#' @export
tapeJacobian <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- pTapeJacobian(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("d", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

#' @describeIn moretapebuilders Tape the Hessian of a tape. The resulting tape returns the Jacobian as a vector (see <https://cppad.readthedocs.io/en/latest/Hessian.html>).
#' @details
#' ## tapeHessian
#' Suppose the function represented by `tape` maps from \eqn{n}-dimensional space to \eqn{1}-dimensional space, then
#' the first \eqn{n} elements of the vector is the gradient of the partial derivative with respect to the first dimension of the function's domain.
#' The next \eqn{n} elements of the vector is the gradient of the partial derivative of the second dimension of the function's domain.
#' The Hessian as a matrix, can be obtained by using [`as.matrix()`] with `ncol = n`.
#' @export
tapeHessian <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- pTapeHessian(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("d^2", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

#' @describeIn moretapebuilders A quadratic function can be written as
#' \deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
#' The function `tapeGradOffset` creates a tape of \eqn{b(\theta)} where \eqn{\theta} is the independent variable.
#' @details
#' ## tapeGradOffset
#' A quadratic function can be written as
#' \deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
#' The gradient of \eqn{f(x; \theta)} with respect to \eqn{x} is
#' \deqn{\Delta f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T)x + b(\theta).}
#' The Hessian is 
#' \deqn{H f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T),}
#' which does not depend on \eqn{x},
#' so the gradient of the function can be rewritten as
#' \deqn{\Delta f(x;\theta) = H f(x; \theta) x + b(\theta)^T.}
#' The tape calculates \eqn{b(\theta)} as
#'  \deqn{b(\theta) = \Delta f(x;\theta) - H f(x; \theta) x,}
#' which does not depend on \eqn{x}.
#' @export
tapeGradOffset <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- pTapeGradOffset(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("doffset:", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

#' @describeIn moretapebuilders
#' Creates a tape of the log of the Jacobian determinant of a function taped in `tape`.
#' @export
tapeLogJacDet <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- ptapelogdetJ(tape$ptr, tape$xtape, tape$dyntape)
  ADFun$new(outptr, 
            name = paste0("logJdet:", tape$name), 
            xtape = tape$xtape, 
            dyntape = tape$dyntape, 
            usertheta = tape$usertheta)
}

#' @describeIn moretapebuilders
#' Convert an ADFun so that the independent values become dynamic parameters
#' and the dynamic parameters become independent values.
#' @export
tapeSwap <- function(tape){
  stopifnot(inherits(tape, "ADFun"))
  outptr <- swapDynamic(tape$ptr, tape$dyntape, tape$xtape)
  ADFun$new(outptr, 
            name = paste0("d", tape$name), 
            xtape = tape$dyntape, 
            dyntape = tape$xtape, 
            usertheta = rep(NA_real_, length(tape$xtape)))
}

