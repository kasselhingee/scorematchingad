#' @title Evaluate the Hessian and gradient offset of a Taped Quadratic Function
#' @description
#' When the score matching objective function is quadratic then the gradient of the score matching objective function can be written using the Hessian and an offset term. This can be useful for solving for the situation when the gradient is zero.
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
#' Tapes of the Hessian and gradient offset are created using [`pTapeHessian()`] and [`pTapeGradOffset()`] respectively.
#' These tapes are then evaluated for every row of `tmat`.
#' When the `tcentres` row is not `NA`, then approximate results are calculated using [`pTaylorApprox()`]
#' @param tape A tape of a quadratic function (such as score matching objective function)
#' @param tmat A matrix of `t` vectors. Each row corresponds to a vector.
#' @param tcentres A matrix of Taylor approximation centres for rows of `tmat` that require approximation. `NA` for rows that do not require approximation.
#' @param approxorder The order of the Taylor approximation to use.
#' @return A list of
#'  + `offset` Array of offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#'  + `Hessian` Array of offsets \eqn{b(t)}, each row corresponding to a row in `tmat`
#' @export
quadratictape_parts <- function(tape, tmat, tcentres = NA * tmat, approxorder = 10){
  stopifnot(nrow(tmat) == nrow(centres))
  stopifnot(testquadratictape(tape))
  toapprox <- !is.na(Yapproxcentres[, 1])

  Hesstape <- pTapeHessian(tape, attr(tape, "xtape"), attr(tape, "dyntape"))
  OffsetTape <- pTapeGradOffset(tape, attr(tape, "xtape"), attr(tape, "dyntape"))
  
  # prepare matrices to be filled
  alloffsets <- matrix(NA, ncol = length(attr(tape, "xtape")), nrow = nrow(tmat))
  allHesss <- matrix(NA, ncol = length(attr(tape, "xtape"))^2, nrow = nrow(tmat))

  # exact evaluations
  if (any(!toapprox)){
    #evaluate Hessians
    Hesss <- apply(tmat[!toapprox, , drop = FALSE], MARGIN = 1, 
       function(x){
         pForward0(Hesstape, 0*attr(tape, "xtape"), x)},
       simplify = FALSE)
    Hesss <- do.call(rbind, Hesss)

    #evaluate offsets
    offsets <- apply(tmat[!toapprox, , drop=FALSE], MARGIN = 1,
        function(x){
          pForward0(OffsetTape, x, vector(mode = "double"))}, 
        simplify = FALSE)
    offsets <- do.call(rbind, offsets)

    alloffsets[!toapprox, ] <- offsets
    allHesss[!toapprox, ] <- Hesss
  } 

  # approximate evaluations
  if (any(toapprox)){
    Hesstape_switched <- swapDynamic(Hesstape, attr(tape, "dyntape"), attr(tape, "xtape")) #Hesstape is wrt to x, but we want it to be wrt to the dynamic parameter like OffsetTape is
    #approximate Hessians
    Hesss <- lapply(1:nrow(tmat), function(i){
      pTaylorApprox(Hesstape_switched, tmat[i, ], centres[i, ], 0*attr(tape, "xtape"), approxorder)
    })
    Hesss <- do.call(rbind, Hesss)

    #approximate Offsets 
    offsets <- lapply(1:nrow(tmat), function(i){
      pTaylorApprox(OffsetTape, tmat[i, ], centres[i, ], 0*attr(tape, "xtape"), approxorder)
    })
    offsets <- do.call(rbind, offsets)
    
    alloffsets[toapprox, ] <- offsets
    allHesss[toapprox, ] <- Hesss
  } 


  return(list(
    offset = alloffsets,
    Hessian = allHesss
  ))
}

