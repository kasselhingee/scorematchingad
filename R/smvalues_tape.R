#' @title Compute score matching discrepancy value, gradient, and Hessian
#' @family tape evaluators
#' @description Computes a range of relevant information for investigating score matching estimators.
#' @inheritParams evaltape
#' @param smotape A taped score matching discrepancy. Most easily created by [`buildsmotape()`].
#' @details The score matching discrepancy values are differ from the Hyv\"arinen Divergence by a constant, see ... .
#' The gradient and Hessian are returned as arrays of row-vectors with each row corresponding to a row in `xmat` and `pmat`. 
#' @examples
#' m <- rppi_egmodel(100)
#' tapes <- buildsmotape("sim", "sqrt", "sph", "ppi",
#'               ytape = rep(1/m$p, m$p),
#'               usertheta = ppi_paramvec(beta = m$beta),
#'               bdryw = "minsq", acut = 0.01)
#' smvalues_tape(tapes$smotape, xmat = m$sample, pmat = m$theta[1:5])
#' @export
smvalues_tape <- function(smotape, xmat, pmat, xcentres = NA * xmat, approxorder = 10){
  stopifnot(inherits(smotape, "ADFun"))
  # prepare tapes
  Jsmofun <- tapeJacobian(smotape)
  Hsmofun <- tapeJacobian(Jsmofun)
  
  smofun_u <- tapeSwap(smotape) #don't use a boundary for taping!
  Jsmofun_u <- tapeSwap(Jsmofun)
  Hsmofun_u <- tapeSwap(Hsmofun)

  smovals <- evaltape(smofun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)

  gradvals <- evaltape(Jsmofun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)
  hessvals <- evaltape(Hsmofun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)

  return(list(
    obj = smovals,
    grad = gradvals,
    hess = hessvals
  ))
}
#' @return
#' A list of 
#'  + `obj` the score matching discrepancy values
#'  + `grad` the gradient of the score matching discrepancy
#'  + `hess` the Hessian of the score matching discrepancy

#' @rdname smvalues_tape
#' @param w Weights to apply to each row of `xmat` for computing the weighted sum. If `NULL` then each row is given a weight of `1`.
smvalues_tape_wsum <- function(tape, xmat, pmat, w=NULL, xcentres = NA * xmat, approxorder = 10){
  evals_l <- smvalues_tape(tape, xmat = xmat, pmat = pmat,
                     xcentres = xcentres, approxorder = approxorder)
  
  # do weight checks afterwards so that eval results can be used to choose weights
  out <- lapply(evals_l, wcolSums, w = w)
  return(out)
}

