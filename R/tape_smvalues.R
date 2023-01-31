#' @title Compute score matching objective value, gradient, and Hessian
#' @description Computes a range of relevant information for investigating score matching estimators.
#' @inheritParams tape_eval
#' @param smotape A taped score matching objective. Most easily created by [`buildsmotape()`].
#' @details The score matching objective values are differ from the Hyvarinen Divergence by a constant, see ... .
#' The gradient and Hessian are returned as arrays of row-vectors with each row corresponding to a row in `xmat` and `pmat`. 
#' @examples
#' m <- ppi_egmodel(100)
#' tapes <- buildsmotape("sphere", "ppi",
#'               utape = rep(1/m$p, m$p),
#'               usertheta = ppi_paramvec(beta = m$beta),
#'               weightname = "minsq", acut = 0.01)
#' tape_smvalues(tapes$smotape, xmat = m$sample, pmat = m$theta[1:5])
#' @export
tape_smvalues <- function(smotape, xmat, pmat, xcentres = NA * xmat, approxorder = 10){
  # prepare tapes
  Jsmofun <- pTapeJacobian(smotape, attr(smotape, "xtape"), attr(smotape, "dyntape"))
  Hsmofun <- pTapeJacobian(Jsmofun, attr(smotape, "xtape"), attr(smotape, "dyntape"))
  
  smofun_u <- swapDynamic(smotape, attr(smotape, "dyntape"), attr(smotape, "xtape")) #don't use a boundary point here!
  Jsmofun_u <- swapDynamic(Jsmofun, attr(smotape, "dyntape"), attr(smotape, "xtape"))
  Hsmofun_u <- swapDynamic(Hsmofun, attr(smotape, "dyntape"), attr(smotape, "xtape"))


  smovals <- tape_eval(smofun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)

  gradvals <- tape_eval(Jsmofun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)
  hessvals <- tape_eval(Hsmofun_u, xmat, pmat, xcentres = xcentres, approxorder = approxorder)

  return(list(
    obj = smovals,
    grad = gradvals,
    hess = hessvals
  ))
}
#' @return
#' A list of 
#'  + `obj` the score matching objective values
#'  + `grad` the gradient of the score matching objective
#'  + `hess` the Hessian of the score matching objective

#' @rdname tape_smvalues
#' @param w Weights to apply to each row of `xmat` for computing the weighted sum. If `NULL` then each row is given a weight of `1`.
tape_smvalues_wsum <- function(tape, xmat, pmat, w=NULL, xcentres = NA * xmat, approxorder = 10){
  evals_l <- tape_smvalues(tape, xmat = xmat, pmat = pmat,
                     xcentres = xcentres, approxorder = approxorder)
  
  # do weight checks afterwards so that eval results can be used to choose weights
  out <- lapply(evals_l, wcolSums, w = w)
  return(out)
}

