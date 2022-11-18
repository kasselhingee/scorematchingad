#' @title Simple wrapper to evaluate CppAD tapes and derivatives many times
#' @param tape A `CppAD` `tape`.
#' @param xmat A matrix of (multivariate) independent variables. Each represents a single independent variable vector.
#' @param pmat A matrix of dynamic parameters. Or a single vector of dynamic parameters to use for all rows of `xmat`.
#' @param xcentres A matrix of approximation for Taylor approximation centres for `xmat`. Use values of `NA` for rows that do not require Taylor approximation.
#' @description Evaluates a tape exactly or approximately for an array of provided variable values and dynamic parameter values.
#' The function `tape_eval_wsum()` computes the column-wise weighted sum of the result.
#' @details
#' Approximation is via Taylor approximation of the independent variable around the approximation centre provided in `xcentres`.
#' @return
#' A matrix, each row corresponding to the evaluation of the same row in `xmat`, `pmat` and the `xcentres`.
#' @export
tape_eval <- function(tape, xmat, pmat, xcentres = NA * xmat, approxorder = 10){
  stopifnot(nrow(xmat) == nrow(xcentres))
  if (is.vector(pmat) | isTRUE(nrow(pmat) == 1)){
    pmat <- matrix(pmat, byrow = TRUE, nrow = nrow(xmat), ncol = length(pmat))
  }
  else {stopifnot(nrow(xmat) == nrow(pmat))}
  toapprox <- !is.na(xcentres[, 1])

  evals_l <- list()
  # exact evaluations
  if (any(!toapprox)){
    evals_l[!toapprox] <- lapply(which(!toapprox), function(i){
      pForward0(tape, xmat[i, ], pmat[i, ])
    })
  }
  if (any(toapprox)){
    evals_l[toapprox] <- lapply(which(toapprox), function(i){
      pTaylorApprox(tape, xmat[i, ], xcentres[i, ], pmat[i, ], approxorder)
    })
  }

  evals <- do.call(rbind, evals_l)
  return(evals)
}

#' @rdname tape_eval
#' @param w Weights to apply to each row of `xmat` for computing the weighted sum.
tape_eval_wsum <- function(tape, xmat, pmat, w=rep(1, nrow(xmat)), xcentres = NA * xmat, approxorder = 10){
  stopifnot(length(w) == nrow(xmat))
  evals <- tape_eval(tape, xmat = xmat, pmat = pmat,
                     xcentres = xcentres, approxorder = approxorder)
  wevals <- evals*w
  return(colSums(wevals))
}


