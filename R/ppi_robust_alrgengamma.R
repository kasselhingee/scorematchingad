#' @title Windham Robustness for Scealy et al 2023
#' @description
#' Performs the Windham robustification algorithm exactly as described in \insertCite{scealy2023ro;textual}{scorecompdir} for score matching via log-ratio transform of the PPI model with \eqn{b_L = 0}. This method gives the same results as the more general implementation in [`Windham()`].
#' @inheritParams ppi_robust
#' @inheritSection ppi_robust return
#' @details
#' This method must fit a PPI model via additive-log ratio transform with \eqn{b_L=0} fixed and the final element of \eqn{\beta} fixed.
#' @references
#' \insertAllCited{}
#' @export
ppi_robust_alrgengamma <- function(Y, cW, ...){
  ellipsis::check_dots_used()
  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    return(drop(dppi(Y, paramvec = theta)))
  }

  est = Windham_raw(Y = Y,
                    estimator = ppi,
                    ldenfun = ldenfun,
                    cW = cW,
                    trans = "alr",
                    ...,
                    multiplicativecorrection = FALSE)

  #make results nicer and consistent with ppi()
  out <- list(
    est = c(list(paramvec = est$theta), fromPPIparamvec(est$theta)),
    SE = "Not calculated.",
    info = est$optim
  )

  return(out)
}
