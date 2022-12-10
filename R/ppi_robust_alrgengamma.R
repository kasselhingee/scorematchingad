#' @title Windham Robustness for Scealy et al 2023
#' @description
#' Performs the Windham robustification algorithm exactly as described in \insertCite{scealy2023ro;textual}{scorecompdir} for score matching via log-ratio transform of the PPI model with \eqn{b_L = 0}. This method gives the same results as the more general implementation in [`Windham()`].
#' @inheritParams ppi_robust
#' @inheritSection ppi_robust return
#' @param ... Passed to a special version of [`Windham()`] and on to [`ppi()`]. The argument `fpcontrol` is not allowed as this is hardcoded into `ppi_robust_alrgengamma()`.
#' @details
#' This method must fit a PPI model via additive-log ratio transform with \eqn{b_L=0} fixed and the final element of \eqn{\beta} fixed.
#' @references
#' \insertAllCited{}
#' @export
ppi_robust_alrgengamma <- function(Y, cW, ...){
  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    return(drop(dppi(Y, paramvec = theta)))
  }

  JSConvergenceMetricThreshold = 1E-6 #used by Janice Scealy in original code
  # The following convergence metric matches the metric used by Janice Scealy in her original code: a threshold on the change in the first element of the beta estimate
  JSConvergenceMetric <- function(residuals){
    p <- (-1 + sqrt(1 + 8 + 8*length(residuals)))/2 #from quadratic formula: qstar = p(p-1)/2 + p-1
    return(abs(residuals[p*(p-1)/2 + 1])) #the first element after the AL is the first element of beta - this is to math
  }

  est = Windham_raw(Y = Y,
                    estimator = ppi,
                    ldenfun = ldenfun,
                    cW = cW,
                    trans = "alr",
                    ...,
                    fpcontrol = list(
                      ConvergenceMetric = JSConvergenceMetric,
                      ConvergenceMetricThreshold = JSConvergenceMetricThreshold),
                    multiplicativecorrection = FALSE)

  #make results nicer and consistent with ppi()
  out <- list(
    est = c(list(paramvec = est$theta), fromPPIparamvec(est$theta)),
    SE = "Not calculated.",
    info = est$optim
  )

  return(out)
}
