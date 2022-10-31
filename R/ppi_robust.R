#' @title Robustly Estimate Parameters of the PPI Distribution
#' @description Uses [WindhamRobust()] and [ppi()] to estimate a PPI distribution robustly.
#' @param Y A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.
#' @param cW A vector of robustness tuning constants. Easy to build using [ppi_cW()] and [ppi_cW_auto()]. See [WindhamRobust()] for more details on `cW`.
#' @param ... Passed to [ppi()] and [WindhamRobust()].
#' @details
#' There are many arguments to the [ppi()] function - I highly recommmend trialling your arguments on [ppi()] first before running `ppi_robust()`.
#' @export
ppi_robust <- function(Y, cW, ...){
  ellipsis::check_dots_used()
  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    return(drop(dppi(Y, paramvec = theta)))
  }

  est = WindhamRobust(Y = Y,
                    estimator = ppi,
                    ldenfun = ldenfun,
                    cW = cW,
                    ...)

  #make results nicer and consistent with ppi()
  out <- list(
    est = c(list(paramvec = est$theta), fromPPIparamvec(est$theta)),
    SE = "Not calculated.",
    info = est$optim
  )

  return(out)
}
