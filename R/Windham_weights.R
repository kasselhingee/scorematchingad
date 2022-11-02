#' @title Windham weights for a given parameter vector
#' @description Computes Windham's robustness weights [ref Windham 1995] for a density proportional to 
#' \eqn{\exp(\eta(\theta) \cdot T(x))} where \eqn{\eta(\theta)} is linear and \eqn{x} is an observation (potentially multivariate).
#' @param cW The robustness tuning constants. One value per element of theta. Easily created for the PPI model using [ppi_cw()] and [ppi_cW_auto()].
#' @param ldenfun A (possibly improper) log density function taking two arguments, `Y` and `theta`.
#' @param theta Parameters for the model
#' @param Y A matrix of measurements. Each row a measurement.
#' @details
#' The elements of theta will be multiplied by cW for calculating the weights.
#' @return
#' For each row of Y, the density at parameter set `cW * theta` is computed using the exponential of `ldenfun`. The returned values are normalised to sum to 1.
#' @family Windham robust functions
#' @export
Windham_weights <- function(ldenfun, Y, theta, cW){
  if (is.null(ldenfun)){stop("ldenfun is NULL")}
  stopifnot(length(cW) == length(theta))
  stopifnot(is.numeric(cW))
  thetaforweights <- cW * theta #the elements of theta with FALSE inWW will be set to zero
  weights <- exp(ldenfun(Y = Y, theta = thetaforweights))
  weights=weights/sum(weights)
  return(weights)
}

