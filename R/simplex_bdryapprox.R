#' @title Generate approximation centres for measurements on the boundary of the simplex
#' @description 
#' Shift points by `shiftsize` towards the centre of the simplex.
#' @details Measurements on the boundary of the simplex can lead to degenerate log-likelihood values, which are avoided in score matching by the divergence weight function. However, `CppAD` still computes the log-likelihood and derivatives, and (to KLH's knowledge) is not able to evaluate the value of \eqn{0/0} or similar.
#' Good approximations of the score matching objective value can be obtained through Taylor approximations, which require a *centre*.
#' The function here generates this approximation centre.
#' @param utable The compositional measurements
#' @param shiftsize The distance to away from `u` and the simplex boundary to create the approximation centre.
#' @return The approximation centres for the measurements (in other words, the points shifted towards the simplex centre).
#' @family Simplex boundary approximation helpers
#' @examples
#' m <- ppi_egmodel(10)
#' simplex_boundaryshift(m$sample, shiftsize = 1E-5)
#' @export
simplex_boundaryshift <- function(utabl, shiftsize = 1E-4){
  p <- ncol(utabl)
  middleofsimplex <- rep(1/p, p)
  shiftdir <- middleofsimplex - utabl
  shiftdir <- shiftdir / sqrt(rowSums(shiftdir^2)) #make a unit vector
  approxcentre <- utabl + shiftsize * shiftdir
  approxcentre <- approxcentre / rowSums(approxcentre) # normalise to get back onto simplex
  return(approxcentre)
}

approxcentre <- function(utabl, shiftsize = 1E-4){
  warning("approxcentre() will become obsolete")
  simplex_boundaryshift(utabl, shiftsize = shiftsize)
}

#' @title Determine whether points are on the simplex boundary
#' @description Tests whether points are within `bdrythreshold` distance of the boundary, where distance is the size of minimum component. 
#' @param utabl Measurement matrix.
#' @param bdrythreshold Measurments closer than `bdrythreshold` to the edge of the simplex will be considered boundary measurements.
#' @param shiftsize Boundary measurements are shifted `shiftsize` towards the simplex centre
#' @family Simplex boundary approximation helpers
#' @return A vector of `TRUE` or `FALSE` values. Rows of `utabl` corresponding to `TRUE` are on the boundary of the simplex.
#' @export
simplex_isboundary <- function(utabl, bdrythreshold = 1E-15){
  isbdrypt <- apply(utabl, MARGIN = 1, min) < bdrythreshold
  return(isbdrypt)
}


#' @title Detect boundary observations and generate Taylor approximation centres
#' @description 
#' Split the matrix of observations `utabl` into boundary points and interior points using [`simplex_isboundary()`].
#' For the boundary points Taylor approximation centres are computed with [`simplex_boundaryshift()`].
#' A weight vector, if supplied, is split consistently with the above split.
#' @details Measurements on the boundary of the simplex can lead to degenerate log-likelihood values, which are avoided in score matching by the divergence weight function. However, `CppAD` still computes the log-likelihood and derivatives, and (to KLH's knowledge) is not able to evaluate the value of \eqn{0/0} or similar.
#' Good approximations of the score matching objective value can be obtained through Taylor approximations, which requires a centre.
#' The function converts the data matrix `utabl` into the necessary sub-data sets.
#' @param utable The compositional measurements
#' @param shiftsize The distance to away from `u` and the simplex boundary to create the approximation centre.
#' @return A list of
#'  * `interior` Interior points.
#'  * `uboundary` Points on the boundary.
#'  * `boundaryapprox` Corresponding approximation centres for `uboundary`.
#'  * `winterior` Weights from `w` corresponding to `interior` (`NULL` if not weights supplied).
#'  * `wboundary` Weights from `w` corresponding to `uboundary` (`NULL` if not weights supplied).
#' @family Simplex boundary approximation helpers
#' @examples
#' m <- ppi_egmodel(10)
#' simplex_boundaryshift(m$sample, shiftsize = 1E-5)
#' @export
simplex_boundarysplit <- function(utabl, bdrythreshold = 1E-15, shiftsize = 1E-10, w = NULL){
  onbdry <- simplex_isboundary(utabl, bdrythreshold = bdrythreshold)
  acentres <- simplex_boundaryshift(utabl[onbdry, , drop = FALSE], shiftsize = shiftsize)
  return(list(
    interior = utabl[!onbdry, , drop = FALSE],
    uboundary = utabl[onbdry, , drop = FALSE],
    boundaryapprox = acentres,
    winterior = w[!onbdry, drop = FALSE],  #when w NULL, then this is NULL too
    wboundary = w[onbdry, drop = FALSE]))
}

