# @title Generate approximation centres for measurements on the boundary of the simplex
#' @param utable The compositional measurements
#' @param shiftsize The distance to away from `u` and the simplex boundary to create the approximation centre.
#' @return The approximation centres for the measurements
#' @examples
#' m <- ppi_egmodel(10)
#' simplex_boundaryshift(m$sample, shiftsize = 1E-5)
# @export
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

# @title Return approximation centres for measurements near the simplex boundary.
#' @param u Measurement
#' @param bdrythreshold Measurments closer than `bdrythreshold` to the edge of the simplex will be considered boundary measurements.
#' @param shiftsize Boundary measurements are shifted `shiftsize` towards the simplex centre
#' @examples
#' u <- matrix(runif(5), nrow = 1)
#' u <- u / sqrt(sum(u^2))
# @export
simplex_isboundary <- function(utabl, bdrythreshold = 1E-15){
  isbdrypt <- apply(utabl, MARGIN = 1, min) < bdrythreshold
  return(isbdrypt)
}


simplex_boundarysplit <- function(utabl, bdrythreshold = 1E-15, shiftsize = 1E-10){
  onbdry <- simplex_isboundary(utabl, bdrythreshold = bdrythreshold)
  acentres <- simplex_boundaryshift(utabl[onbdry, , drop = FALSE], shiftsize = shiftsize)
  return(list(
    interior = utabl[!onbdry, , drop = FALSE],
    uboundary = utabl[onbdry, , drop = FALSE],
    boundaryapprox = acentres))
}

