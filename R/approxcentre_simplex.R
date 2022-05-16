#' @title Generate approximation centres for measurements on the boundary of the simplex
#' @param utable The compositional measurements
#' @param shiftsize The distance to away from `u` and the simplex boundary to create the approximation centre.
#' @return The approximation centres for the measurements
#' @examples
#' m <- sec2_3model(10)
#' approxcentre(m$sample, shiftsize = 1E-5)
#' @export
approxcentre <- function(utabl, shiftsize = 1E-4){
  p <- ncol(utabl)
  middleofsimplex <- rep(1/p, p)
  shiftdir <- middleofsimplex - utabl
  shiftdir <- shiftdir / sqrt(rowSums(shiftdir^2)) #make a unit vector
  approxcentre <- utabl + shiftsize * shiftdir
  approxcentre <- approxcentre / rowSums(approxcentre) # normalise to get back onto simplex
  return(approxcentre)
}
