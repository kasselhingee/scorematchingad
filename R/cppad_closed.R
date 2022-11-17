#' @title The value when the gradient is zero given Hessian and offset
#' @param Yboundarycentres A matrix of Taylor approximation centres for rows of Y that require approximation. `NA` for rows that do not require approximation.

cppad_closed <- function(smotape, Y, Yapproxcentres = NA * Y, 
                         w = rep(1, nrow(Y)),
                         approxorder = 10){
  stopifnot(nrow(Y) == length(w))
  stopifnot(nrow(Y) == nrow(Yapproxcentres))
  
  parts <- quadratictape_parts(smotape, tmat = Y,
                               tcentres = Yapproxcentres,
                               approxorder = approxorder)

  # weight parts
  parts$offset <- parts$offset * w
  parts$Hessian <- parts$Hessian * w
  offset <- colSums(parts$offset)
  Hess <- colSums(parts$Hessian)
  Hess <- matrix(Hess, ncol = sqrt(ncol(parts$Hessian)))
  invHess <- solve(Hess)
  root <- drop(-1 * invHess %*% offset)
  return(list(
    est = root,
    Hessian = Hess,
    offset = offset
  ))
}


