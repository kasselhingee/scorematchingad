#' @title The value when the gradient is zero given Hessian and offset
#' @param Yboundarycentres A matrix of Taylor approximation centres for rows of Y that require approximation. `NA` for rows that do not require approximation.

cppad_closed <- function(smotape, Y, Yapproxcentres = NA * Y, 
                         w = rep(1, nrow(Y)),
                         approxorder = 10){
  stopifnot(nrow(Y) == length(w))
  stopifnot(nrow(Y) == nrow(Yapproxcentres))
  toapprox <- !is.na(Yapproxcentres[, 1])

  parts <- list(offset = matrix(NA, ncol = length(attr(smotape, "xtape")), nrow = nrow(Y)),
                Hessian = matrix(NA, ncol = length(attr(smotape, "xtape"))^2, nrow = nrow(Y)))
  if (any(!toapprox)){
    exactparts <- quadratictape_parts(smotape, Y[!toapprox, , drop = FALSE])
    parts$offset[!toapprox, ] <- exactparts$offset
    parts$Hessian[!toapprox, ] <- exactparts$Hessian
  } 
  # repeat process nearly identically for boundary points - just use approximations instead
  if (any(toapprox)){
    approxparts <- quadratictape_parts_approx(smotape,
                      Y[toapprox, , drop = FALSE],
                      centres = Yapproxcentres[toapprox, , drop = FALSE], 
                      order = approxorder)
    parts$offset[toapprox, ] <- approxparts$offset
    parts$Hessian[toapprox, ] <- approxparts$Hessian
  }

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


