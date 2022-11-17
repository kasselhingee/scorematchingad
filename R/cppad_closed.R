#' @title The value when the gradient is zero given Hessian and offset
#' @param 

cppad_closed <- function(smotape, Y, w = rep(1, nrow(Y))){
  stopifnot(nrow(Y) == length(w))
  parts <- quadratictape_parts(smotape, Y)
  # weight parts
  parts$offset <- parts$offset * w
  parts$Hessian <- parts$Hessian * w
  offset <- colSums(parts$offset)
  Hess <- matrix(colSums(parts$Hessian), ncol = sqrt(ncol(parts$Hessian)))
  invHess <- solve(Hess)
  root <- drop(-1 * invHess %*% offset)
  return(list(
    est = root,
    Hessian = Hess,
    offset = offset
  ))
}


