#' @title Functions for converting between a user friendly theta and other forms of theta
#' @param usertheta A parameter vector. NA elements are to be fitted other elements are fixed.
#' @param isfixed A boolean vector same length as the parameter vector. `TRUE` values are fixed at the value of `starttheta`, `FALSE` are fitted.
#' @description Converts a `usertheta` to `isfixed`.
t_u2i <- function(usertheta){
  isfixed <- !is.na(usertheta)
  return(isfixed)
}


#' @param filler A function that generates start for theta. Must take a single argument, `n` the number for values to generate.
#' @describeIn t_u2i Convert `usertheta` to a `starttheta` by filling in the `NA` elements.
t_u2s <- function(usertheta, filler = function(n){seq(length.out = n)}){
  starttheta <- usertheta
  isfixed <- t_u2i(usertheta)
  starttheta[!isfixed] <- filler(sum(!isfixed))
  return(starttheta)
}

#' @describeIn t_u2i Convert `usertheta` to a `starttheta` by filling in the `NA` elements with numbers between 0 and 1 from `runif()`.
t_u2s_runif <- function(usertheta){
  t_u2s(usertheta, filler = runif)
}

#' @param c A constant
#' @describeIn t_u2i Convert `usertheta` to a `starttheta` by filling in the `NA` elements with numbers between 0 and 1 from `runif()`.
t_u2s_const <- function(usertheta, c){
  t_u2s(usertheta, filler = function(n){rep(c, n)})
}


#' @describeIn t_u2i Convert `starttheta` and `isfixed` back to a `usertheta` by replacing any non-fixed elements with `NA`.
t_sf2u <- function(starttheta, isfixed){
  stopifnot(all(isfixed %in% c(TRUE, FALSE)))
  stopifnot(length(starttheta) == length(isfixed))
  usertheta <- starttheta
  usertheta[!isfixed] <- NA
  return(usertheta)
}

#' @param fitted Only the fitted elements of `theta`. Must be the same number as `NA` values in `usertheta` or `FALSE` in `isfixed`.
#' @describeIn t_u2i Convert fitted values and `usertheta` to a `starttheta` by replacing any non-fixed elements with the fitted values.
t_fu2t <- function(fitted, usertheta){
  isfixed <- t_u2i(usertheta)
  stopifnot(length(fitted) == sum(!isfixed))
  theta <- usertheta
  theta[!isfixed] <- fitted
  return(theta)
}

#' @describeIn t_u2i Convert `fitted`, `starttheta` and `isfixed` to a `theta` by replacing any non-fixed elements with the fitted values.
t_sfi2u <- function(fitted, starttheta, isfixed){
  stopifnot(length(starttheta) == length(isfixed))
  usertheta <- t_sf2u(starttheta, isfixed)
  theta <- t_fu2t(fitted, usertheta)
  return(theta)
}

#' @describeIn t_u2i Convert `starttheta` into a vector of the just the elements to be fitted.
t_si2f <- function(starttheta, isfixed){
  stopifnot(length(starttheta) == length(isfixed))
  stopifnot(all(isfixed %in% c(TRUE, FALSE)))
  return(starttheta[!isfixed])
}
