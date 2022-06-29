#' @title Functions for converting between a user friendly theta and other forms of theta
#' @param usertheta A parameter vector. NA elements are to be fitted other elements are fixed.
#' @param isfixed A boolean vector same length as the parameter vector. `TRUE` values are fixed at the value of `starttheta`, `FALSE` are fitted.
#' @description Converts a `usertheta` to `isfixed`.
t_u2f <- function(usertheta){
  isfixed <- !is.na(usertheta)
  return(isfixed)
}


#' @param filler A function that generates start for theta. Must take a single argument, `n` the number for values to generate.
#' @describeIn t_u2f Convert `usertheta` to a `starttheta` by filling in the `NA` elements.
t_u2s <- function(usertheta, filler = function(n){seq(length.out = n)}){
  starttheta <- usertheta
  isfixed <- t_u2f(usertheta)
  starttheta[!isfixed] <- filler(sum(!isfixed))
  return(starttheta)
}

#' @describeIn t_u2f Convert `usertheta` to a `starttheta` by filling in the `NA` elements with numbers between 0 and 1 from `runif()`.
t_u2s_runif <- function(usertheta){
  t_u2s(usertheta, filler = runif)
}

#' @param c A constant
#' @describeIn t_u2f Convert `usertheta` to a `starttheta` by filling in the `NA` elements with numbers between 0 and 1 from `runif()`.
t_u2s_const <- function(usertheta, c){
  t_u2s(usertheta, filler = function(n){rep(c, n)})
}


#' @describeIn t_u2f Convert `starttheta` and `isfixed` back to a `starttheta` by replacing any non-fixed elements with `NA`.
t_sf2u <- function(starttheta, isfixed){
  stopifnot(length(starttheta) == length(isfixed))
  usertheta <- startttheta
  usertheta[!isfixed] <- NA
  return(usertheta)
}



