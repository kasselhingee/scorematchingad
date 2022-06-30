# Wrappers for numerical optimisers and related info
#' @title Default control parameters for numerical optimisers
#' @description This package uses `Rcgmin()` in many of the score matching estimates and
#' `default_Rcgmin()` returns the default control values used throughout the package.
#' @details
#' The default for `Rcgmin()` is `list(tol = 1E-20)`, which means the optimisation won't end until the squared size of the gradient at the estimate is less than 1E-20.
#' @export
default_Rcgmin <- function(){
  list(tol = 1E-20)
}
#' @describeIn default_Rcgmin Default control parameters for `FixedPoint()`, used in robust estimation.
#' @export
default_FixedPoint <- function(){
  list(Method = "VEA", ConvergenceMetricThreshold = 1E-10)
}


#### Fixed Point Iteration ####
fp <- function(...){
  args = list(...)
  if (is.null(args$control)){args$control <- default_FixedPoint()} #default control
  args <- c(args, args$control) #put control arguments into the highest level of the list
  args2 <- args[names(args) %in% formalArgs(FixedPoint::FixedPoint), drop = FALSE]
  do.call(FixedPoint::FixedPoint, args2)
}

FixedPoint_controlnames <- setdiff(formalArgs(FixedPoint::FixedPoint), c("Function", "Inputs"))


#### Root Finding ####

Rcgmin_controlnames <- c("maxit", "trace", "eps", "dowarn", "tol", "checkgrad", "checkbounds")


#### Function to split control parameters between Rcgmin and FixedPoint ####
splitcontrol <- function(control){
  Rcgmincontrol <- control[names(control) %in% Rcgmin_controlnames, drop = FALSE]
  if (length(Rcgmincontrol) == 0){Rcgmincontrol <- default_Rcgmin()}
  fpcontrol <- control[names(control) %in% FixedPoint_controlnames, drop = FALSE]
  if (length(fpcontrol) == 0){fpcontrol <- default_FixedPoint()}
  return(list(
    fp = fpcontrol,
    Rcgmin = Rcgmincontrol))
}
