# Wrappers for numerical optimisers and related info

#### Fixed Point Iteration ####
fp <- function(...){
  args = list(...)
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
  if (length(Rcgmincontrol) == 0){Rcgmincontrol <- NULL}
  fpcontrol <- control[names(control) %in% FixedPoint_controlnames, drop = FALSE]
  if (length(fpcontrol) == 0){fpcontrol <- NULL}
  return(list(
    fp = fpcontrol,
    Rcgmin = Rcgmincontrol))
}
