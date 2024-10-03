# R6 ADFun class. Making it here because templating of the Cpp class means it is hard to expose via Rcpp
# Using R6 because it has modify in place semantics which match the ADFun objects
# following https://adv-r.hadley.nz/r6.html

#' @name oldADFUNdoc
#' @title A Class for Storing a CppAD Tape (ADFun) Object
#' @description
#' *This is a low level object useful for implementing score matching estimators.*
#' An `R6` class for storing a 'pointer' to a `CppAD` tape in `C++` (also called an `ADFun`) and associated information. 
#' Currently tools for modifying this information are not available in this package, however tools for creating new `ADFun` objects from an existing `ADFun` are available.
#' Typically an `ADFun` object will be created by [`buildsmdtape()`].
#' @field ptr A `Rcpp` external pointer to a `CppAD` `ADFun` object.
#' @param ptr A `Rcpp` external pointer to a `CppAD` `ADFun` object.
#' @field xtape The (numeric) vector of independent variable values used for taping.
#' @param xtape The (numeric) vector of independent variables used for taping.
#' @field dyntape The (numeric) vector of dynamic parameters used for taping.
#' @param dyntape The (numeric) vector of dynamic parameters used for taping.
#' @field usertheta A (numeric) vector of `NA` values and fixed values specifying the parameters of taped function that were considered dynamic parameters or fixed parameters respectively.
#' @param usertheta A (numeric) vector of `NA` values and fixed values specifying the inputs of the taped function that were considered independent variables or dynamic parameters respectively.
#' @field name An easy to read name for the taped function
#' @param name An easy to read name for the taped function
#' @inheritSection buildsmdtape Introduction to CppAD Tapes
#' @inheritSection buildsmdtape Warning: multiple CPU
#' @examples
#' tapes <- buildsmdtape(
#'   "sim", "sqrt", "sph",
#'   ll = "ppi",
#'   ytape =  rep(1/3, 3),
#'   usertheta = ppi_paramvec(p=3),
#'   bdryw = "minsq",
#'   acut = 0.01,
#'   verbose = FALSE)
#' tapes$smdtape$xtape
#' tapes$smdtape$dyntape
#' tapes$smdtape$name
#' tapes$smdtape$ptr
NULL

