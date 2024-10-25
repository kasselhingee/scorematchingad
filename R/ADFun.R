# help for the Rcpp_ADFun objects exposed in src/pADFun.cpp

#' @name ADFun
#' @aliases Rcpp_ADFun
#' @title A Class That Contains CppAD Tapes
#' @description Tapes are a record of operations performed by a function. Tapes can be evaluated, differentiated, and have properties (such as domain and range dimensions). Tapes also have dynamic parameters that can be updated. These classes uses 'reference' semantics, so that changes modify in place and copies all point to the same object (and changes modify that same object).
#'
#'
#' @details `print()` will return some properties of the class. Technically the class name is 'Rcpp_ADFun' (so `inherits(x, "Rcpp_ADFun")` will return `TRUE`) and it is a reference class that connects to `CppAD` tapes in `C++`. Many of the methods available for tapes in `CppAD` are made available here.
#'
#' Tapes cannot be saved from session to session.
#'
#' # Introduction to CppAD Tapes
#' This package uses version 2024000.5 of the algorithmic differentiation library `CppAD` \insertCite{bell2023cp}{scorematchingad} to build score matching estimators.
#' Full help for `CppAD` can be found at <https://cppad.readthedocs.io/>.
#' 
#' Differentiation proceeds by *taping* the basic (*atomic*) operations performed on the independent variables and dynamic parameters. The atomic operations include multiplication, division, addition, sine, cosine, exponential and many more.
#' Example values for the variables and parameters are used to conduct this taping, so care must be taken with any conditional (e.g. if-then) operations, and [`CppAD`](https://cppad.readthedocs.io/) has a special tool for this called `CondExp` (short for `conditional expressions`).
#'
#' The result of taping is an object of class `ADFun` in `CppAD` and is often called a *tape*.
#' This `ADFun` object can be evaluated, differentiated, used for further taping (via `CppAD`'s `base2ad()`), solving differential equations and more.
#' The differentiation is with respect to the independent variables, however the dynamic parameters can be altered which allows for creating a new `ADFun` object where the dynamic parameters become independent variables (see [`tapeSwap()`]).
#' For the purposes of score matching, there are also *fixed* parameters, which are the elements of the model's parameter vector that are given and not estimated.
#' 
#' # Warning: multiple CPU
#' Each time a tape is evaluated the corresponding `C++` object is altered. Parallel use of the same `ADFun` object thus requires care and is not tested. For now I recommend creating a new `ADFun` object for each CPU.
#'
#' Some further help is available by `ADFun$help()`.
#' @param x A vector of independent variables.
#' @param dyn A vector of dynamic parameters.
#' @param q Differentiation order.
#' @param i Index of range result.
#' @field domain The number of independent variables (i.e. dimension of Euclidean domain space)
#' @field eval Evaluation of the function at `x` given new values of the dynamic parameters `dyn`. Has two arguments, `x` and `dyn`.
#' @export ADFun

NULL
