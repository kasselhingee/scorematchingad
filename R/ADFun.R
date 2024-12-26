# help for the Rcpp_ADFun objects exposed in src/pADFun.cpp

#' @name Rcpp_ADFun-class
#' @docType class
#' @aliases Rcpp_ADFun
#' @title A Class That Contains CppAD Tapes
#' @description Tapes are a record of operations performed by a function. Tapes can be evaluated, differentiated, and have properties (such as domain and range dimensions). Tapes also have dynamic parameters that can be updated. These classes uses 'reference' semantics, so that changes modify in place and copies all point to the same object (and changes modify that same object).
#'
#'
#' @details `print()` will return some properties of the class. Technically the class name is 'Rcpp_ADFun' (so `inherits(x, "Rcpp_ADFun")` will return `TRUE`) and it is a reference class that connects to `CppAD` tapes in `C++`. Many of the methods available for tapes in `CppAD` are made available here.
#'
#' Tapes cannot be saved from session to session.
#'
#' # Extends
#' Extends class \linkS4class{C++Object} from the `Rcpp` package ([`Rcpp::C++Object-class`]), which is a `reference class`.
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
#' @field xtape The (numeric) vector of independent variable values used for taping.
#' @field dyntape The (numeric) vector of dynamic parameters used for taping.
#' @param name An easy to read name for the tape
#' @exportClass Rcpp_ADFun
#' 
#' @section Methods accessed via `$`:
#' \describe{
#'   \item{\code{new_dynamic(dyn_params)}}{Specify new values for the dynamic parameters.}
#'   \item{\code{forward(order)}}{Perform forward mode evaluation for the specified Taylor coefficient order.}
#'   \item{\code{Jacobian()}}{Evaluate the Jacobian of the function.}
#'   \item{\code{Hessiani(i)}}{Evaluate the Hessian for the \code{i}-th element of the range (where \code{i = 0, 1, ...}).}
#'   \item{\code{Hessian0()}}{Evaluate the Hessian for the first element of the range.}
#'   \item{\code{Hessianw(weights)}}{Evaluate the Hessian for a weighted sum of the range.}
#'   \item{\code{set_check_for_nan(check)}}{Set whether the tape should check for NaN values during computation (only effective if C++ debugging is enabled).}
#'   \item{\code{get_check_for_nan()}}{Return whether the tape is configured to check for NaN values during computation.}
#'   \item{\code{eval(dyn_params)}}{Evaluate the function with new dynamic parameters.}
#'   \item{\code{Jac(dyn_params)}}{Compute the Jacobian with new dynamic parameters.}
#'   \item{\code{Hes(dyn_params)}}{Compute the Hessian with new dynamic parameters.}
#'   \item{\code{parameter(index)}}{Check if the \code{index}-th component of the range corresponds to a constant parameter.}
#' }
#'
#' @section Properties:
#' \describe{
#'   \item{\code{size_order}}{Number of Taylor coefficient orders, per variable and direction, currently calculated and stored.}
#'   \item{\code{domain}}{Dimension of the domain space (i.e., length of the independent variables vector).}
#'   \item{\code{range}}{Dimension of the range space.}
#'   \item{\code{size_dyn_ind}}{Number of independent dynamic parameters.}
#'   \item{\code{name}}{An optional name for the tape.}
#'   \item{\code{xtape}}{(Read-only) The values of the independent variables used for taping.}
#'   \item{\code{dyntape}}{(Read-only) The values of the dynamic variables used for taping.}
#' }
#'
#' @examples
#' tape <- tape_uld_inbuilt("dirichlet", c(0.1, 0.4, 0.5), c(-0.5, -0.4, -0.2))
#' # Properties
#' tape$domain
#' tape$range
#' tape$size_dyn_ind
#' tape$name
#' tape$xtape
#' tape$dyntape
#' tape$size_order
#'
#' # Convenient evaluation
#' tape$eval(x = c(0.2, 0.3, 0.5), dyn = c(-0.1, -0.1, -0.5))
#' tape$Jac(x = c(0.2, 0.3, 0.5), dyn = c(-0.1, -0.1, -0.5))
#' matrix(tape$Hes(x = c(0.2, 0.3, 0.5), dyn = c(-0.1, -0.1, -0.5)), nrow = tape$domain)
#' ADFun$help(eval)
#' 
#' # Further methods
#' tape$forward(order = 1, x = c(0.2, 0.3, 0.5))
#' tape$Jacobian(x = c(0.2, 0.3, 0.5))
#' tape$Hessiani(x = c(0.2, 0.3, 0.5), index = 0)
#' tape$Hessian0(x = c(0.2, 0.3, 0.5))
#' tape$Hessianw(x = c(0.2, 0.3, 0.5), w = c(2))
#' tape$new_dynamic(dyn = c(-0.1, -0.1, -0.5))
#' tape$parameter(0)
#' tape$set_check_for_nan(FALSE)
#' tape$get_check_for_nan()


NULL

Rcpp::loadModule("ADFun", TRUE)

