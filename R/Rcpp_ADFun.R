# help for the Rcpp_ADFun objects exposed in src/pADFun.cpp

#' @name Rcpp_ADFun-class
#' @docType class
#' @aliases Rcpp_ADFun ADFun
#' @title A Class for CppAD Tapes
#' @description Objects of type `Rcpp_ADFun` contain a tape of a `C++` function (which has class `ADFun` in `CppAD`). These tapes are a record of operations performed by a function. Tapes can be evaluated and differentiated, and have properties (such as domain and range dimensions). Tapes also have dynamic parameters that can be updated. This class, `Rcpp_ADFun` uses reference semantics, so that copies all point to the same object and changes modify in place (i.e. changes modify the same object).
#' Properties and methods of an `Rcpp_ADFun` object are accessed via `$`.
#'
#' @details 
#' An object of class `Rcpp_ADFun` wraps an `ADFun` object from `CppAD`. Many of the properties and behaviour of an `Rcpp_ADFun` object come directly from `ADFun` objects so more details and context can be found by looking at the `ADFun` object help in the `CppAD` [`help`](https://cppad.readthedocs.io).
#' The methods `eval()`, `Jac()` and `Hes()` have been added by `scorematchingad` as there were many cases where this seemed like an easier way to evaluate a tape.
#' 
#' Default printing of an `Rcpp_ADFun` object gives a short summary of the object, see [`print,Rcpp_ADFun`].
#'
#' Tapes cannot be saved from session to session.
#' 
#' # Methods - Tape Properties:
#' + `$size_order` Number of Taylor coefficient orders, per variable and direction, currently calculated and stored in the object.
#' + `$domain` Dimension of the domain space (i.e., length of the independent variables vector `x`).
#' + `$range` Dimension of the range space (i.e., length of the vector returned by `$eval()`).
#' + `$size_dyn_ind` Number of independent dynamic parameters (i.e., length of the vector of dynamic parameters `dyn`).
#' + `$name` A name for the tape (may be empty). This is yet to incorporate the `CppAD` `function_name` property.
#' + `$xtape` The values of the independent variables used for the initial taping.
#' + `$dyntape` The values of the dynamic parameters used for the initial taping.
#' + `$get_check_for_nan()` Debugging: Return whether the tape is configured to check for NaN values during computation. The check for NaN only occurs if the `C++` compilation enables debugging.
#' + `$set_check_for_nan(bool)` Set whether the tape should check for NaN values during computation (only effective if C++ debugging is enabled).
#' + `$parameter(i)` Check if the `i`th component of the range corresponds to a constant parameter. Indexing is by the `C++` default, that is the first component has index `0`, the last component has index `$range - 1`.
#' + `$new_dynamic(dyn)` Specify new values for the dynamic parameters.
#'
#' # Methods - Tape Evaluation:
#' + `$eval(x, dyn)` Evaluate the function at new values of the variables and dynamic parameters. Returns a vector of length `$range`.
#' + `$Jac(x, dyn)` Compute the Jacobian at new values of the variables and dynamic parameters. Returns a vector of length `$range * $domain` arranged so that the first `$domain` elements correspond to the gradient of the first element of the range. The next `$domain` elements correspond to the gradient of the second element of the range, and so on.
#' + `$Hes(x, dyn)` Compute the Hessian of the first element of the range at new values of the variables and dynamic parameters. Returns a vector of length `$domain * $domain` where the `j*n + l` element corresponds to differentiating with respect to the `l`th element of the domain, then with respect to the `j`th element of the domain, with `n` the size of the domain.
#' + `$Jacobian(x)` Evaluate the Jacobian of the function at the current set of dynamic parameters.
#' + `$Hessiani(x, i)` Evaluate the Hessian for the \code{i}-th element of the range (where \code{i = 0, 1, ...}). Returns a vector arranged the same as `$Hes()`.
#' + `$Hessian0(x)` Evaluate the Hessian for the first element of the range (like `$Hes()` but uses the current values of the dynamic parameters). Returns a vector arranged the same as `$Hes()`.
#' + `$Hessianw(x, w)` Evaluate the Hessian for a weighted sum of the range. Returns a vector arranged the same as `$Hes()`.
#' + `$forward(q, x)` Perform forward mode evaluation for the specified Taylor coefficient order `q`. See the [`CppAD` help](https://cppad.readthedocs.io) for more.
#'
#' # Method Arguments
#' + `x` A vector of independent variables.
#' + `dyn` A vector of dynamic parameters.
#' + `q` Taylor coefficient order for evaluating derivatives with `$forward()`.
#' + `i` Index of range result. `i = 0, 1, ..., $range - 1`.
#' + `bool` Either `TRUE` or `FALSE` to set `check_for_nan` behaviour using `$set_check_for_nan()`.
#' + `w` Weights assigned to each element of the range, for use with `$Hessianw()`.
#'
#'
#' # Extends
#' Extends class `C++Object` from the `Rcpp` package ([`Rcpp::C++Object-class`]), which is a `reference class`.
#' For those familiar with `C++`, an object of class `Rcpp_ADFun` contains a pointer to a `CppAD` `ADFun` object. 
#'
#' # Introduction to CppAD Tapes
#' This package uses version 2024000.5 of the algorithmic differentiation library `CppAD` \insertCite{bell2023cp}{scorematchingad} to build score matching estimators.
#' Full help for `CppAD` can be found at <https://cppad.readthedocs.io/>.
#' 
#' When using `CppAD` one first creates a *tape* of the basic (*atomic*) operations of a function.
#' The atomic operations include multiplication, division, addition, sine, cosine, exponential and many more.
#' These tapes can then be used for evaluating the function and its derivatives, and generating further tapes through argument swapping, differentiation and composition (see for example [`tape_swap()`] and [`tape_Jacobian()`]).
#' Tapes can have both *independent* variables and *dynamic* parameters, and the differentiation occurs with respect to the independent variables.
#' The atomic operations within a function are taped by following the function evaluation on example values for the variables and parameters, so care must be taken with any conditional (e.g. if-then) operations, and [`CppAD`](https://cppad.readthedocs.io/) has a special tool for this called `CondExp` (short for `conditional expressions`).
#'
#' The result of taping, called a *tape*, is exposed as an object of class [`Rcpp_ADFun`], which contains a `CppAD` `ADFun` object.
#' Although the algorithmic differentiation is with respect to the independent variables, a new tape (see [`tape_swap()`]) can be created where the dynamic parameters become independent variables.
#' For the purposes of score matching, there are also *fixed* parameters, which are the elements of the model's parameter vector that are given and not estimated.
#'
#' The example values used for taping are saved in the `$xtape` and `$dyntape` properties of [`Rcpp_ADFun`] objects.
#'
#' # Warning: multiple CPU
#' Each time a tape is evaluated the corresponding `C++` object is altered. Parallel use of the same `ADFun` object thus requires care and is not tested. For now I recommend creating a new `ADFun` object for each CPU.
#'
#' # Improvements
#' A few methods for `CppAD` `ADFun` objects are not yet available through `Rcpp_ADFun` objects. These ones would be nice to include:
#' + `optimize()`
#' + `function_name_set()` and `function_name_get()` working with `$name`
#' + `Reverse()`
#'
#' @exportClass Rcpp_ADFun
#' 
#'
#' @examples
#' tape <- tape_uld_inbuilt("dirichlet", c(0.1, 0.4, 0.5), c(-0.5, -0.4, -0.2))
#' # Convenient evaluation
#' tape$eval(x = c(0.2, 0.3, 0.5), dyn = c(-0.1, -0.1, -0.5))
#' tape$Jac(x = c(0.2, 0.3, 0.5), dyn = c(-0.1, -0.1, -0.5))
#' matrix(tape$Hes(x = c(0.2, 0.3, 0.5), dyn = c(-0.1, -0.1, -0.5)), nrow = tape$domain)
#' 
#' # Properties
#' tape$domain
#' tape$range
#' tape$size_dyn_ind
#' tape$name
#' tape$xtape
#' tape$dyntape
#' tape$size_order
#' tape$new_dynamic(dyn = c(-0.1, -0.1, -0.5))
#' tape$parameter(0)
#' tape$set_check_for_nan(FALSE)
#' tape$get_check_for_nan()
#'
#' # Further methods
#' tape$forward(order = 0, x = c(0.2, 0.3, 0.5))
#' tape$Jacobian(x = c(0.2, 0.3, 0.5))
#' tape$Hessiani(x = c(0.2, 0.3, 0.5), i = 0)
#' tape$Hessian0(x = c(0.2, 0.3, 0.5))
#' tape$Hessianw(x = c(0.2, 0.3, 0.5), w = c(2))


NULL

Rcpp::loadModule("ADFun", TRUE)

