#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @importFrom methods formalArgs
#' @importFrom stats cov runif var weighted.mean
#' @importFrom utils installed.packages packageDescription tail
#' @useDynLib scorecompdir
#' @details
#' This package's major contributions are
#'  * A general capacity to implement score matching estimators that use algorithmic differentiation to avoid tedious manual algebra.
#' The package uses `CppAD` and `Eigen` to differentiate model densities and compute the score matching objective (the \eqn{\hat\psi(f, f_0)} in __Score Matching__ section below).
#' The score matching objective is minimised by a modern implementation of conjugate gradient descent available through [`Rcgmin::Rcgmin()`].
#'  * Score matching estimators for the Polynomially-Tilted Pairwise Interaction (PPI) model \insertCite{scealy2022sc}{scorecompdir}. See function [`ppi()`].
#'  * Score matching and hybrid score matching estimators for von Mises Fisher and Bingham directional distributions \insertCite{mardia2016sc}{scorecompdir}.
#'  * Implementation of a modification of Windham's robustifying method \insertCite{windham1995ro}{scorecompdir} for many exponential family distributions. See [`Windham()`].
#' For some models the density approaches infinity at some locations, creating difficulties for the weights in Windham's original method \insertCite{@scealy2022ro}{scorecompdir}.
#' \insertNoCite{*}{scorecompdir}
#' @references
#' \insertAllCited{}
"_PACKAGE"

