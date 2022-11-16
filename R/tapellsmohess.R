#' @title Tape a log-likelihood
#' @param manifoldname Name of manifold and transformation to the manifold
#' @param llname The name of the log-likelihood function to tape
#' @param utape An example measurement value to use for creating the tape. In the natural manifold of the log-likelihood function. `utape` will be converted to the manifold according to `manifoldname` before taping.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements will become *dynamic parameters*. Other elements will be fixed at the provided value.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes.
#' @description Creates a `CppAD` tape of the log-likelihood as a function of values on the provided manifold. The log-likelihood accounts for the change of measure between the natural manifold of the log-likelihood and the manifold with transformation provided by `manifoldname`.
#' @return 
#' An `Rcpp::XPrt` object that points to a `CppAD::ADFun` object.
#' The object will also have attributes:
#'  + `fname` A name of the taped function, derived from `llname`
#'  + `utape` The value of `utape`
#'  + `dyntape` The value of non-fixed elements of `usertheta` used for taping.
#'  + `usertheta` The vector `usertheta`.
#' fname, utape, dyntape
#' @export
tapell <- function(llname,
                   utape,
                   usertheta,
                   manifoldname,
                   thetatape_creator = function(n){seq(length.out = n)}){
  return(NULL) 
}



