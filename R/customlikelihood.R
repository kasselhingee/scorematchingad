# functions for building custom likelihoods: build, evaluate (to check it). (Then use tapell and check it more!)
#' @title Compile a custom log-likelihood function.
#' @param cacheDir passed to [`Rcpp::cppFunction()`].
#' @param rebuild passed to [`Rcpp::cppFunction()`].
#' @param showOutput passed to [`Rcpp::cppFunction()`].
#' @param verbose passed to [`Rcpp::cppFunction()`].
#' @details
#' Can ommit normalising constant.
#' Like TMB, which also uses CppAD with CppAD
#' Must have two inputs and one output. That is first line must look like:
#' 
#' Uses `RcppXPtrUtils::cppXPtr()`
#' 
#' a1type, veca1
#' 
#' It is good practice to check the function is behaving properly using [`evalll()`].
#' @examples
#' 
#' myll <- customll("a1type dirichlet(const veca1 &u, const veca1 &beta) {
#'   size_t d  = u.size();
#'   a1type y(0.);  // initialize summation at 0
#'   for(size_t i = 0; i < d; i++)
#'   {   y   += beta[i] * log(u[i]);
#'   }
#'   return y;
#' }")
#' evalll(myll, rep(1/3, 3), rep(-0.5, 3))

#' @returns An `adloglikelood` object (which is just an `externalptr` with attributes) for the compiled log-likelihood function. The returned object has an attribute `fname`.
#' @export
customll <- function(code, rebuild = FALSE, 
                     cacheDir = getOption("rcpp.cache.dir", tempdir()), 
                     showOutput = verbose, verbose = getOption("verbose")){
  ptr <- RcppXPtrUtils::cppXPtr(code, depends = c("RcppEigen", "scorematchingad"), rebuild = rebuild, cacheDir = cacheDir, showOutput = showOutput, verbose = verbose)
  RcppXPtrUtils::checkXPtr(ptr, type = "a1type", args = c("const veca1&", "const veca1&"))
  ptr <- new_adloglikelihood(ptr)
  return(ptr)
}

# adloglikelood class things
new_adloglikelihood <- function(ptr, fname = NULL){
  if (is.null(attr(ptr, "fname", exact = TRUE)) && !is.null(fname)){
    attr(ptr, "fname") <- fname
  } else if (!is.null(attr(ptr, "fname", exact = TRUE)) && !is.null(fname)){
    warning("overwriting objects fname attribute")
    attr(ptr, "fname") <- fname
  } else if (is.null(attr(ptr, "fname", exact = TRUE)) && is.null(fname)){
    stop("object needs an fname")
  }
  class(ptr) <- c("adloglikelihood", class(ptr))
  return(ptr)
}
validate_adloglikelihood <- function(ll){
  ll_a <- unclass(ll)
  if (typeof(ll_a) != externalptr){
    stop("adloglikelihood is not an 'externalptr'")
  }
  if (is.null(attr(ll_a, "fname", exact = TRUE))){
    stop("adloglikelihood objects must have a fname attribute")
  }
  return(ll)
}
print.adloglikelihood <- function(x, ...){
  print(sprintf("log-likelihood function '%s' for automatic differentiation", attr(x, "fname", exact = TRUE)))
}
