#' @name tape_smd
# @param tranobj A transform object (of type `Rcpp_transform_ad`), typically created by [`manifoldtransform()`].
#' @param ll The name of an inbuilt unnormalised log-density function (which also specifies the parametric model family) or the tape created by [`tape_uld()`] of a custom unnormalised log-density function with respect to the starting manifold metric.
#' @param ytape An example measurement value to use for creating the tapes. In the natural (i.e. `start`) manifold of the log-likelihood function. 
# `ytape` will be converted to the `end` manifold according to the `toM()` method for `tranobj` before taping. 
#' Please ensure that `ytape` is the interior of the manifold and non-zero.
# and it is probably best if all components of `tranobj$toM(ytape)` are non-zero.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements will become *dynamic parameters*. Other elements will be fixed at the provided value. The length of `usertheta` must be the correct length for the log-likelihood - __no checking is conducted__.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes. Please ensure that the values filled by `thetatape_creator` lead to plausible parameter vectors for the chosen log-likelihood.
#' @details
#' Currently available improper log-likelihood functions are:
#'
#' ```{r, results = "asis", echo = FALSE}
#' cat(paste(" +", llnames), sep = "\n")
#' ```
#' @section Warning: There is no checking of the inputs `ytape` and `usertheta`.
NULL

tapell <- function(ll,
                   ytape,
                   usertheta,
                   tranobj,
                   thetatape_creator = function(n){seq(length.out = n)}){
  stopifnot(inherits(tranobj, "Rcpp_transform_ad"))

  starttheta <- t_u2s(usertheta, filler = thetatape_creator)
  ztape <- tranobj$toM(ytape) #the value of ytape transformed to the manifold

  # choose between a canned log-likelihood or a custom log-likelihood
  if (typeof(ll) == "character"){
    ll <- tape_uld_inbuilt(ll, ytape, starttheta)
  } else if (!inherits(ll, "Rcpp_ADFun")){
    stop("ll must be a name or a taped unnormalised log-density")
  }

  uld_dynfixed <- fixdynamic(ll, starttheta, fixedtheta = t_u2i(usertheta)) # fix some of the model parameters if applicable
  uld_reembed <- reembed(uld_dynfixed, tran = tranobj) #change the underlying metric of the manifold by using a different isometric embedding
  return(uld_reembed)
}


llnames <- c(
  "dirichlet",
  "ppi",
  "vMF",
  "Bingham",
  "FB"
)
