#' @title Build CppAD Tapes for Score Matching
#' @family tape builders
#' @param thetatape_creator A function that generates tape values for theta. Must take a single argument, `n` the number for values to generate
#' @param verbose If `TRUE` more details are printed when taping.
#' @description
#' The function `buildsmdtape()` generates `CppAD` tapes (called `ADFun`) for the log-likelihood (without normalising constant) and the score matching discrepancy function \eqn{A(z) + B(z) + C(z)} from `vignette("scorematchingintro")` for a specified model faimly.
#' Three steps are performed by `buildsmdtape()`: first an object that specifies the manifold and any transformation to another manifold is created using [`manifoldtransform()`]; then a tape of the log-likelihood (without normalising constant) is created using [`tapell()`]; finally a tape of \eqn{A(z) + B(z) + C(z)} is created using [`tapesmd()`].
#' @details
#' The model log-likelihood without normalising constant must be implemented in `C++` and is selected by name. Similarly the transforms of the manifold must be implemented in `C++` and selected by name.
#'
#' When using, `CppAD` one first creates *tapes* of functions. These tapes can then be used for evaluating the function and its derivatives, and generating further tapes through argument swapping, differentiation and composition.
#' The taping relies on specifying typical argument values for the functions, so the programming is simplest when the function is defined without conditions.
#' Tapes can have both *independent* variables and *dynamic* variables.
#' The differentiation with `CppAD` occurs with respect to the independent variables.
#' Tapes of tapes are possible, including tapes that swap the independent and dynamic variables - this is how this package differentiates with respect to a dynamic variables (see [`swapDynamic()`]).
#'
#' To build a tape for the score matching discrepancy function, the package first tapes the map from a point \eqn{z} on the `end` manifold to the value of the improper log-likelihood, where the independent variable is the \eqn{z}, and the dynamic variable is a vector of the parameters to estimate, and remaining parameters are fixed and not estimated.
#' This tape is then used to generate a tape for the score matching discrepancy with the parameters to estimate are the independent variable.
#' 
#' # Introduction to CppAD Tapes
#' This package uses version 2022000.2 of the algorithmic differentiation library `CppAD` \insertCite{bell2023cp}{scorecompdir} to build score matching estimators.
#' Full help for `CppAD` can be found at <https://cppad.readthedocs.io/>.
#' 
#' Differentiation proceeds by *taping* the basic (*atomic*) operations performed on the independent variables and dynamic parameters. The atomic operations include multiplication, division, addition, sine, cosine, exponential and many more.
#' Example values for the variables and parameters are used to conduct this taping, so care must be taken with any conditional (e.g. if-then) operations, and `CppAD` has [special tools for this](https://cppad.readthedocs.io/en/latest/CondExp.html).

#' The result of taping is an [`ADFun`](https://cppad.readthedocs.io/en/latest/ADFun.html) object, often called a *tape*.
#' This `ADFun` object can be evaluated, differentiated, used for further taping (see [base2ad](https://cppad.readthedocs.io/en/latest/base2ad.html)), solving differential equations and more.
#' The differentiation is with respect to the independent variables, however the dynamic parameters can be altered which allows for creating a new `ADFun` object where the dynamic parameters have become independent variables (see [`tapeSwap()`]).
#' For the purposes of score matching, there is a third type of parameter: the *fixed* parameters, which are the elements of the model's parameter vector that are given and not estimated.
#' 
#' # Warning: multiple CPU
#' Each time computations such as derivatives are performed the corresponding `C++` object is altered. Parallel use of the same `ADFun` object thus requires care and is not tested. For now I recommend creating a new `ADFun` object for each CPU.

#' @references \insertAllCited{}

#' @return
#' `buildsmdtape()` returns a list of:
#'   + an [`ADFun`] object containing a log-likelihood tape with \eqn{z} on the `end` manifold as the independent variable
#'   + an [`ADFun`] object containing a tape of the score matching discrepancy with the non-fixed parameter values as the independent variables
#'   + some information used to generate the tapes
#'
#' 
#' @examples
#' p <- 3
#' u <- movMF::rmovMF(1, rep(1, p))
#' ltheta <- p #length of vMF parameter vector
#' intheta <- rep(NA, length.out = ltheta)
#' tapes <- buildsmdtape("sph", "identity", "sph", "vMF",
#'               ytape = u,
#'               usertheta = intheta,
#'               "ones", verbose = FALSE
#'               )
#' evaltape(tapes$lltape, u, runif(n = ltheta))
#' evaltape(tapes$smdtape, runif(n = ltheta), u)
#' @export
buildsmdtape <- function(start, tran, end, llname,
                         ytape, usertheta,
                         bdryw = "ones", acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE){

  if(all(!is.na(usertheta))){stop("All elements of theta are fixed")}

  if (!(all(c(tran, end) == c("sqrt", "sph")) | all(c(tran, end) == c("identity", "sim")))){
    if (bdryw != "ones"){warning("Manifold supplied has no boundary. Using bdryw = 'ones' is strongly recommended.")}
  }
  if ((bdryw == "ones") && (abs(acut - 1) > 1E-8)){
    warning("The value of 'acut' is ignored for bdryw == 'ones'")
  }

  tranman <- manifoldtransform(start, tran, end)
  lltape <- tapell(llname = llname,
                    ytape = ytape,
                    usertheta = usertheta, 
                    thetatape_creator = thetatape_creator,
                    tran = tranman$tran)
  stopifnot(is.numeric(acut))
  smdtape <- tapesmd(lltape = lltape,
                        tran = tranman$tran,
                        man = tranman$man,
                        bdryw = bdryw,
                        acut = acut,
                        verbose = verbose)
  return(list(
    lltape = lltape,
    smdtape = smdtape,
    info = list(
      name = llname,
      transform = tran,
      manifold = end,
      ulength = length(ytape),
      bdryw = bdryw,
      acut = acut
    )
  ))
}


