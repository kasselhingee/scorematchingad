#' @title Build a CppAD Tapes for Score Matching
#' @family tape builders
#' @param thetatape_creator A function that generates tape values for theta. Must take a single argument, `n` the number for values to generate
#' @param verbose If `TRUE` more details are printed when taping.
#' @description
#' Generates `CppAD` tapes for the log-likelihood and score matching objective  for a model specified by name and an optional transformation to a manifold.
#' Three steps are performed by `buildsmotape()`, corresponding to each of the functions `manifoldtransform()`, `tapell()` and `tapesmo()`.
#' @details
#' When using, `CppAD` one first creates *tapes* of functions. These tapes can then be used for evaluating the function and its derivatives, and generating further tapes through argument swapping, differentiation and composition.
#' The taping relies on specifying typical argument values for the functions, so the programming is simplest when the function is defined without conditions.
#' Tapes can have both *independent* variables and *dynamic* variables.
#' The differentiation with `CppAD` occurs with respect to the independent variables.
#' Tapes of tapes are possible, including tapes that swap the independent and dynamic variables - this is how this package differentiates with respect to a dynamic variables.
#'
#' To build a tape for the score matching objective function, the package first tapes the map from a point \eqn{z} on the manifold to the value of the log-likelihood, where the independent variable is the \eqn{z}, and the dynamic variable is a vector of the *non*-fixed parameter values.
#' This tape is then used to generate a tape for the score matching objective with the non-fixed parameter values as the independent variable.
#' # Introduction to CppAD Tapes
#' This package uses version 2022000.2 of the algorithmic differentiation library `CppAD` \insertCite{bell2023cp}{scorecompdir} to build score matching estimators.
#' Full help for `CppAD` can be found at <https://cppad.readthedocs.io/>.
#' 
#' Differentiation proceeds by *taping* the basic (*atomic*) operations performed on the independent variables and dynamic parameters. The atomic operations include multiplication, division, addition, sine, cosine, exponential and many more.
#' Example values for the variables and parameters are used to conduct this taping, so care must be taken with any conditional (e.g. if-then) steps, and `CppAD` has [special tools for this](https://cppad.readthedocs.io/en/latest/CondExp.html).

#' Once the taping is complete an [`ADFun`](https://cppad.readthedocs.io/en/latest/ADFun.html) object is created.
#' This tape can be evaluated, differentiated, used for further taping (see [base2ad](https://cppad.readthedocs.io/en/latest/base2ad.html)), solving differential equations and more.
#' The sequence of operations can also be printed.
#' The differentiation is with respect to the independent variables, although the dynamic parameters can be altered, allowing for swapping independent variables and dynamic parameters.
#' 
#' # Warning: multiple CPU
#' Each time computations such as derivatives are performed the corresponding `C++` object is altered. Parallel use of the same `ADFun` object thus requires care and is not tested. For now I recommend creating a new `ADFun` object for each CPU.

#' @references \insertAllCited{}

#' @return
#' `buildsmotape()` returns a list of the log-likelihood tape (note that the *input* for this tape is a function on the manifold), the tape of the score matching objective (the *input* here is the non-fixed parameter values), and some information used to generate the tape.
#'
#' 
#' @examples
#' p <- 3
#' u <- movMF::rmovMF(1, rep(1, p))
#' ltheta <- p #length of vMF parameter vector
#' intheta <- rep(NA, length.out = ltheta)
#' tapes <- buildsmotape("Snative", "vMF",
#'               ytape = u,
#'               usertheta = intheta,
#'               "ones", 1, verbose = FALSE
#'               )
#' pForward0(tapes$lltape, u, runif(n = ltheta))
#' pForward0(tapes$smotape, runif(n = ltheta), u)
#' @export
buildsmotape <- function(start, tran, man, llname,
                         ytape, usertheta,
                         divweight = "ones", acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE){

  if(all(!is.na(usertheta))){stop("All elements of theta are fixed")}

  if (!(all(c(tran, man) == c("sqrt", "sph")) | all(c(tran, man) == c("identity", "sim")))){
    if (divweight != "ones"){warning("Manifold supplied has no boundary. Using divweight = 'ones' is strongly recommended.")}
  }
  if ((divweight == "ones") && (abs(acut - 1) > 1E-8)){
    warning("The value of 'acut' is ignored for divweight == 'ones'")
  }

  tranman <- manifoldtransform(start, tran, man)
  lltape <- tapell(llname = llname,
                    ytape = ytape,
                    usertheta = usertheta, 
                    thetatape_creator = thetatape_creator,
                    tran = tranman$tran)
  stopifnot(is.numeric(acut))
  smotape <- tapesmo(lltape = lltape,
                        tran = tranman$tran,
                        man = tranman$man,
                        divweight = divweight,
                        acut = acut,
                        verbose = verbose)
  return(list(
    lltape = lltape,
    smotape = smotape,
    info = list(
      name = llname,
      transform = tran,
      manifold = man,
      ulength = length(ytape),
      divweight = divweight,
      acut = acut
    )
  ))
}


