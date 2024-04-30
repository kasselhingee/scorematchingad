#' @param tranobj A transform object (of type `Rcpp_transform_ad`), typically created by [`manifoldtransform()`].
#' @param ll The name of an inbuilt improper log-likelihood function to tape (which also specifies the parametric model family) or `C++` code for a custom log-likehood function. The `C++` code has special requirements, similar to `TMB::compile()` - see 'Taping a Custom Log-Likelihood Function' below.
#' @param ytape An example measurement value to use for creating the tapes. In the natural (i.e. `start`) manifold of the log-likelihood function. `ytape` will be converted to the `end` manifold according to the `toM()` method for `tranobj` before taping. 
#' Please ensure that `ytape` is the interior of the manifold, and it is probably best if all components of `tranobj$toM(ytape)` are non-zero.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements will become *dynamic parameters*. Other elements will be fixed at the provided value. The length of `usertheta` must be the correct length for the log-likelihood - __no checking is conducted__.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes. Please ensure that the values filled by `thetatape_creator` lead to plausible parameter vectors for the chosen log-likelihood.
#' @param ... Options passed to [`Rcpp::cppFunction()`] for compiling a custom log-likelihood.
#' @describeIn buildsmdtape Creates a `CppAD` tape of an improper log-likelihood as a function of values on the `end` manifold in `tranobj`. The Jacobian of the associated transformation is used to convert the log-likelihood on the natural manifold `start` of the log-likelihood to the `end` manifold.
#' This conversion is needed to account for the change in measure between the manifolds.
#' @details
#' Inbuilt log-likelihood functions (without normalising constants) are:
#'
#' ```{r, results = "asis", echo = FALSE}
#' cat(paste(" +", llnames), sep = "\n")
#' ```
#'
#' # Taping a Custom Log-Likelihood Function
#' Much care is needed when creating a custom log-likelihood function by passing `C++` code as `ll`.
#' The code must be `C++` that uses only `CppAD` and `Eigen`, which makes it very similar to the requirements of the input to `TMB::compile()` (which also uses `CppAD` and `Eigen`).
#' The start of `code` should always be "`a1type fname(const veca1 &x, const veca1 &theta){`" where `fname` is your chosen name of the log-likelihood function, `x` represents a point in the data space and `theta` is a vector of parameters for the log-likelihood. This specifies that the function will have two vector arguments (of type `veca1`) and will return a single numeric value (`a1type`).
#' 
#' The type `a1type` is a double with special ability for being taped by `CppAD`. The `veca1` type is a vector of `a1type` elements, with the vector wrapping supplied by the `Eigen` C++ package (that is an `Eigen` matrix with 1 column and dynamic number of rows).
#' 
#' The body of the function must use operations from Eigen and/or CppAD, prefixed by `Eigen::` and `CppAD::` respectively. 
#' There are no easy instructions for writing these as it is genuine `C++` code, which can be very opaque to those unfamiliar with `C++`.
#' See the [Eigen documentation](https://eigen.tuxfamily.org/dox/group__QuickRefPage.html) for quick reference to available operations from Eigen. Limited operations are available directly from `CppAD` without `Eigen`: [unary operations](https://cppad.readthedocs.io/latest/unary_standard_math.html) and [binary operations](https://cppad.readthedocs.io/latest/binary_math.html). 
#' For the purposes of score matching the operations should all be smooth to create a smooth log-likelihood and the normalising constant may be omitted.
#' The function [`Rcpp::cppFunction()`] is used to compile the custom likelihood function.
#' Compilation can be very slow.
#'
#' The custom log-likelihood should operate on the starting manifold corresponding to `tran`.
#' It is good practice to check that the taping is giving correct answers by applying [`evaltape()`] to a few known situations.
#' @return 
#' `tapell()` returns an [`ADFun`] object with two additional attributes accessed via `attr()`:  
#'  + `ytape` The value of `ytape`
#'  + `tran` The name of the transform specified in `tranobj`.
#' @examples 
#' maninfo <- manifoldtransform("sim", "sqrt", "sph")
#' ppitape <- tapell(ll = "ppi",
#'                   ytape = c(0.2, 0.3, 0.5),
#'                   usertheta = ppi_paramvec(p = 3), 
#'                   tranobj = maninfo$tran) 
#' evaltape(ppitape, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' evaltape(tapeJacobian(ppitape), 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' 
#' ll <- "a1type dirichlet(const veca1 &u, const veca1 &beta) {
#'   size_t d  = u.size();
#'   a1type y(0.);  // initialize summation at 0
#'   for(size_t i = 0; i < d; i++)
#'   {   y   += beta[i] * log(u[i]);
#'   }
#'   return y;
#' }"
#' dirichlettape <- tapell(ll, 
#'                         ytape = rep(1/3, 3), 
#'                         usertheta = rep(NA, 3),
#'                         tranobj = manifoldtransform("sim", "identity", "sim")$tran)
#' @section Warning: There is limited checking of the inputs.
#' @export
tapell <- function(ll,
                   ytape,
                   usertheta,
                   tranobj,
                   thetatape_creator = function(n){seq(length.out = n)},
                   verbose = FALSE,
                   ...){
  stopifnot(inherits(tranobj, "Rcpp_transform_ad"))
  starttheta <- t_u2s(usertheta, filler = thetatape_creator)
  ztape <- tranobj$toM(ytape) #the value of ytape transformed to the manifold

  # if a canned log-likelihood 
  if (ll %in% llnames){
    llname <- ll
    ll <- getllptr(ll)
  lltape <- ptapell2(ztape, starttheta,
                    llfXPtr = ll, 
                    tran = tranobj,
                    fixedtheta = t_u2i(usertheta),
                    verbose = verbose)
  } else {
  ###### or a custom log-likelihood
    # check provided code
    if (!grepl("^[[:space:]]*a1type", ll)){stop("ll isn't a known log-likelihood function name but C++ code for custom log likelihoods must start with a1type")}
    inputsloc <- regexpr("^[[:space:]]*a1type (?<fname>[^[:space:]]+)\\((?<arg1>[^,]+),[[:space:]]*(?<arg2>[^\\)]+)", ll, perl = TRUE)
    starts <- attr(inputsloc, "capture.start")
    lengths <- attr(inputsloc, "capture.length")
    if (lengths[, "fname"] < 1){stop("Could not find log-likelihood name")}
    if (lengths[, "arg1"] < 1){stop("Could not find first argument")}
    if (lengths[, "arg2"] < 1){stop("Could not find second argument")}
    arg1 <- substr(ll, starts[, "arg1"], starts[, "arg1"] + lengths[, "arg1"] - 1)
    if (!grepl("^const +veca1 *&", arg1)){stop(sprintf("First argument should have type 'const veca1 &', instead it is %s", arg1))}
    arg2 <- substr(ll, starts[, "arg2"], starts[, "arg2"] + lengths[, "arg2"] - 1)
    if (!grepl("^const +veca1 *&", arg2)){stop(sprintf("First argument should have type 'const veca1 &', instead it is %s", arg2))}
    
    # extract function name
    llname <- substr(ll, starts[, "fname"], starts[, "fname"] + lengths[, "fname"] - 1)
    
    # prepare taping code
    code <- paste0("Rcpp::XPtr< CppAD::ADFun<double> > tapecustom(veca1 z, veca1 theta, transform_a1type & tran, Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, bool verbose);",
    "\n\n",
    ll,
    "\n\n",
    "Rcpp::XPtr< CppAD::ADFun<double> > tapecustom(veca1 z, veca1 theta, transform_a1type & tran, Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, bool verbose){
      CppAD::ADFun<double>* out = new CppAD::ADFun<double>;
      *out = tapellcpp(z,
                       theta,", "\n",
"                       *",llname, ",\n",
                       "tran,
                       fixedtheta,
                       verbose);
      
      Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
      return(pout);
    }")
    
    # create taping function, but first check ...
    tapecustom <- Rcpp::cppFunction(
      depends =  c("RcppEigen", "scorematchingad"),
      verbose = verbose,
      code = code,
      ...
    )
    
    # now apply the taping function
    lltape <- tapecustom(ztape, starttheta,
                         tran = tranobj,
                         fixedtheta = t_u2i(usertheta),
                         verbose = verbose)
  }
  
  
  
  
  out <- ADFun$new(ptr = lltape,
                   name = paste(tranobj$name(), llname, sep = "-"),
                   xtape = ztape,
                   dyntape =  as.numeric(starttheta[!t_u2i(usertheta)]),
                   usertheta = as.numeric(usertheta))

  attr(out, "ytape") <- ytape
  attr(out, "tran") <- tranobj$name()
  return(out)
}


llnames <- c(
  "dirichlet",
  "ppi",
  "vMF",
  "Bingham",
  "FB"
)
