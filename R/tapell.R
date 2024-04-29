#' @param tranobj A transform object (of type `Rcpp_transform_ad`), typically created by [`manifoldtransform()`].
#' @param ll The name of an inbuilt improper log-likelihood function to tape (which also specifies the parametric model family) or a custom log-likehood function created by [`customll()`]. The `ll` should operate on the untransformed (i.e. starting) manifold.
#' @param ytape An example measurement value to use for creating the tapes. In the natural (i.e. `start`) manifold of the log-likelihood function. `ytape` will be converted to the `end` manifold according to the `toM()` method for `tranobj` before taping. 
#' Please ensure that `ytape` is the interior of the manifold, and it is probably best if all components of `tranobj$toM(ytape)` are non-zero.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements will become *dynamic parameters*. Other elements will be fixed at the provided value. The length of `usertheta` must be the correct length for the log-likelihood - __no checking is conducted__.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes. Please ensure that the values filled by `thetatape_creator` lead to plausible parameter vectors for the chosen log-likelihood.
#' @describeIn buildsmdtape Creates a `CppAD` tape of an improper log-likelihood as a function of values on the `end` manifold in `tranobj`. The Jacobian of the associated transformation is used to convert the log-likelihood on the natural manifold `start` of the log-likelihood to the `end` manifold.
#' This conversion is needed to account for the change in measure between the manifolds.
#' @details
#' Currently available improper log-likelihood functions are:
#'
#' ```{r, results = "asis", echo = FALSE}
#' cat(paste(" +", llnames), sep = "\n")
#' ```
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
#' @section Warning: There is limited checking of the inputs.
#' @export
tapell <- function(ll,
                   ytape,
                   usertheta,
                   tranobj,
                   thetatape_creator = function(n){seq(length.out = n)},
                   verbose = FALSE){
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
    
    # create taping function
    tapecustom <- Rcpp::cppFunction(
      depends =  c("RcppEigen", "scorematchingad", "Rcpp"),
      includes = c("#include <utils/PrintFor.hpp>", "#include <tapellcpp.h>", "#include <manifoldtransforms/transforms.hpp>", "#include <utils/wrapas.hpp>"),
      verbose = verbose,
      code = code
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
