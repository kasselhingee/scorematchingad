#' Generate a tape of a custom unnormalised log-density
#' @family tape builders
#' @param fileORcode A character string giving the path name of a file containing the unnormalised log-density definition *OR* code. `fileORcode` will be treated as a file name if `fileORcode` contains no new line characters ('\\n' or '\\r\\n') and has a file extension detected by [`tools::file_ext()`].
#' @param Cppopt List of named options passed to `Rcpp::sourceCpp()`
#' @param x Value of independent variables for taping.
#' @param theta Value of the dynamic parameter vector for taping.
#' @description
#' Generate tapes of unnormalised log-densities.
#' Use `tape_ult()` to specify a custom unnormalised log-density using `C++` code much like `TMB::compile()`.
#' Use `tape_uld_inbuilt()` for tapes of inbuilt unnormalised log-densities implemented in this package.
#' @details
#' The function `tape_uld()` uses [`Rcpp::sourceCpp()`] to generate a tape of a function defined in C++. 
#' (An alternative design, where the function is compiled interactively and then taped using a function internal to `scorematchingad`, was not compatible with Windows OS).
#'
#' The result result is NOT safe to save or pass to other CPUs in a parallel operation.
#' 
#' # Writing the `fileORcode` Argument
#' The code (possibly in the file pointed to by `fileORcode`) must be `C++` that uses only `CppAD` and `Eigen`, which makes it very similar to the requirements of the input to `TMB::compile()` (which also uses `CppAD` and `Eigen`).
#' 
#' The start of `code` should always be "`a1type fname(const veca1 &x, const veca1 &theta){`" where `fname` is your chosen name of the log-density function, `x` represents a point in the data space and `theta` is a vector of parameters for the log-density. This specifies that the function will have two vector arguments (of type `veca1`) and will return a single numeric value (`a1type`).
#' 
#' The type `a1type` is a double with special ability for being taped by `CppAD`. The `veca1` type is a vector of `a1type` elements, with the vector wrapping supplied by the `Eigen` C++ package (that is an `Eigen` matrix with 1 column and dynamic number of rows).
#' 
#' The body of the function must use operations from Eigen and/or CppAD, prefixed by `Eigen::` and `CppAD::` respectively. 
#' There are no easy instructions for writing these as it is genuine `C++` code, which can be very opaque to those unfamiliar with `C++`.
#' However, recently ChatGPT and claude.ai have been able to very quickly translating `R` functions to `C++` functions (KLH has been telling these A.I. to use Eigen and CppAD, and giving the definitions of `a1type` and `veca1`).
#' I've found the quick reference pages for for [`Eigen`](https://eigen.tuxfamily.org/dox/) useful. Limited unary and binary operations are available directly from [`CppAD`](https://cppad.readthedocs.io) without `Eigen`. 
#' For the purposes of score matching the operations should all be smooth to create a smooth log-density and the normalising constant may be omitted.
#' @examples
#' \dontrun{
#' out <- tape_uld(system.file("demo_custom_uld.cpp", package = "scorematchingad"), 
#'                 rep(0.2, 5), rep(-0.1, 5))
#' out$fun(c(0.1, 0.2, 0.2, 0.2, 0.2), c(-0.5, -0.5, -0.1, -0.1, 0))
#' out$tape$eval(c(0.1, 0.2, 0.2, 0.2, 0.2), c(-0.5, -0.5, -0.1, -0.1, 0))
#' out$tape$Jac(c(0.1, 0.2, 0.2, 0.2, 0.2), c(-0.5, -0.5, -0.1, -0.1, 0))
#' out$tape$name
#' }
#' @return A list of three objects
#' + `fun` a function that evaluates the function directly
#' + `tape` a tape of the function
#' + `file` the temporary file storing the final source code passed to [`Rcpp::sourceCpp()`]
#' @importFrom RcppEigen fastLmPure
#' @export
tape_uld <- function(fileORcode = "", x, theta, Cppopt = NULL){
  #read code
  filelines <- string_to_lines(fileORcode)
  if ((length(filelines) == 1) && (nchar(tools::file_ext(fileORcode)) > 0)){
    filelines <- readLines(fileORcode)
  }

  #remove empty first lines
  blank_lines <- grepl("^\\s*$", filelines)
  filelines <- filelines[1:length(filelines) >= min(which(!blank_lines))]

  # check provided code
  inputsloc <- regexpr("^[[:space:]]*a1type (?<fname>[^[:space:]]+)\\((?<arg1>[^,]+),[[:space:]]*(?<arg2>[^\\)]+)", filelines[1], perl = TRUE)
  starts <- attr(inputsloc, "capture.start")
  lengths <- attr(inputsloc, "capture.length")
  if (lengths[, "fname"] < 1){stop("Could not find log-density name")}
  if (lengths[, "arg1"] < 1){stop("Could not find first argument")}
  if (lengths[, "arg2"] < 1){stop("Could not find second argument")}
  arg1 <- substr(filelines[1], starts[, "arg1"], starts[, "arg1"] + lengths[, "arg1"] - 1)
  if (!grepl("^const +veca1 *&", arg1)){stop(sprintf("First argument should have type 'const veca1 &', instead it is %s", arg1))}
  arg2 <- substr(filelines[1], starts[, "arg2"], starts[, "arg2"] + lengths[, "arg2"] - 1)
  if (!grepl("^const +veca1 *&", arg2)){stop(sprintf("First argument should have type 'const veca1 &', instead it is %s", arg2))}
  fname <- substr(filelines[1], starts[, "fname"], starts[, "fname"] + lengths[, "fname"] - 1)
  if(fname == "tapeld"){stop("Please choose a name that is not 'tapeld' - tapeld is used internally")}


  # set up C++ code
  expandedfile <- tempfile(fileext = ".cpp")

  ## preamble
  cat(c(
"#include <scorematchingad_forward.h>",
"#include <Rcpp.h>",
"#include <scorematchingad.h>",
"#include <likelihoods/likelihoods.hpp>",
"// [[Rcpp::depends(RcppEigen)]]",
"// [[Rcpp::depends(scorematchingad)]]",
"// [[Rcpp::export]]"
), file = expandedfile, append = FALSE, sep = "\n")
  ## add in file
  cat(filelines, file = expandedfile, append = TRUE, sep = "\n")
  
  ## taping addendum
  cat(gsub("__FNAME__", fname, tapingboilerplate), file = expandedfile, append = TRUE, sep = "\n")

  tryCatch({
  funs <- new.env()
  # compile
  compileout <- do.call(Rcpp::sourceCpp, c(list(file = expandedfile, env = funs), Cppopt))},
    error = function(e){stop("Could not sourceCpp() ", expandedfile, ".", e$message)})

  # execute exported tape-generating function
  tapeptr <- funs$tapeld(x, theta)
  
  # return tape and function
  list(fun = funs[[compileout$functions[[1]]]],
       tape = tapeptr,
       file = expandedfile)
}

    
   
  ## boiler plate for taping
tapingboilerplate <- c(
"// [[Rcpp::export]]",
"pADFun tapeld(veca1 & x, veca1 & theta){",
"  CppAD::ADFun<double> tape;",
"  CppAD::Independent(x, theta);",
"  veca1 y(1);",
"  y(0) = __FNAME__(x, theta);",
"  tape.Dependent(x, y);",
"  pADFun out(tape, x, theta, \"__FNAME__\");",
"  return(out);",
"}"
)

string_to_lines <- function(code_string) {
  strsplit(code_string, "\r?\n")[[1]]
}

