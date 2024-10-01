#' Generate a tape of custom unnormalised log-density
#' @description Use Rcpp::sourceCpp to generate a tape of a function defined in C++. The result is NOT safe to save or pass to other CPUs in a parallel operation.
#' @param file A character string giving the path name of a file containing the unnormalised log-density definition.
#' @param Cppopt List of named options passed to `Rcpp::sourceCpp()`
#' @param x Value of `x` for taping.
#' @param theta Value of the parameter vector for taping.
#' @examples
#' tape_uld(system.file("demo_custom_uld.cpp", package = "scorematchingad"), rep(0.2, 5), rep(-0.1, 5))
#' @return A list of the tape and the function itself
#' @export
tape_uld <- function(fileORcode = "", x, theta, Cppopt = NULL){
  #read code
  filelines <- string_to_lines(fileORcode)
  if ((length(filelines) == 1) && (nchar(tools::file_ext(fileORcode)) == 0)){
    filelines <- readLines(fileORcode)
  }

  #remove empty first lines
  blank_lines <- grepl("^\\s*$", filelines)
  filelines <- filelines[1:length(filelines) >= min(which(!blank_lines))]

  # check provided code
  inputsloc <- regexpr("^[[:space:]]*a1type (?<fname>[^[:space:]]+)\\((?<arg1>[^,]+),[[:space:]]*(?<arg2>[^\\)]+)", filelines[1], perl = TRUE)
  starts <- attr(inputsloc, "capture.start")
  lengths <- attr(inputsloc, "capture.length")
  if (lengths[, "fname"] < 1){stop("Could not find log-likelihood name")}
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

cppad_module <- Rcpp::Module("cppad_module", PACKAGE="scorematchingad")
#need to run something from cppad_module to avoid lazy loading. loadModule should work here, but clashes with use of Module for the manifolds module
#otherwise returns of pADFun get error: Error in .getClassesFromCache(Class) : 
# class should be either a character-string name or a class definition
# so onLoad should include a cppad_module$empty() call. See init.R
