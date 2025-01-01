print.Rcpp_ADFun <- function(x, ...){
  intro <- sprintf("Tape of a function '%s' with %i independent arguments and %i dynamic arguments. Function returns vectors of length %i.\n", x$name, x$domain, x$size_dyn_ind, x$range)

  #print xtape and dyntape
  str_xtape <- paste(format(x$xtape[seq.int(1, length.out = min(x$domain, 3))], ...), collapse = " ")
  if (x$domain > 3){str_xtape <- paste0(str_xtape, "...")}
  str_dyntape <- paste(format(x$dyntape[seq.int(1, length.out = min(x$size_dyn_ind, 3))], ...), collapse = " ")
  if (x$size_dyn_ind > 3){str_dyntape <- paste0(str_dyntape, "...")}

  #help info
  end <- 'For help using this object, run help("Rcpp_ADFun", "scorematchingad").\n'

  allstr <- paste(intro, "xtape is:", str_xtape, "\n dyntape is:", str_dyntape, "\n", end, collapse = " ")
  cat(allstr)
  invisible(allstr)
}

#' @name print,Rcpp_ADFun
#' @title Print or show a summary of an Rcpp_ADFun
#' @description Both `print()` and `show()` will display a summary of a \linkS4class{Rcpp_ADFun} object.
#' @aliases print,Rcpp_ADFun-method
#' @param x An object of class \linkS4class{Rcpp_ADFun}.
#' @param ... Passed to [`format()`].
#' @export
setMethod("print", "Rcpp_ADFun", function(x, ...){
  print.Rcpp_ADFun(x, ...)
})

# need to define show too - this the generic that happens automatically when one types object into the console - and usually it automatically print(), but I guess Rcpp has defined show() for these objects differently.
#' @name print,Rcpp_ADFun
#' @details
#' The `show()` method overrides the default `show()` method for [`Rcpp::C++Object`][Rcpp::C++Object-class] objects from the `Rcpp` package.
#' @aliases show,Rcpp_ADFun-method
#' @param object An object of class \linkS4class{Rcpp_ADFun}.
#' @importFrom methods show
#' @export
setMethod("show", "Rcpp_ADFun", function(object){
  print.Rcpp_ADFun(object)
})
