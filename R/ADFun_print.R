#' @export
print.Rcpp_ADFun <- function(x, ...){
  intro <- sprintf("Tape of a function '%s' with %i independent arguments and %i dynamic arguments. Function returns vectors of length %i.\n", x$name, x$domain, x$size_dyn_ind, x$range)

  #print xtape and dyntape
  str_xtape <- paste(format(x$xtape[seq.int(1, length.out = min(x$domain, 3))], ...), collapse = " ")
  if (x$domain > 3){str_xtape <- paste0(str_xtape, "...")}
  str_dyntape <- paste(format(x$dyntape[seq.int(1, length.out = min(x$size_dyn_ind, 3))], ...), collapse = " ")
  if (x$size_dyn_ind > 3){str_dyntape <- paste0(str_dyntape, "...")}

  allstr <- paste(intro, "xtape is:", str_xtape, "\n dyntape is:", str_dyntape, "\n", collapse = " ")
  cat(allstr)
  invisible(allstr)
}
