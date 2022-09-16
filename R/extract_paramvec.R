#' @title Try to extract estimated parameter vector from the result of an estimator
#' @param estobj The output from an estimator in this package.
#' @description If `estobj` is a list then looks first in the slot `est$paramvec` then the first element of the list. If the output is a numeric vector then it will use this. If none of these return a numeric vector then an error will be flagged.
extract_paramvec <- function(estobj){
    if (is.list(estobj)){
      estparamvec <- try(estobj$est$paramvec)
      if (is.numeric(estparamvec) & is.vector(estparamvec)){return(estparamvec)}
      estparamvec <- estobj[[1]]
      if (is.numeric(estparamvec) & is.vector(estparamvec)){return(estparamvec)}
    }
    estparamvec <- estobj
    if (is.numeric(estparamvec) & is.vector(estparamvec)){return(estparamvec)}
    stop("Could not detect estimated parameter vector")
}
