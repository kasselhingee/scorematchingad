#' @param ... passed directly to estimator for testing purposes

test_estimator2 <- function(estimator, Y, w, ...){
  if (is.null(estimator)){stop("estimator is NULL")}
  if (length(setdiff(c("Y", "w"), formalArgs(estimator))) == 0){
    hasYnw <- TRUE
  } else {
    message("Estimator must have arguments 'Y' and 'w'.")
    hasYnW = FALSE
  }

  if (is.null(w)){
    w <- runif(nrow(Y), 0, 1)
    w <- w / sum(w)
  }

  estargs <- list(...)
  estargs$w <- w
  estargs$Y <- Y
  if ("paramvec" %in% formalArgs(estimator)){
    paramvec = TRUE
    estargs$paramvec <- t_sf2u(starttheta, isfixed)
  } else {
    paramvec = FALSE
  }

  if ("paramvec_start" %in% formalArgs(estimator)){
    paramvec_start = TRUE
  } else {
    paramvec_start = FALSE
  }
  

  estobj <- do.call(estimator, estargs)
  estlocation <- find_paramvec_location(estobj)

  newtheta <- extract_paramvec(estobj)

  if (!isTRUE(length(newtheta) == length(starttheta))){stop("Estimator must return a vector of the same length as the input parameter vector")}
  if (any(abs(newtheta[isfixed] - starttheta[isfixed]) > sqrt(.Machine$double.eps))){
    stop("The fixed elements of theta are altered by estimator.")
  }
  return(list(
    paramvec = paramvec,
    paramvec_start = paramvec_start,
    estlocation = estlocation,
    passes_tests = all(hasYnW)
  ))
}



