# Tests of the ppi() structure

expect_ppi_str <- function(object, p){
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

  #check
  correctfirstnames <- isTRUE(setequal(names(object), c("est", "SE", "info")))
  if (!correctfirstnames){
    testthat::fail("Incorrect slot names")
  }

  correctest <- correcteststr(object$est, p)
  if (object$SE != "Not calculated."){
    correctSE <- correcteststr(object$SE, p)
  } else {correctSE = list(SEnotcalculated = TRUE)}

  correctmethodstr <- isTRUE((class(object$info$method) == "character") & (nchar(object$info$method) > 0))
  passsummary <- list(
    slotnames = correctfirstnames,
    est_structure = correctest,
    SE_structure = correctSE,
    info_method_exists = correctmethodstr
  )
  if (all(unlist(passsummary))){
    succeed()
    return(invisible(act$val))
  }
  else {
    fails <- names(which(!unlist(passsummary)))
    fails <- gsub(".", "$", fails, fixed = TRUE)
    testthat::fail(
      sprintf("The following conditions fail: %s", paste(fails, collapse = ", ")))
  }
}

correcteststr <- function(est, p){
  correctestnames <- setequal(names(est), c("paramvec", "ALs", "bL", "beta"))
  correctestdims <- list()
  if (correctestnames){
    correctestdims$paramvec <- (length(est$paramvec) == cdabyppi:::ppithetalength(p))
    correctestdims$ALs <- isTRUE(all(dim(est$ALs) == c(p-1, p-1)))
    correctestdims$bL <- (length(est$bL) == p-1)
    correctestdims$beta <- (length(est$beta) == p)
  }
  return(c(Names = correctestnames, correctestdims))
}
