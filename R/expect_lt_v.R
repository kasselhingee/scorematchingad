#' Custom expecations for vector comparisons in testthat
#' @examples
#' object = c(1, 3, 5)
#' expected = 3
#' expect_lt_v(c(1,3,5), 3)
expect_lt_v <- function(object, expected){
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

  # check
  if (length(expected) != length(object)){
    if (length(object) %% length(expected) == 0){
      expected <- rep(expected, length(object) / length(expected))
    } else {
      testthat::fail("Length mismatch: %s has length %i but expected has length %i", act$lab, length(object), length(expected))
    }
  }
  comparison <- (object < expected)
  ok = all(comparison)
  if (ok){
    succeed()
    return(invisible(act$val))
  }
  failid <- which(!comparison)
  faildesc <- paste0(sprintf("%s was larger than expected at index c(%s).", act$lab, paste(failid, collapse = ", ")))
  faildetail <- paste(failid, ": ", object[failid], "!<",expected[failid])
  testthat::fail(paste(c(faildesc, faildetail), collapse = "\n"))
}
