test_that("a cppad error creates an R error", {
  # this should error even when NDEBUG is defined because the CppAD Macros that silence the erroring do so before calling the handler
  expect_error(test_Rcpphandler())
})
