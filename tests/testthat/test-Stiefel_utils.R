test_that("The commutation matrix obtains vec(t(A))", {
  A <- matrix(runif(5*3), 5, 3)
  expect_equal(drop(commutation_mat(5, 3)  %*% vec(A)), vec(t(A)))
})