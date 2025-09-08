test_that("The commutation matrix obtains vec(t(A))", {
  A <- matrix(runif(5*3), 5, 3)
  expect_equal(drop(commutation_mat(5, 3)  %*% vec(A)), vec(t(A)))
})

test_that("Stiefel_proj() and Stiefel_projmat() get to same result", {
  A <- rstiefel::rustiefel(5, 3)
  Z <- rstiefel::rustiefel(5, 3)
  expect_equal(Stiefel_proj(A, A), 0*A)
  expect_equal(invvec(Stiefel_projmat(A) %*% vec(Z), 5), Stiefel_proj(Z, A))
})

test_that("Stiefel_projmat_d() matches iterative estimate of derivative", {
  A <- rstiefel::rustiefel(5, 3)
  i <- sample.int(5, 1)
  j <- sample.int(3, 1)
  updateA <- function(A, x){A[i, j] <- x; return(A)} #update element i,j with x
  x <- A[i, j]
  expect_equal(attr(stats::numericDeriv(quote(Stiefel_projmat(updateA(A, x))), "x"), "gradient")[,,1],
               Stiefel_projmat_d(A, i, j),
               tolerance = 1E-7)
})
