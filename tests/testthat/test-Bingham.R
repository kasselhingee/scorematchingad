test_that("taped Bingham log-likelihood gives correct values", {
  p <- 4
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)

  A <- Bingham_theta2Amat(theta)
  expect_equal(Bingham_Amat2theta(A), theta)

  set.seed(123)
  sample <- rBingham(2, A)

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[2,], theta + 1, llname = "Bingham", pman,
                    fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

  u <- t(sample[1, , drop = FALSE])

  expect_equal(pForward0(lltape, u, theta), t(u) %*% A %*% u,
               ignore_attr = TRUE)

  expect_equal(pJacobian(lltape, u, theta), 2 * A %*% u,
               ignore_attr = TRUE)

  # test deriv wrt theta
  llBingham <- function(theta){
    A <- Bingham_theta2Amat(theta)
    return(t(u) %*% A %*% u)
  }
  Rgradt <- numericDeriv(quote(llBingham(theta)), c("theta"))
  lltape_t <- swapDynamic(lltape, theta+1, sample[2, ])
  expect_equal(pJacobian(lltape_t, theta, u), attr(Rgradt, "gradient"),
               tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Bingham_full() optimiser works", {
  p <- 4
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)

  set.seed(123)
  sample <- rBingham(100, A)
  est <- Bingham_full(sample)
  cdabyppi::expect_lt_v(abs(est$sminfo$par - theta), 3 * est$sminfo$SE)
  expect_lt(est$sminfo$sqgradsize, 1E-10)
})
