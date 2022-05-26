test_that("Rivest likelihood runs and matches R code", {
  set.seed(321)
  p <- 3
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  A_es <- eigen(A)
  idx <- 1
  evalorder <- order(abs(A_es$values), decreasing = TRUE)
  m <- A_es$vectors[, evalorder[idx]]
  k <- 2
  sample <- rFB(1000, k, m, A)

  theta <- c(Bingham_Amat2theta(A), k, idx)

  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], seq.int(1, length.out = length(theta)), llname = "Rivest", pman,
                    fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

  expect_equal(pForward0(lltape, sample[1, ], theta), log(qdRivest(sample[1, ], k, A)),
               ignore_attr = TRUE)## very important to check a tape
  #deriv wrt u
  u <- sample[1, ]
  expect_equal(pJacobian(lltape, sample[1, ], theta),
               attr(numericDeriv(quote(log(qdFB(u, k, m, A))), "u"), "gradient"),
               ignore_attr = TRUE, tolerance = 1E-5)
  #deriv wrt theta
  lltape_theta <- swapDynamic(lltape, theta, sample[1, ])
  qdFB_theta <- function(u, theta){
    mats <- FB_theta2mats(theta)
    qdFB(u, mats$k, mats$m, mats$A)
  }
  expect_equal(pJacobian(lltape_theta, theta, u),
               attr(numericDeriv(quote(log(qdFB_theta(u, theta))), "theta"), "gradient"),
               ignore_attr = TRUE, tolerance = 1E-5)
})
