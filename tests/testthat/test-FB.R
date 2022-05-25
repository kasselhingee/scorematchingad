test_that("Fisher-Bingham likelihood runs and fits", {
  set.seed(321)
  p <- 3
  A <- matrix(NA, ncol = p, nrow = p)
  A[upper.tri(A)] <- runif(sum(upper.tri(A)))
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- c(runif(p-1), NA)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  m <- runif(p)
  m <- m / sqrt(sum(m^2))
  k <- 2
  theta <- FB_mats2theta(k, m, A)
  expect_equal(FB_theta2mats(theta),
    list(k = k,
         m = m,
         A = A))


  set.seed(123)
  sample <- Directional::rfb(100, k, m, A)

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], seq.int(1, length.out = length(theta)), llname = "FB", pman,
                    fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

  expect_equal(pForward0(lltape, sample[1, ], theta), log(qdFB(sample[1, ], k, m, A)),
               ignore_attr = TRUE)

  expect_equal(pForward0(lltape, sample[1, ], theta), sum(sample[1, ]  * theta)) ## very important to check a tape

  smotape<- ptapesmo(sample[1,], thetatape,
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function

  out <- smest(smotape, thetatape, sample, control = list(tol = 1E-10))
  cdabyppi:::expect_lt_v(abs(out$par - theta), 3 * out$SE)
})
