test_that("von-Mises Fisher likelihood runs and fits", {
  set.seed(123)
  theta <-  3 * c(1, -1) / sqrt(2)
  sample <- movMF::rmovMF(1000, theta)
  thetatape <- c(1, -1)

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], thetatape, llname = "vMF", pman,
                  fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)
  expect_equal(pForward0(lltape, sample[1, ], theta), sum(sample[1, ]  * theta)) ## very important to check a tape

  smotape<- ptapesmo(sample[1,], thetatape,
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function

  out <- smest(smotape, thetatape, sample, control = list(tol = 1E-10))
  cdabyppi:::expect_lt_v(abs(out$par - theta), 3 * out$SE)
})

test_that("vMF() function works", {
  set.seed(123)
  theta <-  3 * c(1, -1) / sqrt(2)
  sample <- movMF::rmovMF(100, theta)
  out <- vMF(sample)
  cdabyppi:::expect_lt_v(abs(out$par - theta), 3 * out$SE)
})
