test_that("Raw taylor approximation works on single sample points", {
  m = sec2_3model(1)

  psimplex <- pmanifold("simplex") #because above ppill_r is for the simplex
  lltape <- ptapell(m$sample, m$theta, llname = "ppi", pman = psimplex, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_equal(pForward0(lltape, m$sample, m$theta), ptaylorapprox(lltape, m$sample, m$sample + 1E-1, m$theta, 100))

  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psimplex, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, m$sample, m$theta)
  expect_equal(pForward0(smoppi, m$theta, m$sample), pForward0(smoppi_u, m$sample, m$theta)) #sanity check swap of dynamic param and independent params
  expect_equal(pForward0(smoppi_u, m$sample, m$theta),
               ptaylorapprox(smoppi_u, m$sample, m$sample + 1E-1, m$theta, 100))

  # test a lower order approximation
  expect_equal(pForward0(smoppi_u, m$sample, m$theta),
               ptaylorapprox(smoppi_u, m$sample, m$sample + 1E-1, m$theta, 10),
               tolerance = 1E-4)
  # it is much faster, but not as accurate
})
