test_that("Raw taylor approximation works on single sample points", {
  set.seed(123)
  m = sec2_3model(1)
  middleofsimplex <- rep(1/m$p, m$p)
  shiftdir <- middleofsimplex - m$sample
  shiftdir <- shiftdir / sqrt(sum(shiftdir^2)) #make a unit vector
  shiftsize = 1E-2
  approxcentre <- m$sample + shiftsize * shiftdir
  approxcentre <- approxcentre / sum(approxcentre) # normalise to get back onto simplex

  psimplex <- pmanifold("simplex") #because above ppill_r is for the simplex
  lltape <- ptapell(m$sample, m$theta, llname = "ppi", pman = psimplex, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_equal(pForward0(lltape, m$sample, m$theta), pTaylorApprox(lltape, m$sample, approxcentre, m$theta, 100))

  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psimplex, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, m$sample, m$theta)
  expect_equal(pForward0(smoppi, m$theta, m$sample), pForward0(smoppi_u, m$sample, m$theta)) #sanity check swap of dynamic param and independent params
  expect_equal(pForward0(smoppi_u, m$sample, m$theta),
               pTaylorApprox(smoppi_u, m$sample, approxcentre, m$theta, 100))

  # test a lower order approximation
  expect_equal(pForward0(smoppi_u, m$sample, m$theta),
               pTaylorApprox(smoppi_u, m$sample, approxcentre, m$theta, 10),
               tolerance = 1E-4)
  # it is much faster, but not as accurate
})

test_that("Auto approx centre with taylor approximation works on single sample points", {
  set.seed(123)
  m = sec2_3model(1)
  middleofsimplex <- rep(1/m$p, m$p)
  approxcentre <- ((m$sample * (99/100) + middleofsimplex * (1/100)))

  psimplex <- pmanifold("simplex") #because above ppill_r is for the simplex
  lltape <- ptapell(m$sample, m$theta, llname = "ppi", pman = psimplex, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_equal(pForward0(lltape, m$sample, m$theta), pTaylorApprox(lltape, m$sample, approxcentre, m$theta, 100))

  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psimplex, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, m$sample, m$theta)
  expect_equal(pForward0(smoppi, m$theta, m$sample), pForward0(smoppi_u, m$sample, m$theta)) #sanity check swap of dynamic param and independent params
  expect_equal(pForward0(smoppi_u, m$sample, m$theta),
               pTaylorApprox(smoppi_u, m$sample, approxcentre, m$theta, 100))

  # test a lower order approximation
  expect_equal(pForward0(smoppi_u, m$sample, m$theta),
               pTaylorApprox(smoppi_u, m$sample, approxcentre, m$theta, 10),
               tolerance = 1E-4)
  # it is much faster, but not as accurate
})

test_that("approxcentre_simplex gives results equal to shiftsize", {
  m <- sec2_3model(10)
  centres <- approxcentre(m$sample, shiftsize = 1E-5)
  expect_equal(sqrt(rowSums((centres - m$sample)^2)), rep(1E-5, nrow(m$sample)))
})


test_that("Approx taylor with u on boundary generates non-NAN values for sphere for ppi", {
  set.seed(123)
  m <- sec2_3model(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- approxcentre(m$sample, shiftsize = 1E-15)

  psphere <- pmanifold("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_true(!is.nan(pTaylorApprox(lltape, m$sample[1,], acentres[1,], m$theta, 100)))

  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, c(0.1,0.1,0.1), m$theta) #don't use a boundary point here!

  expect_true(is.nan(pForward0(smoppi_u, m$sample[1,], m$theta)))
  approxsmoval <- pTaylorApprox(smoppi_u, m$sample[1, ], acentres[1,], m$theta, 100)
  expect_true(!is.nan(approxsmoval))
  expect_equal(approxsmoval, pForward0(smoppi_u, acentres[1,], m$theta), tolerance = 1E-2)

  #close to the boundary the values gradient should be close flat
  cdabyppi::expect_lt_v(abs(pJacobian(smoppi_u, acentres[1,], m$theta)), 1E-1)

  #close to the boundary the values of Hessian should be close flat
  cdabyppi::expect_lt_v(abs(pHessian(smoppi_u, acentres[1,], m$theta)), 1E-1) #currently there are some HUGE values

  # previously the Jacobian and Hessian gave real answers on the boundary. What has happened?

})

test_that("Approx taylor with u on boundary generates non-NAN values for sphere for dirichlet", {
  set.seed(123)
  m <- sec2_3model(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- approxcentre(m$sample, shiftsize = 1E-15)

  psphere <- pmanifold("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "dirichlet", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_true(!is.nan(pTaylorApprox(lltape, m$sample[1,], acentres[1,], m$theta, 100)))

  smo <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function
  smo_u <- swapDynamic(smo, c(0.1,0.1,0.1), m$theta) #don't use a boundary point here!

  expect_true(is.nan(pForward0(smo_u, m$sample[1,], m$theta)))
  approxsmoval <- pTaylorApprox(smo_u, m$sample[1, ], acentres[1,], m$theta, 100)
  expect_true(!is.nan(approxsmoval))
  expect_equal(approxsmoval, pForward0(smo_u, acentres[1,], m$theta), tolerance = 1E-2)

  #close to the boundary the values gradient should be close flat
  cdabyppi::expect_lt_v(abs(pJacobian(smo_u, acentres[1,], m$theta)), 1E-1)

  #close to the boundary the values of Hessian should be close flat
  cdabyppi::expect_lt_v(abs(pHessian(smo_u, acentres[1,], m$theta)), 1E-1) #currently there are some HUGE values

  # previously the Jacobian and Hessian gave real answers of zero on the boundary. What has happened?
})
