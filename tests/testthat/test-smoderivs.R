test_that("Gradient of smo for ppi wrt u is CLOSE TO CORRECT for interior points, and numerically consistent with cppad smo values", {
  # set.seed(123)
  m <- ppi_egmodel(2)
  acut <- 0.1

  psphere <- manifoldtransform("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, c(0.1,0.1,0.1), m$theta) #don't use a boundary point here!

  testcanntheta <- toPPIcannparam(m$AL + 1, m$bL + 1, m$beta0 + 1)
  testtheta <- toPPIparamvec(m$AL + 1, m$bL + 1, m$beta0 + 1)

  # double check that smo values are equal
  hardcodedsmoval <- estimatorall1_smo(testcanntheta, m$sample[1, , drop = FALSE], acut)
  cppadsmoval <- pForward0(smoppi_u, m$sample[1, ], testtheta)
  stopifnot(abs(hardcodedsmoval -cppadsmoval) < 1E-10)

  # test gradients wrt u
  gradu_cppad <- pJacobian(smoppi_u, m$sample[1, ], testtheta)
  u <- m$sample[1, , drop = FALSE]
  gradu_hardcoded <- numericDeriv(quote(estimatorall1_smo(testcanntheta, u, acut)), c("u"))
  gradu_cppad_numerical1 <- numericDeriv(quote(pForward0(smoppi, testtheta, u)), c("u"))
  gradu_cppad_numerical2 <- numericDeriv(quote(pForward0(smoppi_u, u, testtheta)), c("u"))
  expect_equal(gradu_cppad, attr(gradu_hardcoded,"gradient"), tolerance = 1E-2, ignore_attr = TRUE)
  expect_equal(attr(gradu_cppad_numerical1, "gradient"), attr(gradu_hardcoded,"gradient"),
               tolerance = 1E-2, ignore_attr = TRUE)
  expect_equal(attr(gradu_cppad_numerical1, "gradient"), gradu_cppad, tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Gradient of smo for ppi wrt theta is correct for interior points", {
  set.seed(123)
  m <- ppi_egmodel(2)
  acut <- 0.1

  psphere <- manifoldtransform("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function #this line fails when testing after devtools::load_all()

  testcanntheta <- toPPIcannparam(m$AL + 1, m$bL + 1, m$beta0 + 1)
  testtheta <- toPPIparamvec(m$AL + 1, m$bL + 1, m$beta0 + 1)

  # double check that smo values are equal
  hardcodedsmoval <- estimatorall1_smo(testcanntheta, m$sample[1, , drop = FALSE], acut)
  cppadsmoval <- pForward0(smoppi, testtheta, m$sample[1, ])
  stopifnot(abs(hardcodedsmoval -cppadsmoval) < 1E-10)

  # test gradients wrt theta
  gradt_cppad <- pJacobian(smoppi, testtheta, m$sample[1, ])
  u <- m$sample[1, , drop = FALSE]
  gradt_hardcoded <- numericDeriv(quote(estimatorall1_smo(testcanntheta, u, acut)), c("testcanntheta"))
  gradt_components <- fromPPIparamvec(attr(gradt_hardcoded, "gradient"), m$p)
  gradt_components$beta <- gradt_components$beta * 2 #to account for cannonical exponential form
  gradt_cppad_numerical <- numericDeriv(quote(pForward0(smoppi, testtheta, u)), c("testtheta"))
  expect_equal(gradt_cppad, do.call(toPPIparamvec, gradt_components), tolerance = 1E-5)
  expect_equal(attr(gradt_cppad_numerical, "gradient"), do.call(toPPIparamvec, gradt_components),
               tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Gradient of smo approxcentre for ppi wrt theta is correct", {
  # set.seed(123)
  m <- ppi_egmodel(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- simplex_boundaryshift(m$sample, shiftsize = 1E-15)
  acut <- 0.1

  psphere <- manifoldtransform("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, c(0.1,0.1,0.1), m$theta) #don't use a boundary point here!

  testcanntheta <- toPPIcannparam(m$AL + 1, m$bL + 1, m$beta0 + 1)
  testtheta <- toPPIparamvec(m$AL + 1, m$bL + 1, m$beta0 + 1)

  # double check that smo values are equal
  hardcodedsmoval <- estimatorall1_smo(testcanntheta, m$sample[1, , drop = FALSE], acut)
  cppadsmoval <- pTaylorApprox(smoppi_u, m$sample[1,], acentres[1,], testtheta, 100)
  stopifnot(abs(hardcodedsmoval -cppadsmoval) < 1E-10)

  # test gradients wrt theta
  gradt_cppad <- pJacobian(smoppi, testtheta, m$sample[1, ])
  u <- m$sample[1, , drop = FALSE]
  gradt_hardcoded <- numericDeriv(quote(estimatorall1_smo(testcanntheta, u, acut)), c("testcanntheta"))
  gradt_components <- fromPPIparamvec(attr(gradt_hardcoded, "gradient"), m$p)
  gradt_components$beta <- gradt_components$beta * 2 #to account for cannonical exponential form
  gradt_cppad_numerical_approx <- numericDeriv(quote(pTaylorApprox(smoppi_u, m$sample[1,], acentres[1,], testtheta, 100)),
                                               c("testtheta"))
  expect_false(isTRUE(all.equal(gradt_cppad, do.call(toPPIparamvec, gradt_components), tolerance = 1E-5))) #fails dtheta at boundary is 0 according to CppAD, this is wrong!!
  expect_equal(attr(gradt_cppad_numerical_approx, "gradient"), do.call(toPPIparamvec, gradt_components),
               tolerance = 1E-5, ignore_attr = TRUE)

  # check that gradient is close at the approximation centre
  expect_equal(pJacobian(smoppi, testtheta, acentres[1, ]), do.call(toPPIparamvec, gradt_components),
               tolerance = 1E-5, ignore_attr = TRUE)

  # check that Hessian is close too (and hopefully close to zero)
  hesst_hardcoded <- numDeriv::hessian(function(ctheta) {estimatorall1_smo(ctheta, u, acut)},
                                    testcanntheta)
  expect_equal(pHessian(smoppi, testtheta, acentres[1, ]), hesst_hardcoded,
               ignore_attr = TRUE, tolerance = 1E-1)

  # Hessian matches at acentre:
  hesst_hardcoded_acentre <- numDeriv::hessian(function(ctheta) {estimatorall1_smo(ctheta, acentres[1, , drop = FALSE], acut)},
                                    testcanntheta)
  # convert to usual beta
  betaindx <- (length(m$theta) - m$p + 1):length(m$theta)
  hesst_hardcoded_acentre[betaindx, ] <- hesst_hardcoded_acentre[betaindx, ] * 2
  hesst_hardcoded_acentre[, betaindx] <- hesst_hardcoded_acentre[, betaindx] * 2

  expect_equal(pHessian(smoppi, testtheta, acentres[1, ]), hesst_hardcoded_acentre,
               ignore_attr = TRUE, tolerance = 1E-5)
})

test_that("Gradient of smo approxcentre for ppi wrt u is close", {
  set.seed(123)
  m <- ppi_egmodel(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- simplex_boundaryshift(m$sample, shiftsize = 1E-3)
  acut <- 0.1

  psphere <- manifoldtransform("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, c(0.1,0.1,0.1), m$theta) #don't use a boundary point here!

  testcanntheta <- toPPIcannparam(m$AL + 1, m$bL + 1, m$beta0 + 1)
  testtheta <- toPPIparamvec(m$AL + 1, m$bL + 1, m$beta0 + 1)

  # double check that smo values are equal
  hardcodedsmoval <- estimatorall1_smo(testcanntheta, m$sample[1, , drop = FALSE], acut)
  cppadsmoval <- pTaylorApprox(smoppi_u, m$sample[1,], acentres[1,], testtheta, 100)
  stopifnot(abs(hardcodedsmoval -cppadsmoval) < 1E-10)

  # test gradients wrt u
  gradu_cppad <- pJacobian(smoppi_u, m$sample[1, ], testtheta)
  u <- m$sample[1, , drop = FALSE]
  gradu_hardcoded <- numericDeriv(quote(estimatorall1_smo(testcanntheta, u, acut)), c("u"))

  # check that gradient at boundary is close to the gradient at the approximation centre
  expect_equal(pJacobian(smoppi_u, acentres[1, ], testtheta), attr(gradu_hardcoded, "gradient"),
               tolerance = 1E-1, ignore_attr = TRUE)


  # check that gradient matches at approximation centre
  u <- acentres[1, , drop = FALSE]
  gradu_hardcoded_acentre <- numericDeriv(quote(estimatorall1_smo(testcanntheta, u, acut)), c("u"))
  expect_equal(pJacobian(smoppi_u, acentres[1, ], testtheta), attr(gradu_hardcoded_acentre, "gradient"),
               tolerance = 1E-3, ignore_attr = TRUE)

  # check that Hessian matches at approximation centre
  hessu_hardcoded_acentre <- numDeriv::hessian(function(u) estimatorall1_smo(testcanntheta, u, acut),
                                            acentres[1, , drop = FALSE])
  expect_equal(pJacobian(smoppi_u, acentres[1, ], testtheta), attr(gradu_hardcoded_acentre, "gradient"),
               tolerance = 1E-3, ignore_attr = TRUE)
})

