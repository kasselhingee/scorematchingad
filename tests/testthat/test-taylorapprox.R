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

test_that("Approx taylor with u on boundary generates correct values (excluding gradient) for sphere for ppi", {
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
  pi <- toPPIparamvec(m$ALs, m$bL, beta = 1 + 2 * m$beta0)
  expect_equal(approxsmoval, estimatorall1_smo(pi, m$sample[1,, drop = FALSE], 0.1))
})

test_that("Taylor Approx of Grad SMO gets correct value on interior of simplex", {
  set.seed(123)
  m <- sec2_3model(2)
  # m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- approxcentre(m$sample, shiftsize = 1E-15)

  psphere <- pmanifold("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_true(!is.nan(pTaylorApprox(lltape, m$sample[1,], acentres[1,], m$theta, 100)))

  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function
  smoppi_u <- swapDynamic(smoppi, c(0.1,0.1,0.1), m$theta) #don't use a boundary point here!

  Jsmoppi_u <- pTapeJacobianSwap(smoppi, m$theta, c(0.1,0.1,0.1))
  expect_equal(pForward0(Jsmoppi_u, c(0.1,0.1,0.1), m$theta), pJacobian(smoppi, m$theta, c(0.1,0.1,0.1)))


  approxgrad <- pTaylorApprox(Jsmoppi_u, m$sample[1, ], acentres[1, ], m$theta, 10)
  cppadgrad <- pJacobian(smoppi, m$theta, m$sample[1, ])
  expect_equal(approxgrad, cppadgrad)
})

test_that("Taylor Approx of Grad SMO gets correct value on boundary of simplex", {
  set.seed(123)
  m <- sec2_3model(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- approxcentre(m$sample, shiftsize = 1E-15)
  acut = 0.1

  psphere <- pmanifold("sphere") #because above ppill_r is for the simplex
  lltape <- ptapell(c(0.1,0.1,0.1), m$theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  expect_true(!is.nan(pTaylorApprox(lltape, m$sample[1,], acentres[1,], m$theta, 100)))

  smoppi <- ptapesmo(c(0.1,0.1,0.1), m$theta, pll = lltape, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  Jsmoppi_u <- pTapeJacobianSwap(smoppi, m$theta, c(0.1,0.1,0.1))

  # actual testing
  testcanntheta <- toPPIcannparam(m$ALs + 1, m$bL + 1, m$beta0 + 1)
  testtheta <- toPPIparamvec(m$ALs + 1, m$bL + 1, m$beta0 + 1)
  approxgrad <- pTaylorApprox(Jsmoppi_u, m$sample[1, ], acentres[1, ], testtheta, 10)
  u <- m$sample[1, , drop = FALSE]
  gradt_direct <- numericDeriv(quote(estimatorall1_smo(testcanntheta, u, acut)), c("testcanntheta"))
  gradt_components <- fromPPIparamvec(attr(gradt_direct, "gradient"), m$p)
  gradt_components$beta <- gradt_components$beta * 2 #to account for cannonical exponential form
  expect_equal(approxgrad, do.call(toPPIparamvec, gradt_components),
               tolerance = 1E-5, ignore_attr = TRUE)
})


test_that("Test smest_simplex against direct", {
  set.seed(123)
  m <- sec2_3model(100)
  #add some zeroes
  pushtozero <- function(x){
    if (min(x) > 1E-3){return(x)}
    whichmin <- which.min(x)
    x[whichmin] <- 0
    x <- x / sum(x) #normalise
    return(x)
  }
  newsample <- t(apply(m$sample, MARGIN = 1, pushtozero))
  mean(apply(newsample, 1, min) == 0) #28% have a zero

  acut = 0.1
  direct <- estimatorall1(newsample, acut = acut, betap = m$beta0[3])

  #copied from ppi_cppad()
  theta <- cdabyppi:::ppi_cppad_thetaprocessor(3, betap = m$beta0[3])
  fixedtheta <- !is.na(theta)

  # prepare tapes
  pman <- pmanifold("sphere")
  thetatape <- theta  #must pass the fixed values as the taped value
  thetatape[!fixedtheta] <- 0.73 # any number will do!

  lltape <- ptapell(rep(0.1, m$p), thetatape, llname = "ppi", pman = pman, fixedtheta = fixedtheta, verbose = FALSE)
  smoppi <- ptapesmo(c(0.1,0.1,0.1), thetatape[is.na(theta)], pll = lltape, pman = pman, "minsq", acut = 0.1, verbose = FALSE) #tape of the score function

  est_cppad <- smest_simplex(smoppi, thetatape[is.na(theta)] * 0 - 0.1, newsample,
                             control = list(tol = 1E-10), 1E-5) #takes a long time partially because no gradient supplied
  smestSE_simplex(smoppi, est_cppad$par, newsample, shiftsize = 1E-5)
  cdabyppi:::expect_lt_v(abs(est_cppad$par - direct$estimator1), 0.1*abs(direct$estimator1))
})
