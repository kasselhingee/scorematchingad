test_that("Raw taylor approximation works on single sample points", {
  set.seed(123)
  m = ppi_egmodel(1)
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
  m = ppi_egmodel(1)
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

test_that("Approx taylor with u on boundary generates correct values (excluding gradient) for sphere for ppi", {
  set.seed(123)
  m <- ppi_egmodel(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- simplex_boundaryshift(m$sample, shiftsize = 1E-15)

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
  m <- ppi_egmodel(2)
  # m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- simplex_boundaryshift(m$sample, shiftsize = 1E-15)

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
  m <- ppi_egmodel(2)
  m$sample[1, ] <- c(0, 0.08, 0.92) #make first measurement on boundary
  acentres <- simplex_boundaryshift(m$sample, shiftsize = 1E-15)
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


test_that("Test ppi() against direct when there are boundary points, with and without weights", {
  set.seed(123)
  m <- ppi_egmodel(100)
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
  direct <- ppi(newsample, 
                paramvec = ppi_paramvec(p = ncol(newsample), betap = tail(m$beta0, 1)),
                trans = "sqrt", bdryweight = "minsq", 
                acut = acut)

  est <- ppi(newsample, ppi_paramvec(p = 3, betap = m$beta0[3]), trans = "sqrt", bdryweight = "minsq", acut = acut, method = "cppad",
                   control = list(tol = 1E-10))
  expect_equal(est$est$paramvec, direct$est$paramvec, tolerance = 1E-4)

  # with weights!!
  w <- runif(nrow(newsample))

  direct <- ppi(newsample, 
                paramvec = ppi_paramvec(p = ncol(newsample), betap = tail(m$beta0, 1)),
                trans = "sqrt", bdryweight = "minsq", 
                acut = acut,
                w = w)

  est <- ppi(newsample, ppi_paramvec(p = 3, betap = m$beta0[3]), trans = "sqrt", bdryweight = "minsq", acut = acut, method = "cppad",
                   w = w,
                   control = list(tol = 1E-10))
  expect_equal(est$est$paramvec, direct$est$paramvec, tolerance = 1E-4)

})

test_that("Taylor approx of cppadSE gives suitable SE for estimates", {
  set.seed(123)
  m <- ppi_egmodel(100)
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
  est <- ppi(newsample, ppi_paramvec(p=3, betap = m$beta0[3]), trans = "sqrt", bdryweight = "minsq", acut = acut, method = "cppad",
                   control = list(tol = 1E-10))
  expect_absdiff_lte_v(est$est$paramvec, m$theta, 3 * est$SE$paramvec)
  # note that the SE is a bit hard to test on here because the data has been truncated
})

test_that("Taylor approx of ppi() SE matches on the interior", {
  set.seed(123)
  m <- ppi_egmodel(100)

  acut = 0.1
  direct <- estimatorall1(m$sample, acut = acut, betap = m$beta0[3])

  est_default <- ppi(m$sample, ppi_paramvec(p=3, betap = m$beta0[3]), trans = "sqrt", bdryweight = "minsq", acut = acut, method = "cppad",
                   control = list(tol = 1E-10))

  est_interior <- ppi(m$sample, ppi_paramvec(p=3, betap = m$beta0[3]), trans = "sqrt", bdryweight = "minsq", acut = acut, method = "cppad",
                            bdrythreshold = -1,
                            control = list(tol = 1E-10))

  expect_equal(est_default$est$theta, est_interior$est$theta, tolerance = 1E-3)
  expect_equal(est_default$SE$theta, est_interior$SE$theta, tolerance = 1E-3)
})

test_that("ppi() operates when minimal points in the interior", {
  set.seed(123)
  m <- ppi_egmodel(100)
  #add some zeroes
  pushtozero <- function(x){
    if (min(x) > 10){return(x)}
    whichmin <- which.min(x)
    x[whichmin] <- 0
    x <- x / sum(x) #normalise
    return(x)
  }
  newsample <- t(apply(m$sample, MARGIN = 1, pushtozero))
  newsample <- rbind(newsample, m$sample[1:3, , drop = FALSE])

  acut = 0.1
  direct <- estimator1(newsample, acut = acut, incb = 1, beta = m$beta0)

  est <- ppi(newsample, ppi_paramvec(betaL = m$beta0[1:2], betap = m$beta0[3]), trans = "sqrt", bdryweight = "minsq", acut = acut, method = "cppad",
                            control = list(tol = 1E-10))
  expect_absdiff_lte_v(est$est$paramvec, direct$est$paramvec, 1E-1 * abs(direct$est$paramvec))
})

test_that("Taylor approx of matches estimator1SE with data on the boundary", {
  set.seed(123)
  m <- ppi_egmodel(100)
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
  direct <- estimator1(newsample, acut = acut, incb = 1, beta = m$beta0, computeSE = TRUE)
  directSE <- direct$SE$paramvec

  intheta <- ppi_paramvec(3, betaL = m$beta0[1:2], betap = m$beta0[3])

  # prepare tapes
  tapes <- buildsmotape("sphere", "ppi",
                        rep(0.1, m$p), intheta,
                        weightname = "minsq",
                        acut = acut)

  #prepare data
  datasplit <- simplex_boundarysplit(newsample, bdrythreshold = 1E-15, shiftsize = 1E-15)

  #extra tapes
  Jsmofun_u <- pTapeJacobianSwap(tapes$smotape, m$theta, datasplit$interior[1, ])
  Hsmofun_u <- pTapeHessianSwap(tapes$smotape, m$theta, datasplit$interior[1, ])

  #comparisons isolated from the Rcgmin optimiser
  SE <- cppadSE(
    tapes$smotape, theta = direct$est$paramvec[1:(length(intheta) - m$p)], datasplit$interior,
    Jsmofun_u = Jsmofun_u,
    Hsmofun_u = Hsmofun_u,
    uboundary = datasplit$uboundary, boundaryapprox = datasplit$boundaryapprox,
    approxorder = 100)
  expect_equal(SE, directSE[1:(length(intheta) - m$p)], tolerance = 1E-2)
})
