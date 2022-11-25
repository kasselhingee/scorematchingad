test_that("ppi tape values do not effect ll values", {
  model1 <- ppi_egmodel(1)
  u1 <-  c(0.001, 0.011, 1 - 0.01 - 0.011)
  model0 <- lapply(model1, function(x) x * 0)
  u0 <- rep(0, 3)
  fixedtheta = rep(FALSE, length(model1$theta))

  ueval <- matrix(c(0.4, 0.011, 1 - 0.4 - 0.011), nrow = 1)
  thetaeval <- model1$theta + 1

  psphere <- pmanifold("sphere")
  pppi1 <- ptapell(u1, model1$theta, "ppi", psphere, fixedtheta = fixedtheta, verbose = FALSE)
  pppi2 <- ptapell(u0, model0$theta, "ppi", psphere, fixedtheta = fixedtheta, verbose = FALSE)

  expect_equal(pForward0(pppi1, ueval, thetaeval), pForward0(pppi2, ueval, thetaeval))
  expect_equal(pJacobian(pppi1, ueval, thetaeval), pJacobian(pppi2, ueval, thetaeval))
  expect_equal(pHessian(pppi1, ueval, thetaeval), pHessian(pppi2, ueval, thetaeval))
})


test_that("ppi and dirichlet smo value match when AL and bL is zero and p = 3", {
  beta = c(-0.3, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  utabl <- rppi(10,beta=beta,AL=ALs,bL=bL,maxden=4)

  acut = 0.1
  dirtapes <- buildsmotape("sphere", "dirichlet",
                           rep(0.1, p), rep(0.1, p) * NA,
                           weightname = "minsq", acut = acut)
  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, p), theta * NA,
                           weightname = "minsq", acut = acut)

  ppival <- pForward0(ppitapes$smotape, theta, utabl[2, ])
  dirval <- pForward0(dirtapes$smotape, beta, utabl[2, ])
  expect_equal(ppival, dirval)
})


test_that("cppad ppi estimate works when AL and bL is zero and p = 4", {
  beta = c(-0.3, -0.2, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  set.seed(1234)
  utabl <- rppi(100,beta=beta,AL=ALs,bL=bL,maxden=4)

  acut = 0.1
  dirtapes <- buildsmotape("sphere", "dirichlet",
                           rep(0.1, p), rep(0.1, p) * NA,
                           weightname = "minsq", acut = acut)
  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, p), theta * NA,
                           weightname = "minsq", acut = acut)

  # it looks like the taped function above is not altering bL or beta
  # potentially the ordering of the theta values is wrong??
  out <- cppadest(ppitapes$smotape, theta * 0 + 1, utabl, list(tol = 1E-10))
  stopifnot(out$convergence == 0)

  cppadest <- fromPPIparamvec(out$par, p)

  expect_equal(pForward0(dirtapes$lltape, utabl[2, ], beta), pForward0(ppitapes$lltape, utabl[2, ], theta))
  expect_equal(pJacobian(dirtapes$lltape, utabl[2, ], beta), pJacobian(ppitapes$lltape, utabl[2, ], theta))
  expect_equal(pForward0(dirtapes$smotape, beta, utabl[2, ]), pForward0(ppitapes$smotape, theta, utabl[2, ]))

  hardcodedestimate <- dir_sqrt_minimah(utabl, acut)

  # there is a difference in hardcoded estimates because the hardcoded estimate smobj value is poorer:
  expect_lt(out$value,
             smobj(ppitapes$smotape, c(rep(0, length(theta) - p), hardcodedestimate), utabl) + 1E-3 * abs(out$value)) #larger tolerance some apriori information used
  expect_lt(out$sqgradsize,
               sum(smobjgrad(ppitapes$smotape, c(rep(0, length(theta) - p), hardcodedestimate), utabl)^2))

  cppadSE <- fromPPIparamvec(out$SE, p)
  expect_absdiff_lte_v(cppadest$beta[1:3], hardcodedestimate[1:3], cppadSE$beta[1:3])
  expect_lt(abs(cppadest$beta - hardcodedestimate)[4] / cppadSE$beta[4], 3) #largest beta is hard to estimate well

  expect_lt_v(abs(out$par - theta) / out$SE, 3)#assuming normally distributed with SE given by SE above
})

test_that("ppi with minsq weights match estimator1 with fixed beta for ppi_egmodel", {
  set.seed(123)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(model$sample, paramvec = ppi_paramvec(betaL = model$beta0[1:2], betap = model$beta0[3]),
            method = "cppad",
                   bdrythreshold = 1E-10,
            trans = "sqrt", divweight = "minsq", acut = acut)

  hardcodedestimate <- estimator1(model$sample, acut, incb = TRUE, beta = model$beta0)

  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, model$p), ppi_paramvec(model$p, beta = model$beta0),
                           weightname = "minsq", acut = acut)

  expect_equal(out$info$smval,
               smobj(ppitapes$smotape, hardcodedestimate$est$paramvec[1:(length(model$theta)-3)], model$sample)) #failing now?
  expect_equal(out$info$sqgradsize,
               sum(smobjgrad(ppitapes$smotape, hardcodedestimate$est$paramvec[1:(length(model$theta)-3)], model$sample)^2))

  expect_absdiff_lte_v(out$est$paramvec, hardcodedestimate$est$paramvec, 0.01 * out$SE$paramvec) #proxy for optimisation flatness
  expect_absdiff_lte_v(out$est$paramvec, model$theta, 2 * out$SE$paramvec)
})

test_that("ppi with prodsq weights match estimator1 with fixed beta for ppi_egmodel", {
  set.seed(123)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(model$sample, ppi_paramvec(betaL = model$beta0[1:2], betap = model$beta0[3]),
             method = "cppad",
                   trans = "sqrt", divweight = "prodsq", acut = acut, 
                   control = list(tol = 1E-12))

  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, model$p), ppi_paramvec(model$p, betaL = model$beta0[1:2], betap = model$beta0[3]),
                           weightname = "prodsq", acut = acut)

  hardcodedestimate <- estimator2(model$sample, acut, incb = TRUE, beta0 = model$beta0)

  expect_equal(out$info$smval,
               smobj(ppitapes$smotape, hardcodedestimate$estimator2, model$sample))
  expect_equal(out$info$sqgradsize,
               sum(smobjgrad(ppitapes$smotape, hardcodedestimate$estimator2, model$sample)^2))

  expect_absdiff_lte_v(out$est$paramvec, c(hardcodedestimate$estimator2, model$beta0), 0.01 * out$SE$paramvec) #proxy for optimisation flatness

  expect_absdiff_lte_v(out$est$paramvec, model$theta, 3 * out$SE$paramvec)
})

test_that("ppi with minsq weights match estimatorall1 for p = 4, mostly zero params", {
  acut = 0.1
  #sample size
  n=10000

  #parameters for the PPI model
  beta = c(-0.3, -0.2, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  ALs[1, 2] <- 1
  ALs[2, 1] <- 1
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta <- toPPIparamvec(ALs, bL, beta)

  set.seed(13418)
  utabl <- rppi(n,beta=beta,AL=ALs,bL=bL,maxden=4)
  u <- utabl[2, ]

  out <- ppi(utabl,
             method = "cppad",
                   trans = "sqrt", divweight = "minsq", acut = acut,
                   control = list(tol = 1E-10))

  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, p), ppi_paramvec(p),
                           weightname = "minsq", acut = acut)

  hardcodedestimate <- estimatorall1(utabl, acut)

  expect_equal(out$info$smval,
               smobj(ppitapes$smotape, hardcodedestimate$estimator1, utabl),
               tolerance = 1E-2)
  expect_equal(out$info$sqgradsize,
               sum(smobjgrad(ppitapes$smotape, hardcodedestimate$estimator1, utabl)^2))

  expect_absdiff_lte_v(out$est$paramvec, hardcodedestimate$estimator1, out$SE$paramvec) #proxy for optimisation flatness

  expect_absdiff_lte_v(out$est$paramvec, theta, 3 * out$SE$paramvec)
})

test_that("ppi with minsq weights match estimatorall1 for ppi_egmodel", {
  set.seed(111)
  model <- ppi_egmodel(100, maxden = 4)

  acut = 0.1

  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(model$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:length(model$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- cppadest(smoppi, model$theta * 0 + 1, model$sample, control = list(tol = 1E-12))

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  hardcodedestimate <- estimatorall1(model$sample, acut)

  expect_lt(out$value,
            smobj(smoppi, hardcodedestimate$estimator1, model$sample) + 1E-5 * abs(out$value))

  expect_lt_v(abs(out$par - hardcodedestimate$estimator1) / out$SE, 1) #proxy for optimisation flatness

  expect_lt_v(abs(out$par - model$theta) / out$SE, 3)
})

test_that("ppi with minsq weights match estimatorall1 for ppi_egmodel, fixed final beta", {
  set.seed(123)
  model <- ppi_egmodel(100, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- cppadest(smoppi, model$theta[-length(model$theta)] * 0, model$sample,
               control = list(tol = 1E-10))

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  hardcodedestimate <- estimatorall1(model$sample, acut, betap = model$beta0[model$p])

  expect_lte(out$value,
             smobj(smoppi, hardcodedestimate$estimator1, model$sample) + 1E-3 * abs(out$value))

  expect_lt_v(abs(out$par - hardcodedestimate$estimator1) / out$SE, 0.01) #proxy for optimisation flatness

  expect_lt_v(abs(out$par - model$theta[-length(model$theta)]) / out$SE, 3)
})

test_that("ppi with minsq weights match estimatorall1 for ppi_egmodel, fixed final beta, large n", {
  set.seed(123)
  model <- ppi_egmodel(100000, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- cppadest(smoppi, model$theta[-length(model$theta)] * 0, model$sample,
               control = list(tol = 1E-10, trace = 0)) #very slow at 10000

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  hardcodedestimate <- estimatorall1(model$sample, acut, betap = model$beta0[model$p])

  expect_lte(out$value,
             smobj(smoppi, hardcodedestimate$estimator1, model$sample) + 1E-3 * abs(out$value))

  expect_lt_v(abs(out$par - hardcodedestimate$estimator1) / out$SE, 0.01) #proxy for optimisation flatness

  expect_lt_v(abs(out$par - model$theta[-length(model$theta)]) / out$SE, 3)
})

test_that("ppi with minsq weights performs well on simplex, fixed final beta", {
  set.seed(1234)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  psimplex <- pmanifold("simplex")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psimplex, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = psimplex, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- cppadest(smoppi, model$theta[-length(model$theta)] * 0, model$sample,
               control = list(tol = 1E-12))

  expect_lt_v(abs(out$par - model$theta[-length(model$theta)]) / out$SE, 3)
})

test_that("ppi via cppad matches Score1 for p=5, particularly the order of the off diagonals in ALs", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  prop <- rppi(1000, beta=beta, AL=ALs, bL=bL, maxden=35)

  acut = 0.1
  est_hardcoded <- estimator1(prop, acut, incb = 0, beta = beta)

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL, beta = beta),
                   method = "cppad",
                         trans = "sqrt", acut = acut, divweight = "minsq",
                         control = list(tol = 1E-13))
  expect_equal(est_cppad$est$paramvec, est_hardcoded$est$paramvec, tolerance = 1E-1,
               ignore_attr = TRUE)

  #also it makes sense that the smo and gradient are v low at the hardcoded estimate
  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, p), ppi_paramvec(p, bL = bL, beta = beta),
                           weightname = "minsq", acut = acut)
  expect_lt(sum(smobjgrad(ppitapes$smotape, fromsmatrix(est_hardcoded$est$ALs), prop)^2), 1E-20)
  expect_equal(smobj(ppitapes$smotape, fromsmatrix(est_hardcoded$est$ALs), prop), est_cppad$info$smval, tolerance = 1E-1)

  # check that rearrangement has large gradient
  expect_gt(sum(smobjgrad(ppitapes$smotape, fromsmatrix(est_hardcoded$est$ALs)[c(1:6, 8, 7, 9, 10)], prop)^2), 1E-2)

})
