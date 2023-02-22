test_that("ppi tape values do not effect ll values", {
  model1 <- ppi_egmodel(1)
  u1 <-  c(0.001, 0.011, 1 - 0.01 - 0.011)
  model0 <- lapply(model1, function(x) x * 0)
  u0 <- rep(0, 3)
  fixedtheta = rep(FALSE, length(model1$theta))

  ueval <- matrix(c(0.4, 0.011, 1 - 0.4 - 0.011), nrow = 1)
  thetaeval <- model1$theta + 1

  psphere <- manifoldtransform("sphere")
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
  theta = ppi_paramvec(AL=ALs, bL=bL, beta=beta)

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
  theta = ppi_paramvec(AL=ALs, bL=bL, beta=beta)

  set.seed(1234)
  utabl <- rppi(100,beta=beta,AL=ALs,bL=bL,maxden=4)

  acut = 0.1
  dirtapes <- buildsmotape("sphere", "dirichlet",
                           rep(0.1, p), rep(0.1, p) * NA,
                           weightname = "minsq", acut = acut)
  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, p), ppi_paramvec(AL = 0, bL=0, p = p),
                           weightname = "minsq", acut = acut)

  # it looks like the taped function above is not altering bL or beta
  # potentially the ordering of the theta values is wrong??
  out <- cppad_closed(ppitapes$smotape, Y = utabl)

  expect_equal(pForward0(dirtapes$lltape, utabl[2, ], beta), pForward0(ppitapes$lltape, utabl[2, ], beta))
  expect_equal(pJacobian(dirtapes$lltape, utabl[2, ], beta), pJacobian(ppitapes$lltape, utabl[2, ], beta))
  expect_equal(pForward0(dirtapes$smotape, beta, utabl[2, ]), pForward0(ppitapes$smotape, beta, utabl[2, ]))

  hardcodedestimate <- dir_sqrt_minimah(utabl, acut)

  expect_equal(out$est, hardcodedestimate, ignore_attr = TRUE)

  expect_lt_v(abs(out$est - beta) / out$SE, 3)#assuming normally distributed with SE given by SE above
})

test_that("ppi with minsq weights match estimator1 with fixed beta for ppi_egmodel", {
  set.seed(123)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(model$sample, paramvec = ppi_paramvec(betaL = model$beta[1:2], betap = model$beta[3]),
            method = "closed",
                   bdrythreshold = 1E-10,
            trans = "sqrt", divweight = "minsq", acut = acut)

  hardcodedestimate <- estimator1(model$sample, acut, incb = TRUE, beta = model$beta)

  expect_equal(out$est$paramvec, hardcodedestimate$est$paramvec)
  expect_absdiff_lte_v(out$est$paramvec, model$theta, 2 * out$SE$paramvec)
})

test_that("ppi with prodsq weights match estimator1 with fixed beta for ppi_egmodel", {
  set.seed(123)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(model$sample, ppi_paramvec(betaL = model$beta[1:2], betap = model$beta[3]),
             method = "closed",
                   trans = "sqrt", divweight = "prodsq", acut = acut, 
                   control = list(tol = 1E-12))

  hardcodedestimate <- estimator2(model$sample, acut, incb = TRUE, beta0 = model$beta)

  expect_equal(out$est$paramvec, c(hardcodedestimate$estimator2, model$beta)) 

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
  theta <- ppi_paramvec(AL=ALs, bL=bL, beta=beta)

  set.seed(13418)
  utabl <- rppi(n,beta=beta,AL=ALs,bL=bL,maxden=4)
  u <- utabl[2, ]

  out <- ppi(utabl,
             method = "closed",
                   trans = "sqrt", divweight = "minsq", acut = acut,
                   control = list(tol = 1E-10))

  hardcodedestimate <- estimatorall1(utabl, acut)

  expect_equal(out$est$paramvec, hardcodedestimate$estimator1, ignore_attr = TRUE)

  expect_absdiff_lte_v(out$est$paramvec, theta, 3 * out$SE$paramvec)
})

test_that("ppi with minsq weights match estimatorall1 for ppi_egmodel", {
  set.seed(111)
  model <- ppi_egmodel(100, maxden = 4)
  acut = 0.1
  ppiest <- ppi(Y = model$sample, method = "closed",
                trans = "sqrt", divweight = "minsq", acut = acut)

  hardcodedestimate <- estimatorall1(model$sample, acut)

  expect_equal(ppiest$est$paramvec, hardcodedestimate$estimator1, ignore_attr = TRUE)

  expect_absdiff_lte_v(ppiest$est$paramvec, model$theta, 3 * ppiest$SE$paramvec)
})

test_that("ppi with minsq weights match estimatorall1 for ppi_egmodel, fixed final beta", {
  set.seed(123)
  model <- ppi_egmodel(100, maxden = 4)

  acut = 0.1
  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(betap = tail(model$beta, 1), p = model$p),
             trans = "sqrt", divweight = "minsq", acut = acut)

  hardcodedestimate <- estimatorall1(model$sample, acut, betap = model$beta[model$p])

  expect_equal(out$est$paramvec, t_fu2t(hardcodedestimate$estimator1, ppi_paramvec(betap = tail(model$beta, 1), p = model$p)), ignore_attr = TRUE)

  expect_absdiff_lte_v(out$est$paramvec, model$theta, 3 * out$SE$paramvec)
})

test_that("ppi with minsq weights match estimatorall1 for ppi_egmodel, fixed final beta, large n", {
  skip_on_cran()
  set.seed(123)
  model <- ppi_egmodel(100000, maxden = 4)

  acut = 0.1
  out <- ppi(Y = model$sample, method = "closed",
             paramvec = ppi_paramvec(betap = tail(model$beta, 1), p = model$p),
             trans = "sqrt", divweight = "minsq", acut = acut)

  hardcodedestimate <- estimatorall1(model$sample, acut, betap = model$beta[model$p])

  expect_equal(out$est$paramvec, t_fu2t(hardcodedestimate$estimator1, ppi_paramvec(betap = tail(model$beta, 1), p = model$p)), tolerance = 1E-4)

  expect_absdiff_lte_v(out$est$paramvec, model$theta, 3 * out$SE$paramvec)
})

test_that("ppi with minsq weights performs well on simplex, fixed final beta", {
  set.seed(1234)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(p = model$p, betap = tail(model$beta, 1)),
             trans = "none", divweight = "minsq", acut = acut)

  expect_absdiff_lte_v(out$est$paramvec, model$theta, 3 * out$SE$paramvec)
})

test_that("ppi via cppad matches Score1 for p=5, particularly the order of the off diagonals in ALs", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  expect_warning({prop <- rppi(1000, beta=beta, AL=ALs, bL=bL, maxden=35)}, "maxden")

  acut = 0.1
  est_hardcoded <- estimator1(prop, acut, incb = 0, beta = beta)

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL, beta = beta),
                   method = "closed",
                         trans = "sqrt", acut = acut, divweight = "minsq",
                         control = list(tol = 1E-13))
  expect_equal(est_cppad$est$paramvec, est_hardcoded$est$paramvec,
               ignore_attr = TRUE)

  #also it makes sense that the smo and gradient are v low at the hardcoded estimate
  ppitapes <- buildsmotape("sphere", "ppi",
                           rep(0.1, p), ppi_paramvec(p, bL = bL, beta = beta),
                           weightname = "minsq", acut = acut)
  smvals <- tape_smvalues_wsum(ppitapes$smotape, prop, fromsmatrix(est_hardcoded$est$AL))
  expect_lt(sum(smvals$grad^2), 1E-20)

  # check that rearrangement has large gradient
  smvals <- tape_smvalues_wsum(ppitapes$smotape, prop, fromsmatrix(est_hardcoded$est$AL)[c(1:6, 8, 7, 9, 10)])
  expect_gt(sum(smvals$grad^2), 1E-2)
})
