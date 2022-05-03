test_that("ppi and dirichlet smo value match when AL and bL is zero and p = 3", {
  beta = c(-0.3, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  utabl <- cdabyppi:::rhybrid(10,p,beta,ALs,bL,4)$samp3

  acut = 0.1
  psphere = pmanifold("sphere")
  pdir <- ptapell(c(0.1,0.1,0.1), c(1,2,3), llname = "dirichlet", psphere, fixedtheta = c(FALSE, FALSE, FALSE), verbose = FALSE)
  pppi <- ptapell(c(0.1,0.1,0.1), theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)
  smoppi <- ptapesmo(c(0.1,0.1,0.1), theta, pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  smodir <- ptapesmo(c(0.1,0.1,0.1), c(1,2,3), pll = pdir, pman = psphere, "minsq", acut = acut, verbose = FALSE)

  ppival <- pForward0(smoppi, theta, utabl[2, ])
  dirval <- pForward0(smodir, beta, utabl[2, ])
  expect_equal(ppival, dirval)
})


test_that("cppad ppi estimate works when AL and bL is zero and p = 4", {
  beta = c(-0.3, -0.2, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  set.seed(1234)
  utabl <- cdabyppi:::rhybrid(100,p,beta,ALs,bL,4)$samp3

  acut = 0.1
  psphere = pmanifold("sphere")
  pdir <- ptapell(rep(0.1, p), beta, llname = "dirichlet", psphere, fixedtheta = rep(FALSE, p), verbose = FALSE)
  pppi <- ptapell(rep(0.1, p), theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, p), theta, pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  smodir <- ptapesmo(rep(0.1, p), beta, pll = pdir, pman = psphere, "minsq", acut = acut, verbose = FALSE)

  # it looks like the taped function above is not altering bL or beta
  # potentially the ordering of the theta values is wrong??
  out <- smest(smoppi, theta * 0 + 1, utabl, list(tol = 1E-10))
  stopifnot(out$convergence == 0)

  cppadest <- fromPPIparamvec(out$par, p)

  expect_equal(pForward0(pdir, utabl[2, ], beta), pForward0(pppi, utabl[2, ], theta))
  expect_equal(pJacobian(pdir, utabl[2, ], beta), pJacobian(pppi, utabl[2, ], theta))
  expect_equal(pForward0(smoppi, theta, utabl[2, ]), pForward0(smodir, beta, utabl[2, ]))

  directestimate <- estimator1_dir(utabl, acut)

  # there is a difference in direct estimates because the direct estimate smobj value is poorer:
  expect_lt(out$value,
             smobj(smoppi, c(rep(0, length(theta) - p), directestimate), utabl) + 1E-3 * abs(out$value)) #larger tolerance some apriori information used
  expect_lt(out$gradsize,
               sum(smobjgrad(smoppi, c(rep(0, length(theta) - p), directestimate), utabl)^2))

  cppadestSE <- fromPPIparamvec(out$SE, p)
  cdabyppi:::expect_lt_v(abs(cppadest$beta0 - directestimate)[1:3] / cppadestSE$beta0[1:3], 1)
  expect_lt(abs(cppadest$beta0 - directestimate)[4] / cppadestSE$beta0[4], 3) #largest beta is hard to estimate well

  cdabyppi:::expect_lt_v(abs(out$par - theta) / out$SE, 3)#assuming normally distributed with SE given by SE above
})

test_that("ppi with minsq weights match estimator1 with fixed beta for sec2_3model", {
  set.seed(123)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere,
                  fixedtheta = c(rep(FALSE, length(model$theta) - model$p), rep(TRUE, model$p)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - model$p), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, model$theta[1:(length(model$theta) - model$p)] * 0, model$sample, list(tol = 1E-20))

  directestimate <- estimator1(model$sample, acut, incb = TRUE, beta0 = model$beta0)

  expect_equal(out$value,
               smobj(smoppi, directestimate$estimator1, model$sample))
  expect_equal(out$gradsize,
               sum(smobjgrad(smoppi, directestimate$estimator1, model$sample)^2))

  cdabyppi:::expect_lt_v(abs(out$par - directestimate$estimator1) / out$SE,  0.01) #proxy for optimisation flatness

  cdabyppi:::expect_lt_v(abs(out$par - model$theta[1:(length(model$theta) - model$p)]) / out$SE, 2)
})

test_that("ppi with prodsq weights match estimator1 with fixed beta for sec2_3model", {
  set.seed(123)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere,
                  fixedtheta = c(rep(FALSE, length(model$theta) - model$p), rep(TRUE, model$p)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - model$p), pll = pppi, pman = psphere, "prodsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, model$theta[1:(length(model$theta) - model$p)] * 0, model$sample, list(tol = 1E-20))

  directestimate <- estimator2(model$sample, acut, incb = TRUE, beta0 = model$beta0)

  expect_equal(out$value,
               smobj(smoppi, directestimate$estimator2, model$sample))
  expect_equal(out$gradsize,
               sum(smobjgrad(smoppi, directestimate$estimator2, model$sample)^2))

  cdabyppi:::expect_lt_v(abs(out$par - directestimate$estimator2) / out$SE,  0.01) #proxy for optimisation flatness

  cdabyppi:::expect_lt_v(abs(out$par - model$theta[1:(length(model$theta) - model$p)]) / out$SE, 3)
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
  utabl <- cdabyppi:::rhybrid(n,p,beta,ALs,bL,4)$samp3
  u <- utabl[2, ]
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, p), theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, p), theta, pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, theta * 0, utabl, control = list(tol = 1E-10))

  directestimate <- estimatorall1(utabl, acut)

  expect_equal(out$value,
               smobj(smoppi, directestimate$estimator1, utabl),
               tolerance = 1E-2)
  expect_equal(out$gradsize,
               sum(smobjgrad(smoppi, directestimate$estimator1, utabl)^2))

  cdabyppi:::expect_lt_v(abs(out$par - directestimate$estimator1) / out$SE, 1) #proxy for optimisation flatness

  cdabyppi:::expect_lt_v(abs(out$par - theta) / out$SE, 3)
})

test_that("ppi with minsq weights match estimatorall1 for sec2_3model", {
  set.seed(111)
  model <- sec2_3model(100, maxden = 4)

  acut = 0.1

  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(model$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:length(model$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, model$theta * 0 + 1, model$sample, control = list(tol = 1E-15))

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(model$sample, acut)

  expect_lt(out$value,
            smobj(smoppi, directestimate$estimator1, model$sample) + 1E-5 * abs(out$value))

  cdabyppi:::expect_lt_v(abs(out$par - directestimate$estimator1) / out$SE, 1) #proxy for optimisation flatness

  cdabyppi:::expect_lt_v(abs(out$par - model$theta) / out$SE, 3)
})

test_that("ppi with minsq weights match estimatorall1 for sec2_3model, fixed final beta", {
  set.seed(123)
  model <- sec2_3model(100, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, model$theta[-length(model$theta)] * 0, model$sample,
               control = list(tol = 1E-10))

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(model$sample, acut, betap = model$beta0[model$p])

  expect_lte(out$value,
             smobj(smoppi, directestimate$estimator1, model$sample) + 1E-3 * abs(out$value))

  cdabyppi:::expect_lt_v(abs(out$par - directestimate$estimator1) / out$SE, 0.01) #proxy for optimisation flatness

  cdabyppi:::expect_lt_v(abs(out$par - model$theta[-length(model$theta)]) / out$SE, 3)
})

test_that("ppi with minsq weights performs well on Rpos, fixed final beta", {
  set.seed(1234)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1
  pRpos <- pmanifold("Rpos")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", pRpos, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = pRpos, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, model$theta[-length(model$theta)] * 0, model$sample,
               control = list(tol = 1E-15))

  cdabyppi:::expect_lt_v(abs(out$par - model$theta[-length(model$theta)]) / out$SE, 3)
})

