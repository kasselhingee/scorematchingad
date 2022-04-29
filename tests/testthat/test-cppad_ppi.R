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
}) #passing

test_that("cppad ppi estimate works when AL and bL is zero and p = 4", {
  beta = c(-0.3, -0.2, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  set.seed(123)
  utabl <- cdabyppi:::rhybrid(1000,p,beta,ALs,bL,4)$samp3

  acut = 0.1
  psphere = pmanifold("sphere")
  pdir <- ptapell(rep(0.1, p), beta, llname = "dirichlet", psphere, fixedtheta = rep(FALSE, p), verbose = FALSE)
  pppi <- ptapell(rep(0.1, p), theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, p), theta, pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function
  smodir <- ptapesmo(rep(0.1, p), beta, pll = pdir, pman = psphere, "minsq", acut = acut, verbose = FALSE)

  # it looks like the taped function above is not altering bL or beta
  # potentially the ordering of the theta values is wrong??
  out <- optim(par = theta * 0 + 1,
               fn = function(theta){smobj(smoppi, theta, utabl)},
               gr = function(theta){smobjgrad(smoppi, theta, utabl)},
               method = "BFGS", control =list(maxit = 1000))

  # smoperpt <- apply(utabl, MARGIN = 1, FUN = function(u){pForward0(smoppi, out$par, u)})


  cppadest <- fromPPIparamvec(out$par, p)

  expect_equal(pForward0(pdir, utabl[2, ], beta), pForward0(pppi, utabl[2, ], theta))
  expect_equal(pJacobian(pdir, utabl[2, ], beta), pJacobian(pppi, utabl[2, ], theta))
  expect_equal(pForward0(smoppi, theta, utabl[2, ]), pForward0(smodir, beta, utabl[2, ]))

  directestimate <- estimator1_dir(utabl, acut)

  # there is a difference in direct estimates because the direct estimate smobj value is poorer:
  expect_lte(smobj(smoppi, out$par, utabl),
             smobj(smoppi, c(rep(0, length(theta) - p), directestimate), utabl))

  SE <- smestSE(smoppi, out$par, utabl) #SE is error to true parameters, here using it also as a proxy to optimisation accuracy
  cppadestSE <- fromPPIparamvec(diag(SE), p)
  expect_true(all(abs(cppadest$beta0 - directestimate) / cppadestSE$beta0 < 0.1)) #proxy for optimisation flatness

  expect_true(all(abs(out$par - theta) / diag(SE) < 2)) #assuming normally distributed with SE given by SE above
}) #failing

test_that("ppi with minsq weights match estimator1 with fixed beta for more complex model", {
  set.seed(123)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere,
                  fixedtheta = c(rep(FALSE, length(model$theta) - model$p), rep(TRUE, model$p)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - model$p), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = model$theta[1:(length(model$theta) - model$p)] * 0,
               fn = function(theta){smobj(smoppi, theta, model$sample)},
               gr = function(theta){smobjgrad(smoppi, theta, model$sample)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator1(model$sample, acut, incb = TRUE, beta0 = model$beta0)

  SE <- smestSE(smoppi, out$par, model$sample) #SE is error to true parameters, here using it also as a proxy to optimisation accuracy
  expect_true(all(abs(out$par - directestimate$estimator1) / diag(SE) < 0.5)) #proxy for optimisation flatness

  expect_true(all(abs(out$par - model$theta[1:(length(model$theta) - model$p)]) / diag(SE) < 2)) #assuming normally distributed with SE given by SE above
}) #passing

test_that("ppi with minsq weights match estimatorall1 for p = 4 simple", {
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

  set.seed(134)
  utabl <- cdabyppi:::rhybrid(n,p,beta,ALs,bL,4)$samp3
  u <- utabl[2, ]
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, p), theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, p), theta, pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = theta,
               fn = function(theta){smobj(smoppi, theta, utabl)},
               gr = function(theta){smobjgrad(smoppi, theta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(utabl, acut)

  SE <- smestSE(smoppi, out$par, utabl) #SE is error to true parameters, here using it also as a proxy to optimisation accuracy
  expect_true(all(abs(out$par - directestimate$estimator1) / diag(SE) < 0.5)) #proxy for optimisation flatness

  expect_true(all(abs(out$par - theta) / diag(SE) < 2)) #assuming normally distributed with SE given by SE above
}) #passing

test_that("ppi with minsq weights match estimatorall1 for p = 3, more complex model", {
  set.seed(123)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(model$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:length(model$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = model$theta * 0,
               fn = function(theta){smobj(smoppi, theta, model$sample)},
               gr = function(theta){smobjgrad(smoppi, theta, model$sample)},
               method = "BFGS",
               control = list(maxit = 10000))
  stopifnot(out$convergence == 0)

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(model$sample, acut)

  expect_lte(smobj(smoppi, out$par, model$sample),
             smobj(smoppi, directestimate$estimator1, model$sample))

  SE <- smestSE(smoppi, out$par, model$sample) #SE is error to true parameters, here using it also as a proxy to optimisation accuracy
  expect_true(all(abs(out$par - directestimate$estimator1) / diag(SE) < 0.5)) #proxy for optimisation flatness

  expect_true(all(abs(out$par - model$theta) / diag(SE) < 2)) #assuming normally distributed with SE given by SE above
}) #failing

test_that("ppi with minsq weights match estimatorall1 for p = 3, more complex model, fixed final beta", {
  set.seed(123)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = model$theta[-length(model$theta)] * 0,
               fn = function(theta){smobj(smoppi, theta, model$sample)},
               gr = function(theta){smobjgrad(smoppi, theta, model$sample)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(model$sample, acut, betap = model$beta0[model$p])

  SE <- smestSE(smoppi, out$par, model$sample) #SE is error to true parameters, here using it also as a proxy to optimisation accuracy
  expect_true(all(abs(out$par - directestimate$estimator1) / diag(SE) < 0.5)) #proxy for optimisation flatness

  expect_true(all(abs(out$par - model$theta[-length(model$theta)]) / diag(SE) < 2)) #assuming normally distributed with SE given by SE above
}) #failing

test_that("ppi with minsq weights match estimatorall1 for p = 4, more complex model", {
  set.seed(123)
  model <- sec2_3model_p4(1000, maxden = 8)

  acut = 0.1
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, model$p), model$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(model$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:length(model$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = model$theta * 0,
               fn = function(theta){smobj(smoppi, theta, model$sample)},
               gr = function(theta){smobjgrad(smoppi, theta, model$sample)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(model$sample, acut)

  SE <- smestSE(smoppi, out$par, model$sample) #SE is error to true parameters, here using it also as a proxy to optimisation accuracy
  expect_true(all(abs(out$par - directestimate$estimator1) / diag(SE) < 0.5)) #proxy for optimisation flatness

  expect_true(all(abs(out$par - model$theta) / diag(SE) < 2)) #assuming normally distributed with SE given by SE above
}) #failing
