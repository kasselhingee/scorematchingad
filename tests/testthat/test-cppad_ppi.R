test_that("ppi and dirichlet smo match when AL and bL is zero and p = 3", {
  beta = c(-0.3, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  utabl <- cdabyppi:::rhybrid(10,p,beta,ALs,bL,4)$samp3

  acut = 1
  smoppi <- ptapesmo(c(utabl[2, ], 1:length(theta)), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function
  smodir <- ptapesmo(c(utabl[2, ], 1:length(theta)), p, llname = "dirichlet", manifoldname = "sphere", "minsq", acut = acut)

  ppival <- psmo(smoppi, utabl[2, ], theta)
  dirval <- psmo(smodir, utabl[2, ], beta)
  expect_equal(ppival, dirval)
})

test_that("ppi works when AL and bL is zero and p = 4", {
  beta = c(-0.3, -0.2, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta = toPPIparamvec(ALs, bL, beta)

  utabl <- cdabyppi:::rhybrid(10,p,beta,ALs,bL,4)$samp3

  acut = 1
  smoppi <- ptapesmo(c(utabl[2, ], 1:length(theta)), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function
  smodir <- ptapesmo(c(utabl[2, ], 1:length(theta)), p, llname = "dirichlet", manifoldname = "sphere", "minsq", acut = acut)

  ppival <- psmo(smoppi, utabl[2, ], theta)
  dirval <- psmo(smodir, utabl[2, ], beta)
  expect_equal(ppival, dirval)


  out <- optim(par = theta*0,
               fn = function(theta){smobj(smoppi, theta, utabl)},
               gr = function(theta){smobjgrad(smoppi, theta, utabl)},
               method = "BFGS")
  cppadest <- fromPPIparamvec(out$par, p)

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator1_dir(utabl, acut)
  expect_equal(cppadest$beta0, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})

test_that("ppi with minsq weights match estimator1 for p = 4", {
  acut = 1
  # smofun <- ptapesmo(c(1,1,1,1,3,3,3,3,3,3,3,3,3), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function
  #sample size
  n=1000

  #parameters for the PPI model
  beta = c(-0.3, -0.2, -0.1, 3)
  p = length(beta)
  ALs = matrix(0, nrow = p-1, ncol = p-1)
  ALs[1, 2] <- 1
  ALs[2, 1] <- 1
  bL = matrix(0, nrow = p-1, ncol = 1)
  theta <- toPPIparamvec(ALs, bL, beta)

  set.seed(134)
  utabl <- cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)$samp3
  u <- utabl[2, ]
  smofun <- ptapesmo(c(u, 1:length(theta)), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function

  psmo(smofun, u, theta) #very strange that the 4th component of the gradient of the ll on the sphere is zero
  # deriv of log ppi ll on sphere wrt z[1] is:   deriv(z2 * (z2[2], z2[1], 0, 0)) = deriv(z[1]^2 * z[2]^2 + z[2]^2 + z[1]^2)
  # (1 + 2beta[1]) / z[1] + 2z[1]z[2]^2 + 0
  (1 + 2 * beta[1]) / sqrt(u[1]) + 2 * sqrt(u[1]) * u[2]

  ppill_S_r <- function(z, beta0, ALs, bL){
    p=length(u)
    A = matrix(0, nrow = p, ncol = p)
    A[1:p-1, 1:p-1] = ALs
    b = matrix(c(bL, 0), nrow = p, ncol =1)
    out = 2 + sum((1 + 2 * beta0) * log(z)) + t(z^2) %*% A %*% (z^2) + t(b) %*% (z^2)
    return(out)
  }
  z = sqrt(u)
  numericDeriv(quote(ppill_S_r(z, beta, ALs, bL)), c("z"))

  lltape <- ptapell(p, length(theta), llname = "ppi")
  pJacobian(lltape, utabl[1,], theta)
  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = theta,
               fn = function(theta){smobj(smofun, theta, utabl)},
               gr = function(theta){smobjgrad(smofun, theta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(utabl, acut)
  expect_equal(out$par, directestimate$estimator1, tolerance = 1E-3, ignore_attr = TRUE)
  expect_equal(out$par, theta, tolerance = 1E-3, ignore_attr = TRUE)
})
