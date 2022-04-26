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
  p=4
  # smofun <- ptapesmo(c(1,1,1,1,3,3,3,3,3,3,3,3,3), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function
  #sample size
  n=1000

  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.3,p,1)
  beta0[p]=-0.5
  theta <- toPPIparamvec(ALs, bL, beta0)

  set.seed(134)
  utabl <- cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)$samp3
  smofun <- ptapesmo(c(utabl[2, ], 1:length(theta)), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function

  psmo(smofun, utabl[2, ], theta) #very strange that the 4th component of the gradient of the ll on the sphere is zero
  # deriv of log ppi ll on sphere wrt z[4] is:
  # (1 + 2beta[4]) / z[4] + 0 + 0 which is zero because beta0[4] is -0.5
  lltape <- ptapell(p, length(theta), llname = "ppi")
  pJacobian(lltape, utabl[1,], theta)
  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = theta*0,
               fn = function(theta){smobj(smofun, theta, utabl)},
               # gr = function(theta){smobjgrad(smofun, theta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(utabl, acut)
  expect_equal(out$par, directestimate$estimator1, tolerance = 1E-3, ignore_attr = TRUE)
  expect_equal(out$par, theta, tolerance = 1E-3, ignore_attr = TRUE)
})
