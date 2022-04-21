test_that("ppi with minsq weights match estimator1 for p = 3", {
  acut = 0.1
  p=3
  smofun <- ptapesmo(c(1,1,1,3,3,3,3,3,3), p, llname = "ppi", manifoldname = "sphere", "minsq", acut = acut) #tape of the score function
  #sample size
  n=10

  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  # muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], beta0)

  set.seed(134)
  utabl <- cdabyppi:::rhybrid_singly(n,p,beta0,ALs,bL,4)$samp3

  smobj(smofun, theta, utabl)
  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = theta*0,
               fn = function(theta){smobj(smofun, theta, utabl)},
               gr = function(theta){smobjgrad(smofun, theta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimatorall1(utabl, acut)
  expect_equal(out$par, directestimate$estimator1, tolerance = 1E-3, ignore_attr = TRUE)
})
