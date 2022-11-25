# genscore has a function gen() for simulating data

test_that("toAstar() works on an example matrix", {
  #dimension
  p=3

  #sample size
  n=100

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
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5

  expect_error(Astar <- toAstar(ALs, bL))

  # check conversion
  samp <- c(0.007, 0.05, 1 - 0.007 - 0.05)
  hardcoded <- t(samp[-p]) %*% ALs %*% samp[-p] + t(bL) %*% samp[-p]
  # expect_equal(hardcoded, t(samp) %*% Astar %*% samp)
})

