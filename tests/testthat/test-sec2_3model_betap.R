####model in Section 2.3 with beta=(-0.8,-0.8,?) here a_c is set to 20%####
# edited from 'modelA.R' in the original code

#### Setup ####
#dimension
p=3




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
# beta0=matrix(-0.8,p,1)
# beta0[p]=0
theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)

#### Tests ####
test_that("Score1ac estimator can estimate beta0[0] large, and others correctly", {
  #sample size
  n=1000
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5

  #simulate sample from PPI model
  set.seed(1)
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)
  samp3=samp1$samp3

  #maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler
  #a few times to check that it is an appropriate upper bound.
  # 4 seems to be pretty good (I've run the above rhybrid many times).
  # I.e. the simulation result doesn't suggest changing maxden=4
  stopifnot(samp1$maxden <= 4)
  maxden <- 4

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model
  estimator=cdabyppi:::estimatorall1(samp3,acut, betap = NULL)
  estimate1all=estimator$estimator1
  expect_true(all(abs(theta - estimate1all[1:5]) <= 20)) #weak because not the target of the test
  #invented bounds for beta0 estimates for now, for the first two components
  expect_true(all(abs(beta0[-p] - estimate1all[6:7]) <= 2*3/sqrt(n)))
  # the final (large component) can't be estimated well, so long as larger then I'll be happy
  expect_gt(estimate1all[8], beta0[p])
})

test_that("Score1ac estimator can estimate beta0[1:(p-1)] for beta0[p] larger than -0.5", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p] = 5

  #simulate sample from PPI model
  set.seed(1)
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)
  samp3=samp1$samp3

  stopifnot(samp1$maxden <= 4)
  maxden <- 4

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model, beta[p] correctly fixed
  estimator=cdabyppi:::estimatorall1(samp3,acut, betap = beta0[p])
  estimate1all=estimator$estimator1
  # SE given as if beta0 fixed
  std1=cdabyppi:::estimator1SE(samp3,acut,estimate1all[1:5, , drop = FALSE],estimator$W_est,1, beta0 = c(estimate1all[6:7], beta0[3]))
  #2*SE bounds
  expect_true(all(abs(theta - estimate1all[1:5]) <= 2*std1))
  #invented bounds for beta0 estimates for now, for the first two components
  expect_true(all(abs(beta0[-p] - estimate1all[6:7]) <= 2*3/sqrt(n)))
})

test_that("Score1ac estimator can estimate beta0[1:(p-1)] for beta0[p] large but misspecified", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p]= 5

  #simulate sample from PPI model
  set.seed(1)
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)
  samp3=samp1$samp3

  stopifnot(samp1$maxden <= 4)
  maxden <- 4

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model, beta[p] correctly fixed
  estimator=cdabyppi:::estimatorall1(samp3,acut, betap = -0.5)
  estimate1all=estimator$estimator1
  # SE given as if beta0 fixed
  std1=cdabyppi:::estimator1SE(samp3,acut,estimate1all[1:5, , drop = FALSE],estimator$W_est,1, beta0 = c(estimate1all[6:7], beta0[3]))
  #2*SE bounds
  expect_true(all(abs(theta - estimate1all[1:5]) <= 2*std1))
  #invented bounds for beta0 estimates for now, for the first two components
  expect_true(all(abs(beta0[-p] - estimate1all[6:7]) <= 2*3/sqrt(n)))
})
