####model in Section 2.3 with beta=(-0.8,-0.8,-0.5) here a_c is set to 20%####
# from 'modelA.R' in the original code

test_that("Score1ac estimator works on highly concentrated data, with some components close to the boundary", {
  #dimension
  p=4

  #sample size
  n=200

  #set seed
  set.seed(1)

  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*3
  SigA[2,2]=SigA[2,2]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  SigA[3,1]=SigA[1,3]=cor*sqrt(SigA[3,3]*SigA[1,1])
  SigA[3,2]=SigA[2,3]=cor*sqrt(SigA[3,3]*SigA[2,2])
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5

  #simulate sample from PPI model
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)

  #maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler
  #a few times to check that it is an appropriate upper bound.
  # 4 seems to be pretty good (I've run the above rhybrid many times).
  # I.e. the simulation result doesn't suggest changing maxden=4
  stopifnot(samp1$maxden <= 4)
  maxden <- 4

  #simulated sample:
  samp3=samp1$samp3

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model (only beta[p] fixed at -0.5):
  estimator=cdabyppi:::estimatorall1(samp3,acut,betap = -0.5)
  estimate1all=estimator$estimator1
  # use SE estimates as if beta0 was fixed at the estimate (not estimated)
  std1=cdabyppi:::estimator1SE(samp3,acut,estimate1all[1:9, , drop = FALSE],estimator$W_est,1, beta0 = c(estimate1all[10:12], beta0[p]))
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)
  #2*SE bounds
  expect_true(all(abs(theta - estimate1all[1:9]) <= 2*std1))
  message("Misuse of estimator1SE for estimatorall1 results in this test")
  #invented bounds for beta0 estimates for now
  expect_true(all(abs(beta0[-p] - estimate1all[10:12]) <= 2*3/sqrt(n)))

  #calculate scoring estimate with beta fixed at beta0:
  estimator=cdabyppi:::estimator1(samp3,acut,1, beta0)
  estimate1=estimator$estimator1
  std1=cdabyppi:::estimator1SE(samp3,acut,estimate1,estimator$W_est,1, beta0)
  # check
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)
  #2*SE bounds
  expect_true(all(abs(theta - estimate1) < 2*std1))
})

