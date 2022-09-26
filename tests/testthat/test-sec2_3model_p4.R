####model in Section 2.3 with beta=(-0.8,-0.8,-0.8,-0.5) here a_c is set to 20%####
# from 'modelA.R' in the original code

#dimension
p=4

#sample size
n=200

#set seed
# set.seed(1)

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

test_that("Score1ac estimator of A, b and beta works on highly concentrated data, with some components close to the boundary", {
  #simulate sample from PPI model
  set.seed(1210)
  samp3=rppi(n,beta0,ALs,bL,maxden=4)

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model (only beta[p] fixed at -0.5):
  estimator=cdabyppi:::estimatorall1(samp3,acut,betap = -0.5)
  estimate1all=estimator$estimator1
  # use SE estimates as if beta0 was fixed at the estimate (not estimated)
  std1=cdabyppi:::estimator1SE(samp3,acut,estimate1all[1:9, , drop = FALSE],estimator$W_est,1, beta0 = c(estimate1all[10:12], beta0[p]))
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)
  #2*SE bounds for 75% of parameters
  expect_gte(mean(abs(theta - estimate1all[1:9]) <= 2*std1), 0.75)
  message("Misuse of estimator1SE for estimatorall1 results in this test")
  #invented bounds for beta0 estimates for now
  expect_true(all(abs(beta0[-p] - estimate1all[10:12]) <= 2*3/sqrt(n)))

  #calculate scoring estimate with beta fixed at beta0:
  estimator=cdabyppi:::estimator1(samp3,acut,1, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  std1=estimator$SE$paramvec
  # check
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)
  #2*SE bounds for 75% of parameters
  expect_gte(mean(abs(theta - estimate1[1:length(theta)]) < 2*std1[1:length(theta)]), 0.75)
})

test_that("Score1ac estimator of A and b only (beta fixed) works on highly concentrated data, with some components close to the boundary", {
  #simulate sample from PPI model
  set.seed(124)
  samp3=rppi(n,beta0,ALs,bL,maxden=4)

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate with beta fixed at beta0:
  estimator=cdabyppi:::estimator1(samp3,acut,1, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  std1=estimator$SE$paramvec
  # check
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)
  #2*SE bounds
  expect_true(all(abs(theta - estimate1[1:length(theta)]) < 2*std1[1:length(theta)]))
})
