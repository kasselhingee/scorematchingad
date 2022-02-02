# Model 4 from Section A.11 of the JASA paper.
# Simulated, then estimated using 4 estimators Score1, Score2, Score1ac and Score2ac.
# bL should be estimated in all these tests

#dimension
p=10

#sample size
n=1000

#set seed
set.seed(1)

#parameters for the tGaussian model:
muL=matrix(0,p-1,1)
muL[1:sum(p,-1)]=0.04
aa=matrix(1/10000,p-1,1)
D=diag(as.vector(aa))
SigA=D
ALs=-0.5*solve(SigA)
bL=solve(SigA)%*%muL
beta0=matrix(0,p,1)



#simulate sample from the tGaussian model:
samp3=cdabyppi:::rtGaussian(n,p,muL,SigA)


#the response variable prop is the true proportions samp3
prop=samp3

test_that("Score1ac is within 3 SE for 75% of parameters", {
  #a_c for h function:
  acut=0.1

  #calculate scoring estimate:
  estimator=cdabyppi:::estimator1(prop,acut,1, beta0)
  estimate1=estimator$estimator1

  #estimate of W matrix
  W_est=estimator$W_est

  #standard errors for Score1ac
  std1=cdabyppi:::estimator1SE(prop,acut,estimate1,W_est,1, beta0)
  theta = c(diag(ALs), ALs[upper.tri(ALs)], bL)
  expect_gt(mean(abs(theta - estimate1) <= 3*std1), 0.75)
})

test_that("Score2ac is within 3 SE for 75% of parameters", {
  #a_c for h function:
  acut=2e-07

  #calculate scoring estimate:
  estimator=cdabyppi:::estimator2(prop,acut,1, beta0)
  estimate2=estimator$estimator2

  #estimate of W matrix
  W_est=estimator$W_est

  #standard errors for Score2ac
  std2=cdabyppi:::estimator2SE(prop,acut,estimate2,W_est,1, beta0)
  theta = c(diag(ALs), ALs[upper.tri(ALs)], bL)
  expect_gt(mean(abs(theta - estimate2) <= 3*std2), 0.75)
})

test_that("Score2 estimate with 3 Score2ac-SE for 75% of parameters", {
  acut=10

  #calculate scoring estimate:
  estimator=cdabyppi:::estimator2(prop,acut,1, beta0)
  estimate3=estimator$estimator2

  # Get SE with a hack
  std2=cdabyppi:::estimator2SE(prop,acut, estimator$estimator2, estimator$W_est,1, beta0)
  theta = c(diag(ALs), ALs[upper.tri(ALs)], bL)
  expect_gt(mean(abs(theta - estimate3) <= 3*std2), 0.75)
})

test_that("Score1 estimate with 3 Score1ac-SE for 75% of parameters", {
  acut=10

  #calculate scoring estimate:
  estimator=cdabyppi:::estimator1(prop,acut,1, beta0)
  estimate4=estimator$estimator1

  # Get SE with a hack
  std=cdabyppi:::estimator1SE(prop,acut, estimator$estimator1, estimator$W_est,1, beta0)
  theta = c(diag(ALs), ALs[upper.tri(ALs)], bL)
  expect_gt(mean(abs(theta - estimate3) <= 3*std), 0.75)
})
