# Simulates then estimates for model 7 of section A.11 in the JASA paper


#choose seed
set.seed(1)

#dimension
p=3

#Dirichlet parameters (beta+1)
alpha=matrix(0.50,p,1)
alpha[2]=1.70
alpha[3]=541
alpha-1

#sample size
n=92

#simulate sample from Dirichlet distribution
dirprop=MCMCpack::rdirichlet(n, alpha)

#simulate from Dirichlet-multinomial model
ni=matrix(2000,n,1)
ni=as.vector(ni)
x=matrix(0,n,p)
for (j in 1:n)
{
  x[j,]=rmultinom(1,ni[j],prob=dirprop[j,])
}
prop=x/ni

#use estimated proportions (prop) or true proportions (dirprop) in estimation
#dirfit=prop
dirfit=dirprop

test_that("Dirichlet moment estimator runs", {
  estimate6=dir_moment(dirfit)-1
  expect_gt(mean(abs(estimate6 - alpha) < abs(alpha)) , 0.6)
})

test_that("Score2 runs", {
  #a_c for h function (any large value):
  acut=10


  #calculate scoring estimate:
  estimate3=dir_sqrt_prodh(dirfit,acut)
  expect_gt(mean(abs(estimate3 - alpha) < abs(alpha)) , 0.6)
})


test_that("Score2ac runs", {
  #a_c for h function:
  acut=0.001

  #calculate scoring estimate:
  estimate2=dir_sqrt_prodh(dirfit,acut)
  expect_gt(mean(abs(estimate2 - alpha) < abs(alpha)) , 0.6)
})

test_that("Score1 runs", {
  #a_c for h function (any large value):
  acut=10

  #calculate scoring estimate:
  estimate4=dir_sqrt_minimah(dirfit,acut)
  expect_gt(mean(abs(estimate4 - alpha) < abs(alpha)) , 0.6)
})

test_that("Score1ac runs", {
  #a_c for h function:
  acut=0.01


  #calculate scoring estimate:
  estimate1=dir_sqrt_minimah(dirfit,acut)
  expect_gt(mean(abs(estimate1 - alpha) < abs(alpha)) , 0.6)
})

test_that("Dirichlet score matching estimates are historically correct", {
  acut=0.01

  estimate1=dir_sqrt_minimah(dirfit,acut)
  expect_snapshot_value(signif(estimate1, 8), style = "json2")

  estimate2=dir_sqrt_prodh(dirfit,acut)
  expect_snapshot_value(signif(estimate2, 8), style = "json2")
})



