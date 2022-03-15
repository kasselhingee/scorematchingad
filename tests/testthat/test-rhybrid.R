# test that rhybrid old matches old rhybrid new

#dimension
p=3

#sample size
n=5000

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

test_that("current PPI simulation method gives samples with similar empirical density estimates as the original simulation method", {
  # simulate using old method
  samp2 <- cdabyppi:::rhybrid_singly(n,p,beta0,ALs,bL,4)
  H <- ks::Hpi(samp2$samp3[, -p])
  kde_historic <- ks::kde(samp2$samp3[, -p], H)
  #simulate sample from PPI model
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL,4)
  H <- ks::Hpi(samp1$samp3[, -p])
  kde_current <- ks::kde(samp1$samp3[, -p], H)

  testresult <- ks::kde.test(samp2$samp3[, -p], samp1$samp3[, -p])
  expect_gte(testresult$pvalue, 0.01)
})

