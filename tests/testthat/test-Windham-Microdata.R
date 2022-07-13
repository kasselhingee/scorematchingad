

test_that("windham_diff estimator matches historical results on dataset with Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria and pooled", {

data("microdata", package = "cdabyppi")
countdata=as.matrix(microdata[,12:31])

#sample size
n=94

#dimension
p=20

#calculate totals
tot=matrix(0,n,1)
for (j in 1:p)
{
 tot=tot+countdata[,j]
}
tot=as.vector(tot)

#proportion data
prop=countdata
for (j in 1:n)
{
	prop[j,]=countdata[j,]/tot[j]
}



##########################################################
##Reduce dimensions to p=5
##########################################################


#calculate 5D dataset
comb=matrix(0,n,5)
comb[,1]=prop[,"TM7"]
comb[,2]=prop[,"Cyanobacteria/Chloroplast"]
comb[,3]=prop[,"Actinobacteria"]
comb[,4]=prop[,"Proteobacteria"]
comb[,5]=abs(1-comb[,1]-comb[,2]-comb[,4]-comb[,3])
propreal=comb


#dimension
p=5

#save data
proprealA=propreal



#proportion of zeros in each category
pzero=matrix(0,1,p)
for (j in 1:p)
{
	for (k in 1:n)
	{
		if (propreal[k,j]==0){pzero[1,j]=pzero[1,j]+1}
	}
}
pzero=pzero/n
pzero*100


##############################
##Estimation
####################################


#initial values for robust estimators
ALs=matrix(0,p-1,p-1)
bL=matrix(0,p-1,1)
ALs[1,1]= -127480.0929
ALs[1,2]= 14068.39057
ALs[1,3]= 1782.261826
ALs[1,4]=  -240.076568
ALs[2,2]= -8191.17253
ALs[2,3]=  -8.002680
ALs[2,4]= 374.693979
ALs[3,3]= -46.638659
ALs[3,4]= 9.027633
ALs[4,4]= -39.208915
ALs[lower.tri(ALs)] <- t(ALs)[lower.tri(ALs)]
beta0=matrix(0,p,1)
beta0[1]=-0.80
beta0[2]=-0.85
beta0[3]=0
beta0[4]=-0.2
beta0[5]=0
bL_est=bL
ALs_est=ALs
beta0_est=beta0
sp=p-1

#try simulating to see required maxden
expect_silent(sim <- rppi(1000, p, beta0, ALs, bL, maxden = 0))

#non-robust components (k^*=2 here)
ind_weightA=matrix(0,sp,1)
ind_weightA[3]=1
ind_weightA[4]=1

#calculate robust estimates
cW=0.7
est1=windham_diff(propreal,cW,ALs_est,bL_est,beta0_est, ind_weightA, originalcorrectionmethod = TRUE)
#estimate of A_L:
expect_snapshot_value(signif(est1$est$ALs,6), style = "json2")
#estimate of beta:
expect_snapshot_value(signif(est1$est$beta,6), style = "json2")

#Use Kassel's correction method
est2=windham_diff(propreal,cW,ALs_est,bL_est,beta0_est, ind_weightA, originalcorrectionmethod = FALSE)
expect_equal(est1$est$ALs, est2$est$ALs, tolerance = 1E-5)
expect_equal(est1$est$beta, est2$est$beta, tolerance = 1E-5)
})



test_that("windham_diff estimator matches historical results on dataset with Spirochates, Verrucomicrobia, Cyanobacteria/Chloroplast, TM7 and pooled", {
  data("microdata", package = "cdabyppi")
  countdata=as.matrix(microdata[,12:31])

  #sample size
  n=94

  #dimension
  p=20

  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
    tot=tot+countdata[,j]
  }
  tot=as.vector(tot)

  #proportion data
  prop=countdata
  for (j in 1:n)
  {
    prop[j,]=countdata[j,]/tot[j]
  }

  ##########################################################
  ##Reduce dimensions to p=5
  ##########################################################


  #dimension
  p=5

  #calculate 5D dataset
  comb=matrix(0,n,p)
  comb[,1]=prop[,14]
  comb[,2]=prop[,18]
  comb[,3]=prop[,5]
  comb[,4]=prop[,16]
  for (j in 1:sum(p,-1))
  {
    comb[,p]=comb[,p]+comb[,j]
  }
  comb[,p]=1-comb[,p]

  #save data
  propreal=comb
  proprealA=propreal

  ##############################
  ##Estimation
  ####################################

  #initial values for robust estimators
  ALs=matrix(-10000,p-1,p-1)
  bL=matrix(0,p-1,1)
  beta0=matrix(0,p,1)
  beta0[1]=-0.80
  beta0[2]=-0.80
  beta0[3]=-0.80
  beta0[4]=-0.80
  bL_est=bL
  ALs_est=ALs
  beta0_est=beta0

  #non-robust components (k^*=4 here)
  ind_weightA=matrix(0,p-1,1)

  #calculate robust estimates
  cW=1.25
  est1=windham_diff(propreal,cW,ALs_est,bL_est,beta0_est, ind_weightA = ind_weightA, originalcorrectionmethod = TRUE)
  #estimate of A_L:
  expect_snapshot_value(signif(est1$est$ALs,6), style = "json2")
  #estimate of beta:
  expect_snapshot_value(signif(est1$est$beta,6), style = "json2")

  #Use Kassel's correction method
  est2=windham_diff(propreal,cW,ALs_est,bL_est,beta0_est, ind_weightA, originalcorrectionmethod = FALSE)
  expect_equal(est1$est$ALs, est2$est$ALs, tolerance = 1E-4)
  expect_equal(est1$est$beta, est2$est$beta, tolerance = 1E-5)
})

