#############################################
#preparing microbiome data
#############################################
data("microdata", package = "cdabyppi")
countdata=as.matrix(microdata[,12:31])

#sample size
n=92

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
##Calculate summary statistics
##########################################################


#means (relative abundances):
mprop=matrix(0,1,p)
for (j in 1:p)
{
	mprop[1,j]=mean(prop[,j])
}
#Actinobacteria:
mprop[2]
#Proteobacteria:
mprop[13]
#TM7
mprop[16]
#Cyanobacteria/Chloroplast
mprop[5]

#proportion of zeros in each category
pzero=matrix(0,1,p)
for (j in 1:p)
{
	for (k in 1:n)
	{
		if (countdata[k,j]==0){pzero[1,j]=pzero[1,j]+1}
	}
}
pzero=pzero/n
#Actinobacteria:
pzero[2]
#Proteobacteria:
pzero[13]
#TM7:
pzero[16]
#Cyanobacteria/Chloroplast:
pzero[5]


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



test_that("estimator1 and SE is historically correct with b_L included (article Table 4)", {
  #set beta (this is fixed here)
  beta0=matrix(0,p,1)
  beta0[1]=-0.8
  beta0[2]=-0.85
  beta0[3]=0
  beta0[4]=-0.2
  beta0[5]=0

  #a_c for h function:
  acut=0.01

  #calculate scoring estimate:
  estimator= cdabyppi:::estimator1(propreal,acut,1, beta0)
  expect_snapshot_value(estimator, style = "serialize")
  estimate1=estimator$estimator1

  #estimate of W matrix
  W_est=estimator$W_est

  #standard errors
  std1=estimator1SE(propreal,acut,estimate1,W_est,1, beta0)
  expect_snapshot_value(std1, style = "serialize")

  #estimated parameters
  ALs=matrix(0,p-1,p-1)
  bL=matrix(0,p-1,1)
  sp=p-1
  xx=c(1:sp)
  ind=combn(xx, 2, FUN = NULL, simplify = TRUE)
  qind=length(ind[1,])
  for (j in 1:sp)
  {
  	ALs[j,j]=estimate1[j]
  }
  true2=estimate1[p:sum(qind,sp)]
  for (j in 1:qind)
  {
  	ALs[ind[1,j],ind[2,j]]=true2[j]
  	ALs[ind[2,j],ind[1,j]]=true2[j]
  }
  bnum=sum(qind,sp,1)
  btot=sum(qind,sp,sp)
  k=1
  for (j in bnum:btot)
  {
  	bL[k]=estimate1[j]
  	k=k+1
  }

})

# test_that("estimator1 and SE is historically correct with b_L ommitted (article table 3)", {

