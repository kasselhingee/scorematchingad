#### Preparing microbiome data #############################################
data("microdata", package = "cdabyppi")
microdata <- microdata[!microdata$IndividualID %in% c(2079, 2280), ] #remove two outlying measurements
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


####Calculate summary statistics####


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


####Reduce dimensions to p=5####


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

#set beta (this is fixed here)
beta0=matrix(0,p,1)
beta0[1]=-0.8
beta0[2]=-0.85
beta0[3]=0
beta0[4]=-0.2
beta0[5]=0

#a_c for h function:
acut=0.01

#### Test including b_L ####
test_that("estimator1 and SE is historically correct with b_L included (article Table 4)", {


  #calculate scoring estimate:
  estimator= cdabyppi:::estimator1(propreal,acut,1, beta0)
  estimate1=estimator$estimator1
  #rearrange to historical ordering
  ordindx <- order(combparam2uppertriorder(length(estimate1) + p)) #the plus p is for the beta that isn't estimated
  ordindx <- ordindx[1:length(estimate1)]
  estimate1 <- estimate1[ordindx]
  dim(estimate1) <- c(length(estimate1), 1)
  #check historically
  expect_snapshot_value(signif(estimate1, 8), style = "json2") #8 is the default number of digits for jsonlite::serializeJSON

  #estimate of W matrix
  W_est=estimator$W_est
  expect_snapshot_value(round(max(W_est), 8), style = "json2") #have to use round here because the json conversion doesn't necessarily show it in scientific notation
  expect_snapshot_value(signif(mean(W_est), 8), style = "json2")
  expect_snapshot_value(signif(which.max(W_est), 8), style = "json")

  #standard errors
  std1=cdabyppi:::estimator1SE(propreal,acut,estimator$estimator1,estimator$W_est,1, beta0)
  #rearrange back to combn ordering
  std1 <- std1[ordindx]
  expect_snapshot_value(signif(std1, 8), style = "json2")

  #estimated parameters
  thetamats <- fromPPIparamvec(c(estimator$estimator1, rep(NA, p)))
  ALs <- thetamats$ALs
  bL <- thetamats$bL
  dim(bL) <- c(length(bL), 1)
  #values in Table 4 in the article:
  expect_snapshot_value(signif(ALs[upper.tri(ALs, diag = TRUE)], 8), style = "json2")
  expect_snapshot_value(signif(bL, 8), style = "json2")
  expect_snapshot_value(signif(estimate1/std1, 8), style = "json2")

})

#### Test omitting b_L ####
test_that("estimator1 and SE is historically correct with b_L ommitted (article table 3)", {

  #calculate scoring estimate:
  estimator=cdabyppi:::estimator1(propreal,acut,0, beta0)
  estimate1=estimator$estimator1
  #rearrange to historical ordering
  ordindx <- order(combparam2uppertriorder_matrix(length(estimate1)))
  estimate1 <- estimate1[ordindx]
  dim(estimate1) <- c(length(estimate1), 1)
  expect_snapshot_value(signif(estimate1, 8), style = "json2")

  #estimate of W matrix
  W_est=estimator$W_est
  expect_snapshot_value(round(max(W_est), 8), style = "json2") #have to use round here because the json conversion doesn't necessarily show it in scientific notation
  expect_snapshot_value(signif(mean(W_est), 8), style = "json2")
  expect_snapshot_value(signif(which.max(W_est), 8), style = "json")

  #standard errors
  std1=estimator1SE(propreal,acut,estimator$estimator1,estimator$W_est,0, beta0)
  std1 <- std1[ordindx] #rearrange to combn ordering for historical comparison

  #estimated parameters
  thetamats <- fromPPIparamvec(c(estimator$estimator1, rep(NA, p-1), rep(NA, p)))
  ALs <- thetamats$ALs

  #values in Table 3 in the article:
  expect_snapshot_value(signif(ALs[upper.tri(ALs, diag = TRUE)], 8), style = "json2")
  expect_snapshot_value(round(estimate1/std1, 8), style = "json2")

})

#### Dirchlet Model ####
test_that("Dirchlet moment fitting is historically correct", {
  alpha_dir=cdabyppi:::dir_moment(propreal)
  expect_snapshot_value(signif(alpha_dir, 8), style = "json2")
})
