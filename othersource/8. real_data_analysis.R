##########################################################
##Input data
##########################################################

library(readxl)
microdata <- read_excel("microbiome_data.xlsx")
View(microdata)
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

##########################################################
##Histograms: Figure 3 in the article
#########################################################

par(mfrow = c(2,3),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.0)
hist(propreal[,1],nclass=100,main="",xlab="(a)")
hist(propreal[,2],nclass=100,main="",xlab="(b)")
hist(propreal[,3],nclass=100,main="",xlab="(c)")
hist(propreal[,4],nclass=100,main="",xlab="(d)")
hist(propreal[,5],nclass=100,main="",xlab="(e)")
#dev.off()



####################################################################
##Score matching estimator plus standard errors: Table 4
####################################################################

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
estimator=estimator1(propreal,acut,1)
estimate1=estimator$estimator1
print(estimate1)

#estimate of W matrix
W_est=estimator$W_est

#standard errors
std1=estimator1SE(propreal,acut,estimate1,W_est,1)
print(std1)

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

#values in Table 4 in the article:
signif(ALs,6)
signif(bL,6)
signif(estimate1/std1,3)


####################################################################
##Score matching estimator plus standard errors: Table 3
####################################################################


#calculate scoring estimate:
estimator=estimator1(propreal,acut,0)
estimate1=estimator$estimator1
print(estimate1)

#estimate of W matrix
W_est=estimator$W_est

#standard errors
std1=estimator1SE(propreal,acut,estimate1,W_est,0)
print(std1)

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


#values in Table 3 in the article:
signif(ALs,6)
signif(estimate1/std1,3)



####################################################################
##Model diagnostics
####################################################################

#sample size
n=100000

#set seed
set.seed(1)


#simulate sample from the fitted PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,0.2)

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)

#simulated sample:
samp3=samp1$samp3


#rounding simulated data
samp4=round(2000*samp3,0)/2000
samp3=samp4


#Figure 4 in article
par(mfrow = c(2,3),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.5)
qqplot(propreal[,1],samp3[,1],xlab="(a)",ylab="")
qqplot(propreal[,2],samp3[,2],xlab="(b)",ylab="")
qqplot(propreal[,3],samp3[,3],xlab="(c)",ylab="")
qqplot(propreal[,4],samp3[,4],xlab="(d)",ylab="")
plot(propreal[,1],propreal[,3],xlab="(e)",ylab="")
plot(samp3[1:1000,1],samp3[1:1000,3],cex=1,xlab="(f)",ylab="")
#dev.off()


#fit a 5D Dirichlet, simulate and then round
alpha_dir=dirichmom(propreal)
dirfit=rdirichlet(n, alpha_dir)
dirfit2=round(2000*dirfit,0)/2000
dirfit=dirfit2


#Table 5 in the article
ks.test(samp3[,1],propreal[,1])
ks.test(dirfit[,1],propreal[,1])
ks.test(samp3[,2],propreal[,2])
ks.test(dirfit[,2],propreal[,2])
ks.test(samp3[,3],propreal[,3])
ks.test(dirfit[,3],propreal[,3])
ks.test(samp3[,4],propreal[,4])
ks.test(dirfit[,4],propreal[,4])
ks.test(samp3[,5],propreal[,5])
ks.test(dirfit[,5],propreal[,5])


#Table 6 in the article
mean(propreal[,1])
mean(samp3[,1])
mean(dirfit[,1])
mean(propreal[,2])
mean(samp3[,2])
mean(dirfit[,2])
mean(propreal[,3])
mean(samp3[,3])
mean(dirfit[,3])
mean(propreal[,4])
mean(samp3[,4])
mean(dirfit[,4])
mean(propreal[,5])
mean(samp3[,5])
mean(dirfit[,5])
sd(propreal[,1])
sd(samp3[,1])
sd(dirfit[,1])
sd(propreal[,2])
sd(samp3[,2])
sd(dirfit[,2])
sd(propreal[,3])
sd(samp3[,3])
sd(dirfit[,3])
sd(propreal[,4])
sd(samp3[,4])
sd(dirfit[,4])
sd(propreal[,5])
sd(samp3[,5])
sd(dirfit[,5])


