####################################################
##Model with beta=(-0.8,-0.8,-0.5) in Section 2.3
####################################################

#dimension
p=3

#sample size
n=100000


#set seed
set.seed(1)


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
print(ALs)
print(bL)
print(SigA)
print(muL)

#simulate sample from PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,4)

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)


#simulated sample:
prop=samp1$samp3


#data on square root scale
z=sqrt(prop)

#a_c value 
acut=1000


########################
## h for Score2 
#######################



#h function
h=matrix(1,n,1)
for (j in 1:p)
{
	h=h*z[,j]
}


hist(h,nclass=100)
#ecdf(h)(5e-11)


#calculate quantiles
quantile(h,probs = seq(0, 1, 0.1))
quantile(h,probs = seq(0, 0.05, 0.01))


########################
## h for Score1 
#######################


indh=matrix(0,n,1)
h2=h
for (j in 1:n)
{
	indh[j]=1
	zmin=z[j,1]
	
	for (i in 2:p)
	{
		zmin_prev=zmin
		zmin=min(zmin,z[j,i])
		if (zmin_prev > zmin){indh[j]=i}
	}
	zmin_prev=zmin
	zmin=min(zmin,acut)
	if (zmin_prev > zmin){indh[j]=0}
	h2[j]=zmin
	
}
h=h2


hist(h,nclass=100)
#ecdf(h)(0.01)

#calculate quantiles
quantile(h,probs = seq(0, 1, 0.1))
quantile(h,probs = seq(0, 0.05, 0.01))










