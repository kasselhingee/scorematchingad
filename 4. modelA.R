###############################################################
##model in Section 2.3 with beta=(-0.8,-0.8,-0.5)
##here a_c is set to 20%
###############################################################


#dimension
p=3

#sample size
n=100


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
samp3=samp1$samp3





####################################################################3
##Score1ac estimator
##########################################################3

#a_c for h function:
acut=1.6e-02


#calculate scoring estimate for full model (only beta[p] fixed at -0.5):
estimator=estimatorall1(samp3,acut,0)
estimate1all=estimator$estimator1
print(estimate1all)

#calculate scoring estimate with beta fixed at beta0:
estimator=estimator1(samp3,acut,1)
estimate1=estimator$estimator1
print(estimate1)


#combine all estimates together
estimates=rbind(estimate1all,estimate1)


#output parameter estimates
write.table(t(estimates),"estimates.csv",sep=",",row.names=FALSE,col.names=FALSE)



