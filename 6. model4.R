

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
ALs
bL
beta0



#simulate sample from the tGaussian model:
samp3=rtGaussian(n,p,muL,SigA)


#the response variable prop is the true proportions samp3 
prop=samp3




####################################################################3
##Score1ac plus standard errors
##########################################################3

#a_c for h function:
acut=0.1

#calculate scoring estimate:
estimator=estimator1(prop,acut,1)
estimate1=estimator$estimator1
print(estimate1)

#estimate of W matrix
W_est=estimator$W_est

#standard errors for Score1ac
std1=estimator1SE(prop,acut,estimate1,W_est,1)
print(std1)

#################################################################
##Score2ac plus standard errors
##############################################################

#a_c for h function:
acut=2e-07

#calculate scoring estimate:
estimator=estimator2(prop,acut,1)
estimate2=estimator$estimator2
print(estimate2)

#estimate of W matrix
W_est=estimator$W_est

#standard errors for Score2ac
std2=estimator2SE(prop,acut,estimate2,W_est,1)
print(std2)



#################################################################
##Score2
##############################################################

#a_c for h function (any large value):
acut=10

#calculate scoring estimate:
estimator=estimator2(prop,acut,1)
estimate3=estimator$estimator2
print(estimate3)


####################################################################3
##Score1
##########################################################3

#a_c for h function (any large value):
acut=10

#calculate scoring estimate:
estimator=estimator1(prop,acut,1)
estimate4=estimator$estimator1
print(estimate4)


#combine all estimates together
estimates=t(cbind(t(estimate1),t(std1),t(estimate2),t(std2),t(estimate3),t(estimate4)))

#output parameter estimates
write.table(t(estimates),"estimates.csv",sep=",",row.names=FALSE,col.names=FALSE)



