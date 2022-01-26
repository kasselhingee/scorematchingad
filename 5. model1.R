

#dimension
p=5

#sample size
n=92

#set seed
set.seed(1)

#parameters for the PPI model:
ALs=matrix(0,p-1,p-1)
bL=matrix(0,p-1,1)
ALs[1,1]= -127480.0929
ALs[1,2]= 14068.39057
ALs[1,3]= 1782.261826 
ALs[1,4]=  -240.076568 
ALs[2,1]= 14068.3906 
ALs[2,2]= -8191.17253 
ALs[2,3]=  -8.002680
ALs[2,4]= 374.693979
ALs[3,1]=1782.2618 
ALs[3,2]= -8.00268
ALs[3,3]= -46.638659 
ALs[3,4]= 9.027633
ALs[4,1]= -240.0766 
ALs[4,2]=  374.69398  
ALs[4,3]=  9.027633
ALs[4,4]= -39.208915
beta0=matrix(0,p,1)
beta0[1]=-0.80
beta0[2]=-0.85
beta0[3]=0
beta0[4]=-0.2
beta0[5]=0
ALs
bL
beta0


#simulate sample from the PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,0)

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)

#simulated sample:
samp3=samp1$samp3


#simulate sample from the multinomial PPI model:
ni=matrix(2000,n,1)
ni=as.vector(ni)
x=matrix(0,n,p)
for (j in 1:n)
{
	x[j,]=rmultinom(1,ni[j],prob=samp3[j,])
}
prop1=x/ni

#the response variable prop is either true proportions samp3 
#or estimated proportions prop1
prop=samp3
#prop=prop1



####################################################################3
##Score1ac plus standard errors
##########################################################3

#a_c for h function:
acut=0.01

#calculate scoring estimate:
estimator=estimator1(prop,acut,0)
estimate1=estimator$estimator1
print(estimate1)

#estimate of W matrix
W_est=estimator$W_est

#standard errors for Score1ac
std1=estimator1SE(prop,acut,estimate1,W_est,0)
print(std1)

#################################################################
##Score2ac plus standard errors
##############################################################

#a_c for h function:
acut=0.0001

#calculate scoring estimate:
estimator=estimator2(prop,acut,0)
estimate2=estimator$estimator2
print(estimate2)

#estimate of W matrix
W_est=estimator$W_est

#standard errors for Score2ac
std2=estimator2SE(prop,acut,estimate2,W_est,0)
print(std2)



#################################################################
##Score2
##############################################################

#a_c for h function (any large value):
acut=10

#calculate scoring estimate:
estimator=estimator2(prop,acut,0)
estimate3=estimator$estimator2
print(estimate3)


####################################################################3
##Score1
##########################################################3

#a_c for h function (any large value):
acut=10

#calculate scoring estimate:
estimator=estimator1(prop,acut,0)
estimate4=estimator$estimator1
print(estimate4)


#combine all estimates together
estimates=t(cbind(t(estimate1),t(std1),t(estimate2),t(std2),t(estimate3),t(estimate4)))

#output parameter estimates
write.table(t(estimates),"estimates.csv",sep=",",row.names=FALSE,col.names=FALSE)



####################################################################
##ScoreMult p=3 only with all elements of beta equal
##########################################################3


#mult=multestimator(x)
#print(mult)

#combine all estimates together
#estimates=t(cbind(t(estimate1),t(std1),t(estimate2),t(std2),t(mult),t(estimate3),t(estimate4)))

#output parameter estimates
#write.table(t(estimates),"estimates.csv",sep=",",row.names=FALSE,col.names=FALSE)


