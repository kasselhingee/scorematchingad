
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

library(MCMCpack)
library(sirt)

#simulate sample from Dirichlet distribution
dirprop=rdirichlet(n, alpha)

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


####################################################################
##MLDir 
##########################################################

estimate7=dirichlet.mle(dirfit,eps=10^(-15),convcrit=1e-10,maxit=10000,progress=FALSE)$alpha-1
print(estimate7)

####################################################################
##MLDirMult 
##########################################################

library(dirmult)
estimate8=dirmult(x,epsilon=10^(-10),trace=FALSE)$gamma-1


####################################################################
##MomDir
##########################################################

estimate6=dirichmom(dirfit)-1



####################################################################
##Score2
##########################################################


#a_c for h function (any large value):
acut=10


#calculate scoring estimate:
estimate3=estimator2_dir(dirfit,acut)
print(estimate3)



####################################################################
##Score2ac
##########################################################


#a_c for h function:
acut=0.001


#calculate scoring estimate:
estimate2=estimator2_dir(dirfit,acut)
print(estimate2)


####################################################################
##Score1
##########################################################


#a_c for h function (any large value):
acut=10


#calculate scoring estimate:
estimate4=estimator1_dir(dirfit,acut)
print(estimate4)


####################################################################
##Score1ac
##########################################################

#a_c for h function:
acut=0.01


#calculate scoring estimate:
estimate1=estimator1_dir(dirfit,acut)
print(estimate1)


#output estimates
estimates2=t(cbind(t(estimate7),t(estimate6),t(estimate3),t(estimate2),t(estimate4),t(estimate1)))
write.table(t(estimates2),"estimates.csv",sep=",",row.names=FALSE,col.names=FALSE)




