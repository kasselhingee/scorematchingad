############################################################
##PPI model with beta=(-0.8,-0.8,-0.5)
############################################################

#dimension
p=3

#sample size
n=10000

#set seed
set.seed(1)

#parameters for the PPI model:
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
beta0[3]=-0.5
print(ALs)
print(bL)
print(SigA)
print(muL)

#simulate sample from the PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,4)

#save sample
sampA=samp1$samp3

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)


#plot the marginal distributions (part of Figure 2 in the article)
par(mfrow = c(3,3),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.0)
hist(samp1$samp3[,1],nclass=100,main="",xlab="(a)")
hist(samp1$samp3[,2],nclass=100,main="",xlab="(b)")
hist(samp1$samp3[,3],nclass=100,main="",xlab="(c)")
#dev.off()


############################################################
##PPI model with beta=(-0.5,-0.5,-0.5)
############################################################


#set seed
set.seed(1)

#beta parameter
beta0=matrix(-0.5,p,1)
beta0[3]=-0.5

#simulate sample from the PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,4)

#save sample
sampB=samp1$samp3

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)

#plot the marginal distributions (part of Figure 2 in the article)
hist(samp1$samp3[,1],nclass=100,main="",xlab="(d)")
hist(samp1$samp3[,2],nclass=100,main="",xlab="(e)")
hist(samp1$samp3[,3],nclass=100,main="",xlab="(f)")


############################################################
##PPI model with beta=(0,0,-0.5)
############################################################


#set seed
set.seed(1)

#beta parameter
beta0=matrix(0,p,1)
beta0[3]=-0.5

#simulate sample from the PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,4)

#save sample
sampC=samp1$samp3

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)


#plot the marginal distributions (part of Figure 2 in the article)
hist(samp1$samp3[,1],nclass=100,main="",xlab="(g)")
hist(samp1$samp3[,2],nclass=100,main="",xlab="(h)")
hist(samp1$samp3[,3],nclass=100,main="",xlab="(i)")


############################################################
##PPI model with beta=(0,0,0)
############################################################


#set seed
set.seed(1)

#beta parameter
beta0=matrix(0,p,1)
beta0[3]=0

#simulate sample from the PPI model
samp1=rhybrid(n,p,beta0,ALs,bL,4)

#save sample
sampD=samp1$samp3

#maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler 
#a few times to check that it is an appropriate upper bound. 
maxden=samp1$maxden
print(maxden)


#############################################################
##plot of bivariate distributions 
############################################################

#Figure 1 in the article
par(mfrow = c(2,2),mar=c(4.5,4.5,0.5,0.5),cex.axis=1.0)
plot(sampA[,1],sampA[,2],cex=0.1,xlim=c(0,0.4),ylim=c(0,0.4),xlab="(a)",ylab="")
plot(sampB[,1],sampB[,2],cex=0.1,xlim=c(0,0.4),ylim=c(0,0.4),xlab="(b)",ylab="")
plot(sampC[,1],sampC[,2],cex=0.1,xlim=c(0,0.4),ylim=c(0,0.4),xlab="(c)",ylab="")
plot(sampD[,1],sampD[,2],cex=0.1,xlim=c(0,0.4),ylim=c(0,0.4),xlab="(d)",ylab="")



   
