#Getting 1000 MCMC samples of theta1 and theta2 for fixed phi1
library(msm)
a<-0.4
#the prior of xi is uniform(-a,a)
prop.s<-c(0.05,0.05) # sd of proposal normal for parameters.
para0<-c(0.3,0.3) #initial value for parameters.
phi <- 0.08
m<-3000
##  MH algorithm within Gibbs.
burnin <- 1000    #burn-in time

#para: Bayesian estimate
para <-matrix(nrow=m, ncol=2) 
acc.prob <-rep(0,2) 
current.para<-para0
lb <- -1/log(min(y1,y2))

for (tt in 1:m){ 
  
  prop.para<- current.para  
  
  #for xi_1      
  j<-1
  prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower=-a,upper=a ) 
  
  # quasi-random walk metropolis
  #fn...-log,the jumping distribution is not symmetric.
  loga <- -L(current.para[1],current.para[2],y1,y2,phi,delta1,delta2,n)+L(prop.para[1],prop.para[2],y1,y2,phi,delta1,delta2,n)
  log_acc_prob<-NULL
  log_acc_prob<-loga
  +dtnorm(current.para[j], prop.para[j], prop.s[j],-a,lb,log=TRUE) 
  -dtnorm(prop.para[j], current.para[j], prop.s[j],-a,lb,log=TRUE)
  u<-runif(1)  
  u1<-log(u)  
  if(u1<log_acc_prob) { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  #for xi_2      
  j<-2
  prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower=-a,upper=a ) 
  
  # quasi-random walk metropolis
  #fn...-log,the jumping distribution is not symmetric.
  loga <- -L(current.para[1],current.para[2],y1,y2,phi,delta1,delta2,n)+L(prop.para[1],prop.para[2],y1,y2,phi,delta1,delta2,n)
  log_acc_prob<-NULL
  log_acc_prob<-loga
  +dtnorm(current.para[j], prop.para[j], prop.s[j],-a,lb,log=TRUE) 
  -dtnorm(prop.para[j], current.para[j], prop.s[j],-a,lb,log=TRUE)
  u<-runif(1)  
  u1<-log(u)  
  if(u1<log_acc_prob) { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  para[tt,]<-current.para  
  
}

#Checking Diagnostics##############
library(coda)
plot(mcmc(para))
plot(mcmc(PARA))
##############################
burnin <-1000
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
samp<-matrix(0,400,2)
g<-5
samp<-PARA[1:g==g,]
library(coda)
plot(mcmc(para))
plot(mcmc(samp))
autocorr.plot(samp)
acceptance.rate <- 1 - rejectionRate(mcmc(samp))
acceptance.rate

###########################################
file <- sprintf("C:/Users/Reseearcher Mode/Desktop/phi1.RData")
save(PARA,file=file)
##############################
burnin <-10000
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
samp<-matrix(0,8000,3)
g<-5
samp<-PARA[1:g==g,]

###########################################
#posterior mean and sd for beta0, beta1, alpha.
para.b<-apply(samp,2,mean) #parameter estimated by Bayes
sd.b<-apply(samp,2,sd) #its SD
#95% HPD interval.
library(lme4)
samp1<-mcmc(samp)
HPDinterval(samp1)