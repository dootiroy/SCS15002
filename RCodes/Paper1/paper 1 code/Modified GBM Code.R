gbm<-read.csv("C:/Users/Dooti/Desktop/SEER/GBM.csv")

# removing patients with zero survival months
gbm <- gbm[gbm$survival_months > 0, ]
library(survival)
##Change all the covariates included in the model as binary factors if they were categorical##


##For young adults age: 15-39####
gbm1 <- gbm[gbm$age_at_diagnosis >= 16 & gbm$age_at_diagnosis <= 39,]
gbm2 <- gbm1[gbm1$year_of_diagnosis >= 1970 & gbm1$year_of_diagnosis <=2004,]


y2 <-gbm2$survival_months
delta2 <- gbm2$vstatus

fit3 <- survfit(Surv(y2,delta2)~1,conf.int=.95,conf.type="log-log",type="kaplan-meier",se.fit=T)
summary(fit3)

plot(fit3,xlab="Survival Time  (months)",ylab="Survival function estimate",xlim=c(0,120))
abline(h=0.12,col=4,lty=2,)


n<-nrow(gbm2)
t<-y2
delta<-delta2
X<-as.matrix(gbm2[c(6,23,24,25)])

fn <- function(para){  
  #use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2,3)]
  mu<-para[1]
  sigma<-para[2]
  xi<-para[3]
  
  s<-0
  if(xi==0){  
    S <-exp(-exp((log(t)-mu)/sigma))
    logf <- log(S)+(log(t)-mu)/sigma-log(sigma*t)
  }	
  else { 
    cond <- (t > exp(mu - sigma/xi) & xi > 0) | (t < exp(mu - sigma/xi) & xi < 0)
    S <- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
    f <- (S)*(1/(sigma*t))*((1+xi*(log(t)-mu)/sigma)^(1/xi-1)) #density function
    # logf <- log(f)                                # the log density function		
    logf <- ifelse(cond, log(f), Inf)
  }		
  
  s <- sum(exp(X %*% beta)*(1 - S) - (((X %*% beta) + logf)*delta))
  return(s)
  
}


para0<-c(3,1,.32,rep(0,4))
mle<- nlm(fn,para0,hessian=T)

mle$estimate

fn_prior_mu<-function(para){
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2,3)]
  mu<-para[1]
  sigma<-para[2]
  xi<-para[3]
  return(mu^2/(2*sig_mu^2))
  
}

fn_prior_sig<-function(para){
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2,3)]
  mu<-para[1]
  sigma<-para[2]
  sigma2<-sigma^2
  xi<-para[3]
  return((a_sig+1)*log(sigma2)+b_sig/sigma2-log(sigma))
  
}


fn_prior_bet<-function(para){
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2,3)]
  mu<-para[1]
  sigma<-para[2]
  xi<-para[3]
  return(t(beta)%*%diag(sig_bet^(-2))%*%beta)
  
}



####### get Bayesian estimation: mean, sd
#use the uniform prior for xi, nonformative prior for beta

library(msm)

#Choose hyperparameters

sig_mu<-4
sig_bet<-c(5,6,6,6)
a_sig<-.01 ; b_sig<-2
a<-1

#the prior of xi is uniform(-a,a)

prop.s<-c(.01,.01,0.02,0.01,0.01,0.01,0.01) # sd of proposal normal for parameters.
para0<-mle$estimate #initial value for parameters.

m<-10000

##  MH algorithm within Gibbs.

burnin <- 10     #burn-in time

#para: Bayesian estimate

para <-matrix(nrow=m, ncol=7) 
acc.prob <-rep(0,7) 
current.para<-para0


for (tt in 1:m){ 
  
  prop.para <- current.para  
  
  #for mu
  j<-1 
  if(current.para[3]>0)
  {
    ub <- log(min(t)) + current.para[2]/current.para[3]
    prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower = -Inf, upper = ub) 
  }
  
  if(current.para[3]<0)
  {
    lb <- log(max(t)) + current.para[2]/current.para[3]
    prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower = lb, upper = Inf) 
  }
  
  #random walk metropolis
  #fn...-log
  loga <- fn(current.para)-fn(prop.para)-fn_prior_mu(current.para)+fn_prior_mu(prop.para)
  
  u<-runif(1)  
  u1<-log(u)  
  if(u1<loga) { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  
  
  #for sigma, sigma^2 uses Inverse Gamma distribution     
  j<-2 
  if(current.para[3]>0)
  {
    lb <- (current.para[1] - log(min(t)))*current.para[3]
    prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower =lb, upper = Inf) 
  }
  
  if(current.para[3]<0)
  {
    ub <- (current.para[1] - log(max(t)))*current.para[3]
    prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower = -Inf, upper = ub) 
  }
  #random walk metropolis
  #fn...-log
  loga <- fn(current.para)-fn(prop.para)-fn_prior_sig(current.para)+fn_prior_sig(prop.para)
  
  u<-runif(1)  
  u2<-log(u)  
  if(u2<loga) { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  
  #for xi      
  j<-3
  lb <- min(-1, (current.para[2]/(current.para[1]-log(max(t)))))
  ub <- max(1, (current.para[2]/(current.para[1]-log(min(t)))))
  prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower=lb, upper=ub ) 
  
  # quasi-random walk metropolis
  #fn...-log,the jumping distribution is not symmetric.
  loga <- fn(current.para)-fn(prop.para)
  log_acc_prob<-NULL
  log_acc_prob<-loga
  +dnorm(current.para[j], prop.para[j], prop.s[j],log=TRUE) 
  -dnorm(prop.para[j], current.para[j], prop.s[j],log=TRUE)
  u<-runif(1)  
  u3<-log(u)  
  if(u3<log_acc_prob) { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  
  #for beta's
  for (j in 4:7){
    prop.para<- current.para   
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    loga <- fn(current.para)-fn(prop.para)-fn_prior_bet(current.para)+fn_prior_bet(prop.para)
    
    u<-runif(1)  
    u4<-log(u)  
    if(u4<loga) { 
      current.para[j]<-prop.para[j]   
      acc.prob[j] <- acc.prob[j]+1 
    }    
  }    
  
  
  para[tt,]<-current.para  
  
  
}


library(coda)
plot(mcmc(para))
plot(mcmc(PARA))
##############################
burnin <-1000
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
samp<-matrix(0,1800,7)
g<-5
samp<-PARA[1:g==g,]
library(coda)
plot(mcmc(para))
plot(mcmc(samp))
autocorr.plot(samp)
acceptance.rate <- 1 - rejectionRate(mcmc(samp))
acceptance.rate

###########################################
#posterior mean and sd for beta0, beta1, alpha.
para.b<-apply(samp,2,mean) #parameter estimated by Bayes
sd.b<-apply(samp,2,sd) #its SD
#MLE posterior mean and sd for beta0, beta1, alpha.
para.mle<-laplace(fnwei,c(2,0.6,3))$mode
sd.mle<-sqrt(diag(laplace(fnwei,c(2,0.6,3))$var))
#################################################

#95% HPD interval.
library(lme4)
HPDinterval(mcmc(samp))



#Calculation for CPO
#likelihood for each observation
fi<-function(para,i){
  beta<-para[-c(1,2,3)]
  mu<-para[1]
  sigma<-para[2]
  xi<-para[3]
  
  s<-0
  
  if(xi==0){  
    S<-exp(-exp((log(t[i])-mu)/sigma))
    f<- S*exp((log(t[i])-mu)/sigma)*(1/(sigma*t[i]))
  }	
  else { 
    S<- exp(-(1+xi*(log(t[i])-mu)/sigma)^(1/xi)) #the survival function
    f<-(S)*(1/(sigma*t[i]))*((1+xi*(log(t[i])-mu)/sigma)^(1/xi-1)) #density function
    
  }		
  
  
  s<-exp(-exp(t(X[i,])%*%beta)*(1-S))*((exp(t(X[i,])%*%beta)*f)^delta[i])
  return(s)
  
}
#CPO estimator
MM<-length(samp[,1])
CPO.i<-numeric(n)
FF<-numeric(MM)
for(j in 1:n){
  for(k in 1:MM){
    FF[k]<-fi(samp[k,],j)}
  CPO.i[j]<-1/mean(1/FF)
}
LPML<-sum(log(CPO.i))

###############
#DIC ##########
MM<-length(samp[,1])
FF<-numeric(n)
d.bar<-numeric(MM)
for(k in 1:MM){
  for(j in 1:n){
    FF[j]<-log(fi(samp[k,],j))
  }
  d.bar[k]<-(-2)*sum(FF)
}
D.bar<-mean(d.bar) #firt part of DIC
FF<-numeric(n)

for(j in 1:n){
  FF[j]<-log(fi(para.b,j))
}
D.hat<-(-2)*sum(FF)
DIC<-2*D.bar-D.hat 

LPML
DIC

M<-read.table(file="C:/Users/Dooti/Desktop/SEER/CPOwei.txt",head=TRUE)
CG<-M$CPO_GEV
lCPO.i<-log(CPO.i)
d<-numeric(1725)
d = lCPO.i - CG

pdf("C:/Users/Dooti/Desktop/SEER/mCPOPLOT1.pdf")
plot(d,xlab="Observation",ylab="Difference of log CPO",ylim=c(-1,1),col = ifelse(d < 0,'red','blue'))
dev.off()



































