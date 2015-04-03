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


fn1<-function(para){  
  #use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2)]
  lambda<-para[1]
  alpha<-para[2]
  
  
  
  s<-0
  
  S<-exp(-(t/lambda)^alpha)
  f<-exp(-(t/lambda)^alpha)*(alpha/lambda)*((t/lambda)^(alpha-1))
  logf<-log(f)		
  
  s <- sum(exp(X %*% beta)*(1 - S) - (((X %*% beta) + logf)*delta))
  return(s)
  
}


para0<-c(5,0.1,0,0,0,0)
mle<- nlm(fn1,para0,hessian=T)




fn_prior_lam<-function(para){
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2)]
  lambda<-para[1]
  alpha<-para[2]
  lambda2<-lambda^2
  return((a_lam+1)*log(lambda2)+b_lam/lambda2-log(lambda))
  
}


fn_prior_bet<-function(para){
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2)]
  lambda<-para[1]
  alpha<-para[2]
  return(t(beta)%*%diag(sig_bet^(-2))%*%beta)
  
}

fn_prior_alp<-function(para){
  #define a function,return the minus of the log likelihood
  beta<-para[-c(1,2)]
  lambda<-para[1]
  alpha<-para[2]
  return((1-a_alp)*log(alpha)+(alpha/b_alp))
  
}

####### get Bayesian estimation: mean, sd
#use the uniform prior for xi, nonformative prior for beta

#Choosing prior parameters
sig_bet<-rep(5,4)
a_lam<-.1 ; b_lam<-2
a_alp<-0.01 ; b_alp<-0.01
a<-5
prop.s<-c(.01,.01,0.01,0.01,0.01,0.01) # sd of proposal normal for parameters.
para0<-mle$estimate #initial value for parameters.

m<-10000
##  MH algorithm within Gibbs.
burnin <- 5000     #burn-in time

#para: Bayesian estimate
para <-matrix(nrow=m, ncol=6) 
acc.prob <-rep(0,6) 
current.para<-para0


for (tt in 1:m){ 
  
  prop.para<- current.para  
  
  
  #for lambda, lambda^2 uses Inverse Gamma distribution     
  j<-1 
  prop.para[j]<- rnorm(1,current.para[j],prop.s[j])
  
  #random walk metropolis
  #fn...-log
  loga <- fn1(current.para)-fn1(prop.para)-fn_prior_lam(current.para)+fn_prior_lam(prop.para)
  
  u<-runif(1)  
  u1<-log(u)  
  
  if(u1<loga) 
  { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  
  #for alpha      
  j<-2
  prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower=1/a,upper=a ) 
  
  #random walk metropolis
  #fn...-log
  loga <- fn1(current.para)-fn1(prop.para)-fn_prior_alp(current.para)+fn_prior_alp(prop.para)
  
  u<-runif(1)  
  u1<-log(u)  
  
  if(u1<loga) 
  { 
    current.para[j]<-prop.para[j]   
    acc.prob[j] <- acc.prob[j]+1 
  }  
  
  
  #for beta's
  
  for (j in 3:6)
  {
    prop.para<- current.para   
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    loga <- fn1(current.para)-fn1(prop.para)-fn_prior_bet(current.para)+fn_prior_bet(prop.para)
    
    u<-runif(1)  
    u1<-log(u)  
    if(u1<loga) 
    { 
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
samp<-matrix(0,1800,6)
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
  beta<-para[-c(1,2)]
  lambda<-para[1]
  alpha<-para[2]
  
  
  
  s<-0
  
  S<-exp(-(t[i]/lambda)^alpha)
  f<-exp(-(t[i]/lambda)^alpha)*(alpha/lambda)*((t[i]/lambda)^(alpha-1))
  
  
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

lCPO.i<-log(CPO.i)
data <- data.frame(CPO_GEV=lCPO.i)
write.table(data, file = "C:/Users/Dooti/Desktop/SEER/CPOwei.txt")




































