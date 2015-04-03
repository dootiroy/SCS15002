###### generate random variables(simulation.r)
library(survival)
library(msm)
####Install new package for reaing csv
library(gdata)
library(gplots)
############### READING DATA ##################################################
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

# fn <- function(para){  
#   #use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
#   #define a function,return the minus of the log likelihood
#   theta<-para[4]
#   mu<-para[1]
#   sigma<-para[2]
#   xi<-para[3]
#   
#   s<-0
#   if(xi==0){  
#     S <-exp(-exp((log(t)-mu)/sigma))
#     logf <- log(S)+(log(t)-mu)/sigma-log(sigma*t)
#   }  
#   else { 
#     cond <- (t > exp(mu - sigma/xi) & xi > 0) | (t < exp(mu - sigma/xi) & xi < 0)
#     S <- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
#     f <- (S)*(1/(sigma*t))*((1+xi*(log(t)-mu)/sigma)^(1/xi-1)) #density function
#     # logf <- log(f)                                # the log density function		
#     logf <- ifelse(cond, log(f), Inf)
#   }		
#   
#   s <- sum(-theta*(1 - S) + ((log(theta) + logf)*delta))
#   return(s)
#   
# }
# 
# 
# para0<-c(1,1,0.5,0.1)
# mle<- nlm(fn,para0,hessian=T)
# 
# mle$estimate
fn0<-function(para){
  
  theta<-para[4]
  mu<-para[1]
  sigma<-para[2]
  xi<-para[3]
  s<-0
  
  if(xi==0){  
    S<-exp(-exp((log(t)-mu)/sigma))
    logf<- log(S)+(log(t)-mu)/sigma-log(sigma*t)
  }	
  else{ 
    S<- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
    logf<- log(S)+(1/xi-1)*log(1+xi*(log(t)-mu)/sigma)-log(sigma*t) # the log density function
    
  }		
  
  for (i in 1:n){
    if(delta[i]==0)
      s<-s+theta*(1-S[i])
    else s<-s+theta*(1-S[i])-log(theta) - logf[i]
  }
  
  return(s)
  
}

############### MLE(without covariates) ##############################
ss<- nlm(fn0,c(1,1,0.5,1),hessian=T)
ss$estimate

surv<-function(para,t){
  theta<-para[4]
  mu<-para[1]
  sigma<-para[2]
  xi<-para[3]
  return(exp(-theta*(1-exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi))) ))
  
}

surv_t<-surv(ss$estimate,t)

library(survival)
##################use the Kaplan-Meier estimator.###############################

fit1 <- survfit(Surv(t,delta)~1,conf.int=.95,conf.type="log-log",type="kaplan-meier",se.fit=TRUE)
summary(fit1)
plot(fit1,xlab="Survival Time in years",ylab="Survival function estimate",ylim=c(0,1),xlim=c(0,120))

###################### GEV FIT #################################################

# surv<-function(para,t){
#   theta<-para[4]
#   mu<-para[1]
#   sigma<-para[2]
#   xi<-para[3]
#   return(exp(-theta*(1-exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi))) ))
#   
# }
# 
# 
# surv_t<-surv(mle$estimate,t)

#################  WEIBULL FIT  ####################################
fnwei<-function(para){
  lambda<-para[1]
  alpha<-para[2]
  theta<-para[3]
  S<-exp(-(t/lambda)^alpha)
  f<-exp(-(t/lambda)^alpha)*(alpha/lambda)*((t/lambda)^(alpha-1))
  return(sum(-theta*(1-S)+delta*log(theta*f)))
}
library(LearnBayes)  
LAMBDA<-laplace(fnwei,c(1,0.6,1))$mode[1]
ALPHA<-laplace(fnwei,c(1,0.6,1))$mode[2]
THETA<-laplace(fnwei,c(1,0.6,1))$mode[3]


survwei<-function(lambda,alpha,theta,t){
  S<-exp(-(t/lambda)^alpha)
  SS<-exp(-theta*(1-S))
  return(SS)}

surv_w<-survwei(LAMBDA,ALPHA,THETA,t)
#####################################################################
#########Plotting Survival Curve for model vs. Kaplan Meier##########

pdf("C:/Users/Dooti/Desktop/SEER/KMCurveApply.pdf")
plot(fit1,xlab="Survival Time (in months)",ylab="Survival function estimate",ylim=c(0,1), xlim=c(0,150))
abline(h=0.12,col=4,lty=2,)
text(6, 0.2, "12%", col = "red") 
points(t,surv_t,type="p",cex=0.8,pch=19,col="red")
points(t,surv_w,type="p",cex=0.8,pch=19,col="blue")
legend("topright",c("By Kaplan-Meier method","Fitted mGEV model","Fitted Weibull Model" ),lty=1,col=c(1,2,4))
dev.off()
#####################################################################
# 
# ########################################
# #MCMC ##################################
# #Define prior for mu xi sigma beta#######################
# 
# fn_prior_mu<-function(para){
#   #define a function,return the minus of the log likelihood
#   beta<-para[-c(1,2,3)]
#   mu<-para[1]
#   sigma<-para[2]
#   xi<-para[3]
#   return(mu^2/(2*sig_mu^2))
#   
# }
# 
# fn_prior_sig<-function(para){
#   #define a function,return the minus of the log likelihood
#   beta<-para[-c(1,2,3)]
#   mu<-para[1]
#   sigma<-para[2]
#   sigma2<-sigma^2
#   xi<-para[3]
#   return((a_sig+1)*log(sigma2)+b_sig/sigma2-log(sigma))
#   
# }
# 
# 
# fn_prior_bet<-function(para){
#   #define a function,return the minus of the log likelihood
#   beta<-para[-c(1,2,3)]
#   mu<-para[1]
#   sigma<-para[2]
#   xi<-para[3]
#   return(t(beta)%*%diag(sig_bet^(-2))%*%beta)
#   
# }
# 
# 
# 
# 
# ####### get Bayesian estimation: mean, sd
# #use the uniform prior for xi, nonformative prior for beta
# #how to choose???????
# sig_mu<-10
# sig_bet<-rep(10,5)
# a_sig<-.001 ; b_sig<-.001
# a<-0.6
# #the prior of xi is uniform(-a,a)
# prop.s<-c(.01,.01,0.05,rep(.01,5)) # sd of proposal normal for parameters.
# para0<-c(3.04785156,0.01171875,0.36035156,-1.38964844,0.31035156,10.61035156,2.90722656,-0.76855469) #initial value for parameters.
# 
# m<-400
# ##  MH algorithm within Gibbs.
# burnin <- 50     #burn-in time
# 
# #para: Bayesian estimate
# para <-matrix(nrow=m, ncol=8) 
# acc.prob <-rep(0,8) 
# current.para<-para0
# 
# 
# for (tt in 1:m){ 
#   
#   prop.para<- current.para  
#   
#   #for mu
#   j<-1 
#   prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
#   #random walk metropolis
#   #fn...-log
#   loga <- fn(current.para)-fn(prop.para)+fn_prior_mu(current.para)-fn_prior_mu(prop.para)
#   
#   u<-runif(1)  
#   u<-log(u)  
#   if(u<loga) { 
#     current.para[j]<-prop.para[j]   
#     acc.prob[j] <- acc.prob[j]+1 
#   }  
#   
#   
#   
#   #for sigma, sigma^2 uses Inverse Gamma distribution     
#   j<-2 
#   prop.para[j]<- rnorm(1,current.para[j],prop.s[j])
#   #random walk metropolis
#   #fn...-log
#   loga <- fn(current.para)-fn(prop.para)+fn_prior_sig(current.para)-fn_prior_sig(prop.para)
#   
#   u<-runif(1)  
#   u<-log(u)  
#   if(u<loga) { 
#     current.para[j]<-prop.para[j]   
#     acc.prob[j] <- acc.prob[j]+1 
#   }  
#   
#   
#   #for xi      
#   j<-3
#   prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower=-a,upper=a ) 
#   
#   # quasi-random walk metropolis
#   #fn...-log,the jumping distribution is not symmetric.
#   loga <- fn(current.para)-fn(prop.para)
#   log_acc_prob<-NULL
#   log_acc_prob<-loga+pnorm((a-current.para[j])/prop.s[j],log.p=TRUE)-pnorm((a-prop.para[j])/prop.s[j],log.p=TRUE)
#   u<-runif(1)  
#   u<-log(u)  
#   if(u<log_acc_prob) { 
#     current.para[j]<-prop.para[j]   
#     acc.prob[j] <- acc.prob[j]+1 
#   }  
#   
#   
#   #for beta's
#   for (j in 4:8){
#     prop.para<- current.para   
#     prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
#     #random walk metropolis
#     loga <- fn(current.para)-fn(prop.para)+fn_prior_bet(current.para)-fn_prior_bet(prop.para)
#     
#     u<-runif(1)  
#     u<-log(u)  
#     if(u<loga) { 
#       current.para[j]<-prop.para[j]   
#       acc.prob[j] <- acc.prob[j]+1 
#     }    
#   }    
#   
#   
#   para[tt,]<-current.para  
#   
#   
# } 









































