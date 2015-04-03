library(coda)
library(msm)
library(lme4)

#Define Clayton Copula structure
Cphi <- function(u1,u2,phi)
{
  return(((u1^(-phi))+(u2^(-phi))-1)^(-1/phi))
}
#Define partial derivative structure
dCphi_du1 <- function(u1,u2,phi)
{
  b <- (u1^(-phi)+u2^(-phi)-1)^(-1/phi)
  a <- u1^(-phi-1)*((u1^(-phi)+u2^(-phi)-1)^(-1-(1/phi)))
  return(a)
}
dCphi_du2 <- function(u1,u2,phi)
{
  b <- (u1^(-phi)+u2^(-phi)-1)^(-1/phi)
  a <- u2^(-phi-1)*((u1^(-phi)+u2^(-phi)-1)^(-1-(1/phi)))
  
  return(a)
}
dCphi_du1du2 <- function(u1,u2,phi)
{
  b <- (u1^(-phi)+u2^(-phi)-1)^(-1/phi)
  a <- (1+phi)*((u1*u2)^(-phi-1))*((u1^(-phi)+u2^(-phi)-1)^(-2-(1/phi)))
  return(a)
}
#Define pdf of mgev(mu,sigma,xi)
mgev_pdf <- function(t,mu,sigma,xi)
{
  
  if(xi == 0)
  {
    s = exp(-exp((log(t)-mu)/sigma))
    f = s*exp((log(t) - mu)/sigma)*(1/(sigma*t))
  }
  else
  {
    # check <- rep(exp(mu -sigma/xi),length(t))
    cond <- (t > exp(mu - sigma/xi) & xi > 0) | (t < exp(mu - sigma/xi) & xi < 0)
    S <- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
    r <- (S)*(1/(sigma*t))*((1+xi*(log(t)-mu)/sigma)^(1/xi-1)) #density function
    f = ifelse(cond,r,0)
  }
  
  return(f)
}

#Define survival function of mgev(0,1,xi)
mgev_surv <- function(t,mu,sigma,xi) {
  if(xi != 0){
    cond <- (t > exp(mu - sigma/xi) & xi > 0) | (t < exp(mu - sigma/xi) & xi < 0)
    s1 <- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
    s <- ifelse(cond,s1,1)
  }
  else{
    s = exp(-exp((log(t)-mu)/sigma))
  }
  return(s)
}

L <- function(mu1,sigma1, xi1, mu2,sigma2,xi2,y1,y2,phi,delta1,delta2,n){
  
  
  logl <- delta1*delta2*log(dCphi_du1du2(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi)*mgev_pdf(y1,mu1,sigma1,xi1)*mgev_pdf(y2,mu2,sigma2,xi2))
  + delta1*(1-delta2)*log(dCphi_du1(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi)*(mgev_pdf(y1,mu1,sigma1,xi1)))
  + delta2*(1-delta1)*log(dCphi_du2(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi)*(mgev_pdf(y2,mu2,sigma2,xi2)))   
  + (1-delta1)*(1-delta2)*log(Cphi(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi))
  
  
  return(sum(logl))
  
}


# wrapper function for L for optim - MLE of mu's/sigma's/xi's
# theta will be vector of mu's/sigma's/xi's
# L_wrap1 <- function(theta, y1,y2,phi,delta1,delta2,n) {
#   mu1 <- theta[1]
#   sigma1 <- theta[2]
#   xi1 <- theta[3]
#   mu2 <- theta[4]
#   sigma2 <- theta[5]
#   xi2 <- theta[6]
#   L(mu1,sigma1, xi1, mu2,sigma2,xi2,y1,y2,phi,delta1,delta2,n)
# }
# 
# # testing MLE 
# theta0 <- c(5, 2, 0.1, 5, 2, 0.1)
# res_mle <- optim(theta0, L_wrap1, , control = list(fnscale=-1), y1=y1, y2=y2, phi=5, delta1=delta1, delta2=delta2, n=n)

# Wrapper to incorporate beta's
# theta is vector of beta's/sigma's/xi's
L_wrap2 <- function(theta, x, y1,y2,phi,delta1,delta2,n) {

  sigma1 <- theta[1]
  xi1 <- theta[2]
  sigma2 <- theta[3]
  xi2 <- theta[4]
  beta01 <- theta[5]
  beta11 <- theta[6]
  beta02 <- theta[7]
  beta12 <- theta[8]
  # mu's are now vectors, changing with each subject
  mu1 <- beta01 + beta11*x
  mu2 <- beta02 + beta12*x
  L(mu1,sigma1, xi1, mu2,sigma2,xi2,y1,y2,phi,delta1,delta2,n)
}

# testing MLE 
theta0 <- c(2, 0.1, 2, 0.1, -4, -.4, -4, .4)
res_mle <- optim(theta0, L_wrap2, , control = list(fnscale=-1), x=X, y1=y1, y2=y2, phi=5, delta1=delta1, delta2=delta2, n=n)

# #Load the simulated data set from Desktop
test <-read.csv("C:/Desktop stuff/Copula GEV/DRS.csv")

y1 <- test[,6][test$treat==1]
y2 <- test[,6][test$treat==0]
delta1 <- test[,7][test$treat==1]
delta2 <- test[,7][test$treat==0]
n <- length(y1)
X <- test[,5][test$treat==1]
X[X==2] <- 0


#Define prior distributions
# para0<-c(3,1,.32,3,1,.32,y1,y2,delta1,delta2,3,n)
# mle<- nlm(L,para0,hessian=T)
# mle$estimate

fn_prior_sig1<-function(a_sig, b_sig, sigma1){
  #define a function,return the minus of the log likelihood
  sigma1_sq <- sigma1^2
  
  return((a_sig+1)*log(sigma1_sq)+b_sig/sigma1_sq-log(sigma1))
  
}
fn_prior_sig2<-function(a_sig, b_sig, sigma2){
  #define a function,return the minus of the log likelihood
  sigma2_sq <- sigma2^2
  
  return((a_sig+1)*log(sigma2_sq)+b_sig/sigma2_sq-log(sigma2))
  
}

fn_prior_beta01<-function(beta01,sig_beta01){
  #define a function,return the minus of the log likelihood
  return(beta01^2/sig_beta01^2)
}

fn_prior_beta11<-function(beta11,sig_beta11){
  #define a function,return the minus of the log likelihood
  return(beta11^2/sig_beta11^2)
}
fn_prior_beta02<-function(beta02,sig_beta02){
  #define a function,return the minus of the log likelihood
  return(beta02^2/sig_beta02^2)
}
fn_prior_beta12<-function(beta12,sig_beta12){
  #define a function,return the minus of the log likelihood
  return(beta12^2/sig_beta12^2)
}
#Getting 1000 MCMC samples of theta1 and theta2 for fixed phi1

#phi here is phi1
phi <- 1
phi1 <- phi
 
#Choose hyperparameters

#variances of beta01, beta11, beta02, beta12
sig_beta01 <- 10
sig_beta11 <- 10
sig_beta02 <- 10
sig_beta12 <- 10

#for sigma1 and sigma2 inverse gamma distribution
a_sig <- .01 ; b_sig <- 2
a <- 1

#the prior of xi is uniform(-a,a)
propsd.mu1 <- 0.1
propsd.mu2 <- 0.1
propsd.sigma1 <- 0.1
propsd.sigma2 <- 0.1
propsd.xi1 <- 0.1
propsd.xi2 <- 0.1
propsd.beta01 <- 0.1
propsd.beta11 <- 0.1
propsd.beta02 <- 0.1
propsd.beta12 <- 0.1


# prop.s <- c(.01,.01,0.02,0.01,0.01,0.01,0.01,0.01) # sd of proposal normal for parameters.

#Calculate the initial estimate of each of the parameters except phi and put those in current para
#initial value for parameters.
# para0.mu1 <- mle$estimate 
# para0.mu2 <- 
para0.sigma1 <- res_mle$par[1]
para0.sigma2 <- res_mle$par[3]
para0.xi1 <- res_mle$par[2]
para0.xi2 <- res_mle$par[4]
para0.beta01 <- res_mle$par[5]
para0.beta11 <- res_mle$par[6]
para0.beta02 <- res_mle$par[7]
para0.beta12 <- res_mle$par[8]



#para: Bayesian estimate

# para <-matrix(nrow=m, ncol=8) 
#Define accidental probability
# acc.mu1 <- 0
# acc.mu2 <- 0
acc.sigma1 <- 0
acc.sigma2 <- 0
acc.xi1 <- 0
acc.xi2 <- 0
acc.beta01 <- 0
acc.beta11 <- 0
acc.beta02 <- 0
acc.beta12 <- 0

# current.para<-para0

  
  current.xi1 <- para0.xi1
  current.xi2 <- para0.xi2
  current.sigma1 <- para0.sigma1
  current.sigma2 <- para0.sigma2
  current.beta01 <- para0.beta01
  current.beta11 <- para0.beta11
  current.beta02 <- para0.beta02
  current.beta12 <- para0.beta12
  current.mu1 <-  para0.beta01 + para0.beta11*X
  current.mu2 <-  para0.beta02 + para0.beta12*X

# check updates
c1 <- c2 <- c3 <- c4 <- c5 <- c6 <- c7 <- c8 <- 0

##  MH algorithm within Gibbs.
#Choose total number of simulations
  m<-100
  para <- matrix(nrow=m, ncol=8) 
#Choose burn-in time
  burnin <- 10     
  
for (tt in 1:m){ 
  
#   prop.para <- current.para  
  

  prop.xi1 <- current.xi1
  prop.xi2 <- current.xi2
  prop.sigma1 <- current.sigma1
  prop.sigma2 <- current.sigma2
  prop.beta01 <- current.beta01
  prop.beta11 <- current.beta11
  prop.beta02 <- current.beta02
  prop.beta12 <- current.beta12
  prop.mu1 <- current.beta01 + current.beta11*X
  prop.mu2 <- current.beta02 + current.beta12*X
  
  #for sigma1, sigma1^2 uses Inverse Gamma distribution     
  
  if(current.xi1 > 0)
  {
    lb <- (max(current.mu1) - log(min(y1)))*current.xi1
    prop.sigma1 <- rtnorm(1,current.sigma1,propsd.sigma1 ,lower =lb, upper = Inf) 
  }
  
  if(current.xi1 < 0)
  {
    ub <- (min(current.mu1) - log(max(y1)))*current.xi1
    prop.sigma1 <- rtnorm(1,current.sigma1 ,propsd.sigma1,lower = -Inf, upper = ub) 
  }
  #random walk metropolis
  #fn...-log
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          fn_prior_sig1(a_sig,b_sig,current.sigma1) + fn_prior_sig1(a_sig,b_sig,prop.sigma1)
  
  u<-runif(1)  
  u2<-log(u)  
  if(u2<loga) { 
    c1 <- c1 + 1
    current.sigma1 <- prop.sigma1   
    acc.sigma1 <- acc.sigma1 +1 
  }  
  
  
  #for xi1      

  lb <- min(-1, (current.sigma1/(min(current.mu1)-log(max(y1)))))
  ub <- max(1, (current.sigma1/(max(current.mu1)-log(min(y1)))))
  prop.xi1<- rtnorm(1,current.xi1,propsd.xi1,lower=lb, upper=ub ) 
  
  # quasi-random walk metropolis
  #fn...-log,the jumping distribution is not symmetric.
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) + 
          dnorm(current.xi1, prop.xi1, propsd.xi1,log=TRUE) - dnorm(prop.xi1, current.xi1, propsd.xi1,log=TRUE)
  
  log_acc_prob<-NULL
  log_acc_prob<-loga
  
  u<-runif(1)  
  u3<-log(u)  
  if(u3<log_acc_prob) { 
    current.xi1 <- prop.xi1  
    acc.xi1 <- acc.xi1 + 1 
  }  
  
  #for sigma2, sigma2^2 uses Inverse Gamma distribution     
  
  if(current.xi2 > 0)
  {
    lb <- (max(current.mu2) - log(min(y2)))*current.xi2
    prop.sigma2 <- rtnorm(1,current.sigma2,propsd.sigma2 ,lower =lb, upper = Inf) 
  }
  
  if(current.xi2 < 0)
  {
    ub <- (min(current.mu2) - log(max(y2)))*current.xi2
    prop.sigma2 <- rtnorm(1,current.sigma2 ,propsd.sigma2,lower = -Inf, upper = ub) 
  }
  #random walk metropolis
  #fn...-log
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          fn_prior_sig2(a_sig,b_sig,current.sigma2) + fn_prior_sig2(a_sig,b_sig,prop.sigma2)
  
  u<-runif(1)  
  u2<-log(u)  
  if(u2<loga) { 
    current.sigma2 <- prop.sigma2   
    acc.sigma2 <- acc.sigma2 +1 
  }  
  
  
  #for xi2      
  
  lb <- min(-1, (current.sigma2/(min(current.mu2)-log(max(y2)))))
  ub <- max(1, (current.sigma2/(max(current.mu2)-log(min(y2)))))
  prop.xi2<- rtnorm(1,current.xi2,propsd.xi2,lower=lb, upper=ub ) 
  
  # quasi-random walk metropolis
  #fn...-log,the jumping distribution is not symmetric.
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) + 
          dnorm(current.xi2, prop.xi2, propsd.xi2,log=TRUE) - dnorm(prop.xi2, current.xi2, propsd.xi2,log=TRUE)
  
  log_acc_prob<-NULL
  log_acc_prob<-loga
  
  u<-runif(1)  
  u3<-log(u)  
  if(u3<log_acc_prob) { 
    current.xi2 <- prop.xi2  
    acc.xi2 <- acc.xi2 + 1 
  }  
  
  #for beta01
#     prop.beta01 <- current.beta01  
    if(current.xi1 > 0){
      ub <- log(min(y1)) + current.sigma1/current.xi1 - current.beta11*max(X)
      prop.beta01 <- rtnorm(1,current.beta01,propsd.beta01,lower = -Inf, upper = ub  ) 
    }
    if(current.xi1 < 0){
      lb <- log(max(y1)) + current.sigma1/current.xi1 - current.beta11*min(X)
      prop.beta01 <- rtnorm(1,current.beta01,propsd.beta01,lower = lb, upper = Inf  ) 
    }    
    #random walk metropolis
    loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
            L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
            fn_prior_beta01(current.beta01,sig_beta01) + fn_prior_beta01(prop.beta01,sig_beta01)
    
    u<-runif(1)  
    u4<-log(u)  
    if(u4<loga) { 
      current.beta01 <- prop.beta01  
      current.mu1 <- current.beta01 + current.beta11*X
      acc.beta01 <- acc.beta01 + 1 
    }    
      
  #for beta11
#   prop.beta11 <- current.beta11  
  if(current.xi1 > 0){
    ub <- (log(min(y1)) + current.sigma1/current.xi1 - current.beta01)/max(X)
    prop.beta11 <- rtnorm(1,current.beta11,propsd.beta11,lower = -Inf, upper = ub  ) 
  }
  if(current.xi1 < 0){
    lb <- (log(max(y1)) + current.sigma1/current.xi1 - current.beta01)/min(X)
    prop.beta01 <- rtnorm(1,current.beta01,propsd.beta01,lower = lb, upper = Inf  ) 
  }    
#   prop.beta11 <- rnorm(1,current.beta11,propsd.beta11 )   
  #random walk metropolis
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) -
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) -
          fn_prior_beta11(current.beta11,sig_beta11) + fn_prior_beta11(prop.beta11,sig_beta11)
  
  u<-runif(1)  
  u4<-log(u)  
  if(u4<loga) { 
    current.beta11 <- prop.beta11   
    current.mu1 <- current.beta01 + current.beta11*X
    acc.beta11 <- acc.beta11 + 1 
  }   
  
  
  #for beta02
#   prop.beta02 <- current.beta02  
if(current.xi2 > 0){
  ub <- log(min(y2)) + current.sigma2/current.xi2 - current.beta12*max(X)
  prop.beta02 <- rtnorm(1,current.beta02,propsd.beta02,lower = -Inf, upper = ub  ) 
}
if(current.xi2 < 0){
  lb <- log(max(y2)) + current.sigma2/current.xi2 - current.beta12*min(X)
  prop.beta02 <- rtnorm(1,current.beta02,propsd.beta02,lower = lb, upper = Inf ) 
}    
#   prop.beta02 <- rnorm(1,current.beta02,propsd.beta02 )   
  #random walk metropolis
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) - 
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) -
          fn_prior_beta02(current.beta02,sig_beta02) + fn_prior_beta02(prop.beta02,sig_beta02)
  
  u<-runif(1)  
  u4<-log(u)  
  if(u4<loga) { 
    current.beta02 <- prop.beta02  
    current.mu2 <- current.beta02 + current.beta12*X
    acc.beta02 <- acc.beta02 + 1 
  }   

  #for beta12
#   prop.beta12 <- current.beta12  
#   prop.beta12 <- rnorm(1,current.beta12,propsd.beta12 )   
if(current.xi2 > 0){
  ub <- (log(min(y2)) + current.sigma2/current.xi2 - current.beta02)/max(X)
  prop.beta12 <- rtnorm(1,current.beta12,propsd.beta12,lower = -Inf, upper = ub  ) 
}
if(current.xi2 < 0){
  lb <- (log(max(y2)) + current.sigma2/current.xi2 - current.beta02)/min(X)
  prop.beta12 <- rtnorm(1,current.beta12,propsd.beta12,lower = lb, upper = Inf  ) 
}
  #random walk metropolis
  loga <- L(current.mu1,current.sigma1,current.xi1,current.mu2,current.sigma2,current.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) -
          L(prop.mu1,prop.sigma1,prop.xi1,prop.mu2,prop.sigma2,prop.xi2,y1=y1,y2=y2,phi=phi1,delta1=delta1,delta2=delta2,n=n) -
          fn_prior_beta12(current.beta12,sig_beta12)+fn_prior_beta12(prop.beta12,sig_beta12)
  
  u<-runif(1)  
  u4<-log(u)  
  if(u4<loga) { 
    current.beta12 <- prop.beta12 
    current.mu2 <- current.beta02 + current.beta12*X
    acc.beta12 <- acc.beta12 + 1 
  }   
  
  para[tt,1] <- current.sigma1
  para[tt,2] <- current.xi1
  para[tt,3] <- current.sigma2
  para[tt,4] <- current.xi2
  para[tt,5] <- current.beta01
  para[tt,6] <- current.beta11
  para[tt,7] <- current.beta02
  para[tt,8] <- current.beta12
  
}  

#Checking Diagnostics##############
plot(mcmc(para))
plot(mcmc(PARA))
##############################
burnin <-1000
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
samp<-matrix(0,400,2)
g<-5
samp<-PARA[1:g==g,]

autocorr.plot(samp)
acceptance.rate <- 1 - rejectionRate(mcmc(samp))
acceptance.rate

###########################################
file <- sprintf("C:/Desktop stuff/Copula GEV/phi_8.RData")
save(PARA,file=file)
###########################################
#posterior mean and sd for beta0, beta1, alpha.
para.b<-apply(samp,2,mean) #parameter estimated by Bayes
sd.b<-apply(samp,2,sd) #its SD
#95% HPD interval.
samp1<-mcmc(samp)
HPDinterval(samp1)
#################################################
#Define B(phi,phi1)

B <-function(phi,phi1,PARA,y1,y2,delta1,delta2,n)
{
  s <- 0
  N = nrow(PARA)
  for(i in 1:N)
  {
    num = exp(L(PARA[i,1],PARA[i,2],y1,y2,phi,delta1,delta2,n)-L(PARA[i,1],PARA[i,2],y1,y2,phi1,delta1,delta2,n))
    s = s + num
  }
  ratio = s/N
  return(ratio)
}

phi_hat <-optimize(f = B ,phi1=phi1,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n,lower = 4,upper = 9, maximum = TRUE)
