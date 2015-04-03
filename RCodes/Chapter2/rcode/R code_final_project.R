
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
#Define pdf of MGEV(0,1,xi)
Mgev_pdf <- function(t,xi)
{
  if(xi==0)
    return((1/t^2)*exp(-1/t))
  else 
    #return(exp(-(1+xi*log(t))^(-1/xi))*((1/t)*((1+xi*log(t))^(-1 - (1/xi)))))
    denom = t * ((1+xi*log(t))^(1+1/xi))
    num = (1+xi*log(t))^(-1/xi)
    num = exp(-num)
    return(num/denom)
  
}

#Define survival function of MGEV(0,1,xi)
Mgev_surv <- function(t,xi)
{
  if(xi==0)
    return(1-exp(-1/t))
  else 
    return(1-exp(-(1+xi*log(t))^(-1/xi)))
  
}
#Making sure the log likelihood is always positive
# mylog <-function(x){
#   if(x>0) return(log(x))
#   else return(-10000)
# }
# #Define the joint likelihood
# L_mod <- function(theta1,theta2,y1,y2,phi,delta1,delta2,n)
# {
#   logl <- 0
#   for(i in 1:n)
#   {
#     
#     logl <- logl + delta1[i]*delta2[i]*log(dCphi_du1du2(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi)*Mgev_pdf(y1[i],theta1)*Mgev_pdf(y2[i],theta2))
#                  + delta1[i]*(1-delta2[i])*log((dCphi_du1(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi))*(Mgev_pdf(y1[i],theta1)))
#                  + delta2[i]*(1-delta1[i])*log((dCphi_du2(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi))*(Mgev_pdf(y2[i],theta2)))   
#                  + (1-delta1[i])*(1-delta2[i])*log(-1+Cphi(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi)+Mgev_surv(y1[i],theta1)+Mgev_surv(y2[i],theta2))
#   }
#   
#   return(logl)
#   
# }
# L <- function(theta1,theta2,y1,y2,phi,delta1,delta2,n)
# {
#   logl <- 0
#   for(i in 1:n)
#   {
#     
#     logl <- logl + delta1[i]*delta2[i]*log(dCphi_du1du2(Mgev_surv(y1[i],theta1),Mgev_surv(y2[i],theta2),phi)*Mgev_pdf(y1[i],theta1)*Mgev_pdf(y2[i],theta2))
#     + delta1[i]*(1-delta2[i])*log(dCphi_du1(Mgev_surv(y1[i],theta1),Mgev_surv(y2[i],theta2),phi)*(Mgev_pdf(y1[i],theta1)))
#     + delta2[i]*(1-delta1[i])*log(dCphi_du2(Mgev_surv(y1[i],theta1),Mgev_surv(y2[i],theta2),phi)*(Mgev_pdf(y2[i],theta2)))   
#     + (1-delta1[i])*(1-delta2[i])*log(Cphi(Mgev_surv(y1[i],theta1),Mgev_surv(y2[i],theta2),phi))
#   }
#   
#   return(logl)
#   
# }
L <- function(theta1,theta2,y1,y2,phi,delta1,delta2,n)
{

    
    logl <- delta1*delta2*log(dCphi_du1du2(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi)*Mgev_pdf(y1,theta1)*Mgev_pdf(y2,theta2))
    + delta1*(1-delta2)*log(dCphi_du1(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi)*(Mgev_pdf(y1,theta1)))
    + delta2*(1-delta1)*log(dCphi_du2(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi)*(Mgev_pdf(y2,theta2)))   
    + (1-delta1)*(1-delta2)*log(Cphi(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi))

  
  return(sum(logl))
  
}
# L2 <- function(theta1,theta2,y1,y2,phi,delta1,delta2,n)
# {
#   logl <- 1
#   for(i in 1:n)
#   {
#     
#     logl <- logl * ((dCphi_du1du2(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi)*Mgev_pdf(y1[i],theta1)*Mgev_pdf(y2[i],theta2))^(delta1[i]*delta2[i]))*
#     (((1-dCphi_du1(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi))*(-Mgev_pdf(y1[i],theta1)))^(delta1[i]*(1-delta2[i])))*
#     (((1-dCphi_du2(1-Mgev_surv(y1[i],theta1),1-Mgev_surv(y2[i],theta2),phi))*(-Mgev_pdf(y2[i],theta2)))^(delta2[i]*(1-delta1[i])))*  
#     ((Cphi(Mgev_surv(y1[i],theta1),Mgev_surv(y2[i],theta2),phi))^((1-delta1[i])*(1-delta2[i])))
#   }
#   
#   return(logl)
#   
# }
#Load the simulated data set from Desktop
data <-load("C:/Users/Reseearcher Mode/Desktop/mydata.RData")
y1 <- tc[,1]
y2 <- tc[,2]
delta1 <- delta[,1]
delta2 <- delta[,2]
n <- length(y1)
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
  #################################################
#Define B(phi,phi1)
# B <-function(phi,phi1,PARA,y1,y2,delta1,delta2,n)
# {
#   s <- 0
#   N = nrow(PARA)
#   for(i in 1:N)
#   {
#     num = exp(L(PARA[i,1],PARA[i,2],y1,y2,phi,delta1,delta2,n))
#     den = exp(L(PARA[i,1],PARA[i,2],y1,y2,phi1,delta1,delta2,n))
#     s = s + (num/den)
#   }
#   ratio = s/N
#    return(ratio)
# }
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

phi_hat <-optimize(f = B ,phi1=1,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n,lower = 1,upper = 5, maximum = TRUE)
phi_hat <-optimize(f = B ,phi1=0.5,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n,lower = 0.1,upper = 5, maximum = TRUE)


phi<-seq(0.1,7,0.1)
Lvalues <- unlist(lapply(phi, L, theta1=0.55, theta2=0.55, y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n))

plot(phi, Lvalues)

