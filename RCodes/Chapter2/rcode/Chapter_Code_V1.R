library(copula)

#Generate Bivariate Data using Clayton Copula from mGEV distribution#
cop2 <- claytonCopula(5, dim = 2)
v <- rCopula(2000, cop2)

pairs(v)
v1 <- 1-v

#Take inverse survival transformation on v to get mgev(0,1,xi) samples.
xi <- 0.3

#set original value: xi=0.4
if(xi != 0) {
  t =  exp((-1 + (-log(v1))^xi)/xi)    #the survival time
} else {
  t = -log(v1)
}

# Criteria to reject samples with t in "bad" range
if(xi > 0){
  crit <- apply((t > exp(-1/xi)), 1, all)
} else if (xi < 0){
  crit <- apply((t < exp(-1/xi)), 1, all)
} else {
  crit <- rep(TRUE,2000)
}

all(crit)

# Select first 1000 TRUEs 
t <- head(t[crit,],1000)

pairs(t)


#Get your censored y and delta
c <- 2 #change this to adjust the censoring rate.
tc <- t
tc[tc > c] <- c
delta <- (t < c)*1
sum(tc[,1]==c)/length(tc[,1]) #censoring rate of y1
sum(tc[,2]==c)/length(tc[,2]) #censoring rate of y2
pairs(tc)

save(tc,delta,file="C:/Desktop stuff/Copula GEV/mydata.RData")


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
#Define pdf of mgev(0,1,xi)
mgev_pdf <- function(t,xi)
{
  
  if(xi == 0)
  {
    f = exp(-t)
  }
  else
  {
    check <- rep(exp(-1/xi),length(t))
    cond <- ((t > check & xi > 0) | (t < check & xi < 0))
    r = (1/t)*((1 + xi*log(t))^(1/xi - 1))*(exp(-(1 + xi*log(t))^(1/xi)))
    f = ifelse(cond,r,0)
  }
  
  return(f)
}

#Define survival function of mgev(0,1,xi)
mgev_surv <- function(t,xi)
{
  if(xi != 0)
  {
    v <- 1 + xi*log(t)
    v[v < 0] <- 0
    s =  exp(-(v^(1/xi)))    #the survival time
  }
  else
  {
    s = exp(-t)
  }
  return(s)
}

L <- function(xi1,xi2,y1,y2,phi,delta1,delta2,n)
{
  
  
  logl <- delta1*delta2*log(dCphi_du1du2(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi)*mgev_pdf(y1,xi1)*mgev_pdf(y2,xi2))
  + delta1*(1-delta2)*log(dCphi_du1(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi)*(mgev_pdf(y1,xi1)))
  + delta2*(1-delta1)*log(dCphi_du2(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi)*(mgev_pdf(y2,xi2)))   
  + (1-delta1)*(1-delta2)*log(Cphi(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi))
  
  
  return(sum(logl))
  
}


# #Load the simulated data set from Desktop
data <-load("C:/Desktop stuff/Copula GEV/mydata.RData")

y1 <- tc[,1]
y2 <- tc[,2]
delta1 <- delta[,1]
delta2 <- delta[,2]
n <- length(y1)

#Getting 1000 MCMC samples of theta1 and theta2 for fixed phi1
library(msm)
a<-0.7
#the prior of xi is uniform(-a,a)
prop.s<-c(0.05,0.05) # sd of proposal normal for parameters.
para0<-c(0.1,0.1) #initial value for parameters.
#phi here is phi1
phi <- 5.211642
phi1 <- phi
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
library(lme4)
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


#Part 2 

load("C:/Desktop stuff/Copula GEV/phi_2.RData")
data1 <- PARA
load("C:/Desktop stuff/Copula GEV/phi_4.RData")
data2 <- PARA
load("C:/Desktop stuff/Copula GEV/phi_6.RData")
data3 <- PARA
load("C:/Desktop stuff/Copula GEV/phi_8.RData")
data4 <- PARA

phi.skel <- c(2,4,6,8)
PARA <- list(data1,data2,data3,data4)
N <- length(data1[,1])
J <- length(phi.skel)
H <- matrix(nrow=N,ncol=J)
for (i in 1:N){
  for (j in 1:J){
    H[i,j] <- L(PARA[[j]][i,1],PARA[[j]][i,2],y1,y2,phi.skel[j],delta1,delta2)
  }
}
H <- exp(H)
bigf <- function(r1,H){ #r1 is a vector of length J
#   r <- c(1,r1) r is just r1
  r <- r1
  J <- dim(H)[2]
  t1 <- sum(log(t(H)*exp(r)))
  t2 <- sum(log(colSums(t(H)*exp(r))))*J
  return(t1-t2)
}

init <- c(.1,.1,.1,.1)
r <- optim(init,bigf,method="L-BFGS-B",lower=0.001,control=list(fnscale=-1),H=H)
r.est <- r$par
rfin <- exp(-(r.est - r.est[1]))



mat1 <- matrix(nrow=N,ncol=J)
biglist <- list()
for (j in 1:J){
  biglist[[j]] <- mat1
}

for (l in 1:N){
  for (j in 1:J){
    for(i in 1:J){
      biglist[[i]][l,j] <- N*exp(L(PARA[[j]][l,1],PARA[[j]][l,2],y1,y2,phi.skel[i],delta1,delta2))/rfin[i]
    }
  }
}

biglistmat <- matrix(data=0,nrow=N,ncol=J)
for (i in 1:J){
  biglistmat <- biglistmat + biglist[[i]]
}

B.big <-function(phi,PARA,y1,y2,delta1,delta2,biglistmat)
{
  s = 0
  N <- dim(biglistmat)[1]
  J <- dim(biglistmat)[2]
  for (i in 1:N){
    for (j in 1:J){
      temp <- exp(L(PARA[[j]][i,1],PARA[[j]][i,2],y1,y2,phi,delta1,delta2))/biglistmat[i,j]
      s <- s + temp
    }
  }
  return(s)
}

phi.est.final <- optimize(f = B.big ,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,biglistmat=biglistmat,lower = 0.1,upper = 10, maximum = TRUE)

