#########Generating the data of GEV distribution from simulation#############
n<-1000 #sample size
#the true value of the parameter
beta<-c(2,0.6) #true parameter
xi<-0.5
set.seed(123)
#create the covariate
age<-sample(1:1000,n,replace=T)
ages<-scale(age)
X<-cbind(1,ages)
theta<-exp(X%*%beta)
N<-rep(NA,n)
t<-rep(NA,n)
delta<-rep(NA,n)
y<-rep(0,n)
c<-0.8 #change this to adjust the censoring rate.

for(i in 1:n){
 N[i]<-rpois(1,theta[i])
 if(N[i]==0)
      {y[i]<-c
      delta[i]<-0
      }
    else{      
    u<-runif(N[i])
    #set original value: u=0,sigma=1,xi=0.5 for GEV distribution
    xx<- (-1+(-log(1-u))^xi)/xi 
    z<-exp(xx)  #the survival time t[i]<-min(z)
    t[i]<-min(z)
  y[i]<-min(t[i],c)
 delta[i]<- as.numeric(t[i]<=c)

     }
 
}

sum(y==c)/length(y) #censoring rate

data <- data.frame(X=X, y=y,delta=delta)
#####write.table(x, file = "C:/Users/dooti/Desktop/simulation-gev-weibull.txt")

####################################################
###### generate random variables(simulation.r)
######data<-read.table("simulationmin.txt")
######data<-read.table("C:/Users/dooti/Desktop/simulation-gev-weibull.txt")

n<-nrow(data)
y<-data$y
X<-cbind(data$X.1,data$X.2)
delta<-data$delta

#########Drawing Kaplan Meier Curve#############
library(survival)
#use the Kaplan-Meier estimator.

fit1 <- survfit(Surv(y,delta)~1,conf.int=.95,conf.type="log-log",type="kaplan-meier",se.fit=T)
summary(fit1)
plot(fit1,xlab="Survival Time",ylab="Survival function estimate",xlim=c(0.10,0.40))


#plot the survival curve based on the Weibull model with parameter alpha and lambda=1.
fn0<-function(para){
 t<-y
 beta<-para[1:2]
 alpha<-para[3]
 theta<-exp(X%*%beta) 
  S<-exp(-t^alpha)
  f<-exp(-t^alpha)*alpha*(t^(alpha-1))
  return(sum(-theta*(1-S)+delta*log(theta*f)))
}
#########Plotting Survival Curve for model vs. Kaplan Meier###########

library(LearnBayes)  
BETA<-laplace(fn0,c(2,0.6,3))$mode[1:2]
alpha<-laplace(fn0,c(2,0.6,3))$mode[3]
THETA<-exp(X%*%BETA)
surv<-function(xi,theta,y){
t<-y
S<-exp(-t^alpha)
SS<-exp(-theta*(1-S))
return(SS)}
###########################################
#We need to take average ##################
###########################################
YYY<-matrix(0,1000,1000)
for(i in 1:1000){
YYY[i,]<-surv(XI,THETA[i],y)
}
IND<-sample(1:1000,10)
sy<-apply(YYY,2,mean)
pdf("C:/Users/dooti/Desktop/SurvivalCurveByMean2.pdf")
plot(fit1,xlab="Survival Time",ylab="Survival function estimate",xlim=c(0.15,0.40))
points(y,sy,type="p",cex=0.8,pch=19,col="red")
###for(j in 1:10){
###Q<-IND[j]
###points(y,YYY[Q,],type="p",cex=0.1,pch=19,col=j+3)
###}
legend("topright",c("By Kaplan-Meier method","By Mean Proposed model"),lty=1,col=c(1,2))
dev.off()

###### get MLE for data using Weibull distribution#############

#use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
#define a function,return the minus of the log likelihood
fn<-function(para){
	
	beta<-para[c(1,2)]
	alpha<-para[3]
	s<-0
	S<- exp(-y^alpha) #the survival function
	logf<- log(alpha)+(alpha-1)*log(y)-(y^alpha) # the log density function
		

		 
	for (i in 1:n){
		if(delta[i]==0)
		s<-s+exp(t(X[i,])%*%beta)*(1-S[i])
		else s<-s+exp(t(X[i,])%*%beta)*(1-S[i])-t(X[i,])%*%beta - logf[i]
		}
	return(s)
			
	}



##  MH algorithm
burnin <-3000     #burn-in time
prop.s<-c(1.75,0.2,0.001) # sd of proposal normal for beta0, beta1, alpha.
m<-10000
para0<-c(2.7,.7,1.5) #initial value for beta0, beta1, alpha.
#para:beta0,beta1,alpha
para <-matrix(nrow=m, ncol=3) 
acc.prob <-c(0,0,0) 
current.para<-para0

for (t in 1:m){  
    for (j in 1:3){
    prop.para<- current.para   
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    loga <- fn(prop.para)-fn(current.para)
    
    u<-runif(1)  
    u<-log(u)  
    if(u<loga) { 
            current.para<-prop.para   
            acc.prob[j] <- acc.prob[j]+1 
            }    
    }                   
    para[t,]<-current.para  
} 

###########################################
#posterior mean and sd for beta0, beta1, alpha.
para.b<-apply(para[(burnin+1):m,],2,mean) #parameter estimated by Bayes
sd.b<-apply(para[(burnin+1):m,],2,sd) #its SD
#MLE posterior mean and sd for beta0, beta1, alpha.
para.mle<-laplace(fn0,c(2,0.6,3))$mode
sd.mle<-sqrt(diag(laplace(fn0,c(2,0.6,3))$var))
#################################################

#95% HPD interval.
library(lme4)
HPDinterval(para)

##############################
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
##############################
#Calculation for CPO
#likelihood for each observation
fi<-function(para,i){	
	beta<-para[c(1,2)]
	alpha<-para[3]
	
	S<- exp(-y[i]^alpha) #the survival function
	f<- alpha*((y[i])^(alpha-1))*exp(-(y[i]^alpha)) # the log density function		 
	if(delta[i]==0){s<-exp(-exp(t(X[i,])%*%beta)*(1-S))}
	else {s<-exp(-exp(t(X[i,])%*%beta)*(1-S))+exp(t(X[i,])%*%beta)*f}
	return(s)			
	}
#CPO estimator
MM<-length(PARA[,1])
CPO.i<-numeric(n)
FF<-numeric(MM)
for(j in 1:n){
for(k in 1:MM){
FF[k]<-fi(PARA[k,],j)}
CPO.i[j]<-1/mean(1/FF)
}
LPML<-sum(log(CPO.i))

###############
#DIC ##########
MM<-length(PARA[,1])
FF<-numeric(n)
d.bar<-numeric(MM)
for(k in 1:MM){
for(j in 1:n){
FF[j]<-log(fi(PARA[k,],j))
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

