
#########Genarting the data from simulation#############
n<-1000 #sample size
#the true value of the parameter
beta<-c(1,0.2) #true parameter
xi<-0.3 #GEV parameter
set.seed(123)
#creat the covariate
age<-sample(1:100,n,replace=T)
ages<-scale(age)
X<-cbind(1,ages)
theta<-exp(X%*%beta)
N<-rep(NA,n)
t<-rep(NA,n)
delta<-rep(NA,n)
y<-rep(0,n)
c<-1.7 #change this to adjust the censoring rate.

for(i in 1:n){
 N[i]<-rpois(1,theta[i])
 if(N[i]==0)
      {y[i]<-c
      delta[i]<-0
      }
    else{      
    u<-runif(N[i])
    #set original value: u=0,sigma=1,xi=0.3 for GEV distribution
    xx<- ((-log(u))^(-xi)-1)/xi 
    z<-exp(xx)  #the survival time t[i]<-min(z)
    t[i]<-min(z)
  y[i]<-min(t[i],c)
 delta[i]<- as.numeric(t[i]<=c)

     }
 
}

sum(y==c)/length(y) #censoring rate

x <- data.frame(X=X, y=y,delta=delta)
write.table(x, file = "C:/Users/dor11002/Desktop/simulationgev.txt")
#We just generated Data set and save in desktop!!!##
####################################################
###### generate random variables(simulation.r)
#data<-read.table("simulationmin.txt")
data<-read.table("C:/Users/dor11002/Desktop/simulationgev.txt")

n<-nrow(data)
y<-data$y
X<-cbind(data$X.1,data$X.2)
delta<-data$delta

#########Drawing Kaplan Meier Curve#############
library(survival)
#use the Kaplan-Meier estimator.

fit1 <- survfit(Surv(y,delta)~1,conf.int=.95,conf.type="log-log",type="kaplan-meier",se.fit=T)
summary(fit1)
plot(fit1,xlab="Survival Time",ylab="Survival function estimate",ylim=c(0.2,1),xlim=c(0.2,1))


#plot the survival curve based on the model.
fn0<-function(para){
 t<-y
 beta<-para[1:2]
 xi<-para[3]
 theta<-exp(X%*%beta) 
  S<-1-exp(-(1+xi*log(t))^(-1/xi))
  f<-exp(-(1+xi*log(t))^(-1/xi))/(t*(1+xi*log(t))^(1+1/xi))
 return(sum(-theta*(1-S)+delta*log(theta*f)))
}
#########Plotting Survival Curve for model vs. Kaplan Meier###########

library(LearnBayes)  
BETA<-laplace(fn0,c(2,0.6,0.3))$mode[1:2]
XI<-laplace(fn0,c(2,0.6,0.3))$mode[3]
THETA<-exp(X%*%BETA)
surv<-function(xi,theta,y){
t<-y
S<-1-exp(-(1+xi*log(t))^(-1/xi))
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
pdf("C:/Users/dor11002/Desktop/SurvivalCurveByMean.pdf")
plot(fit1,xlab="Survival Time",ylab="Survival function estimate",ylim=c(0.2,1), xlim=c(0.2,1))
points(y,sy,type="p",cex=0.8,pch=19,col="blue")
legend("topright",c("By Kaplan-Meier method","By Proposed MGEV model"),lty=1,col=c(1,4))
dev.off()

##################################Bayesian Estimation###############
library(mvtnorm)

###### get MLE for data

#use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
#define a function,return the minus of the log likelihood
fn<-function(para){
	
	beta<-para[c(1,2)]
	xi<-para[3]
	s<-0
	if(xi==0){  
		S<-1-exp(-1/y)
		logf<- -1/y-2*log(y)
		 }
	
	else{ 
		S<- 1-exp(-(1+xi*(log(y)))^(-1/xi)) #the survival function
	logf<- (-(1+xi*(log(y)))^(-1/xi))-log(y)-(1/xi+1)*log(1+xi*(log(y))) # the log density function
		
		 }
		 
	for (i in 1:n){
		if(delta[i]==0)
		s<-s+exp(t(X[i,])%*%beta)*(1-S[i])
		else s<-s+exp(t(X[i,])%*%beta)*(1-S[i])-t(X[i,])%*%beta - logf[i]
		}
	return(s)
			
	}


##  MH algorithm
burnin <-3000     #burn-in time
prop.s<-c(1.75,0.2,0.01) # sd of proposal normal for beta0, beta1, xi.
m<-10000
para0<-c(0,0,0) #initial value for beta0, beta1, xi.
#para:beta0,beta1,xi
para <-matrix(nrow=m, ncol=3) 
acc.prob <-c(0,0,0) 
current.para<-para0

for (t in 1:m){  
    for (j in 1:3){
    prop.para<- current.para   
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    loga <- fn(current.para)-fn(prop.para)
    
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
#posterior mean and sd for beta0, beta1, xi.
para.b<-apply(para[(burnin+1):m,],2,mean) #parameter estimated by Bayes
sd.b<-apply(para[(burnin+1):m,],2,sd) #its SD
#MLE posterior mean and sd for beta0, beta1, xi.
para.mle<-laplace(fn0,c(2,0.6,0.3))$mode
sd.mle<-sqrt(diag(laplace(fn0,c(2,0.6,0.3))$var))
#################################################

#95% HPD interval.
library(coda)
library(lme4)
para1 <- mcmc(para)
HPDinterval(para1)

#posterior means for cure rates.
pi.b <- exp(- exp(X%*%para.b[c(1,2)]))
pi.mle <- exp(- exp(X%*%para.mle[c(1,2)]))

#compare the MLEs and Posterior means for cure rates using boxplots.
pdf("C:/Users/dor11002/Desktop/Boxplots.pdf")
pi.estimate <- cbind(c(pi.mle,pi.b),c(rep(1,n),rep(2,n)))
boxplot(pi.estimate[,1]~pi.estimate[,2],ylab="Cure Rates")
dev.off()
 


##############################
#############################
# ###########

#Calculation for CPO
#likelihood for each observation
fi<-function(para,i){	
	beta<-para[c(1,2)]
	alpha<-para[3]
	s<-0
        if(xi==0){  
		  S<-1-exp(-1/y)
		  logf<- exp(-1/y	)*(y^(-2))	 
		 }
	
	else{ 
	S<- 1-exp(-(1+xi*(log(y[i])))^(-1/xi)) #the survival function
	logf<- (-(1+xi*log(y[i]))^(-1/xi))-log(y[i])-(1/xi+1)*log(1+xi*(log(y[i]))) # the log density function	 
	if(delta[i]==0){s<-exp(t(X[i,])%*%beta)*(1-S[i])}
	else {s<-exp(t(X[i,])%*%beta)*(1-S[i])-t(X[i,])%*%beta - logf[i]}
	return(s)			
	}

}
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
##############################

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
 