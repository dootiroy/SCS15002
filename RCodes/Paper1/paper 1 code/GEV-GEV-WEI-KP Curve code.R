n<-2000 #sample size
#the true value of the parameter
beta<-c(2,0.6) #true parameter
xi<-0.5
set.seed(222)
#create the covariate
age<-sample(1:2000,n,replace=T)
ages<-scale(age)
X<-cbind(1,ages)
theta<-exp(X%*%beta)
N<-rep(NA,n)
t<-rep(NA,n)
delta<-rep(NA,n)
y<-rep(0,n)
c<-0.5 #change this to adjust the censoring rate.

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

n<-nrow(data)
y<-data$y
X<-cbind(data$X.1,data$X.2)
delta<-data$delta

#########Drawing Kaplan Meier Curve#############(
library(survival)
#use the Kaplan-Meier estimator.

fit1 <- survfit(Surv(y,delta)~1,conf.int=.95,conf.type="log-log",type="kaplan-meier",se.fit=T)
summary(fit1)
plot(fit1,xlab="Survival Time",ylab="Survival function estimate",xlim=c(0.15,0.40))



####Plot the survival curve with GEV Minima distribution######
#plot the survival curve based on the Weibull model with parameter alpha and lambda=1.
fngev<-function(para){
 t<-y
 beta<-para[1:2]
 xi<-para[3]
 theta<-exp(X%*%beta) 
  S<-exp(-(1+xi*log(t))^(1/xi))
  f<-exp(-(1+xi*log(t))^(1/xi))*(1/t)*(1+xi*log(t))^((1/xi)-1)
 return(sum(-theta*(1-S)+delta*log(theta*f)))
}
#########Plotting Survival Curve for model vs. Kaplan Meier###########
####MLE is calculated by Learn Bayes package###
library(LearnBayes)  
BETAG<-laplace(fngev,c(2,0.6,0.0003))$mode[1:2]
XIG<-laplace(fngev,c(2,0.6,0.0003))$mode[3]
THETAG<-exp(X%*%BETAG)
survgev<-function(xi,theta,y){
t<-y
S<-exp(-(1+xi*log(t))^(1/xi))
SS<-exp(-theta*(1-S))
return(SS)}
###########################################
#We need to take average ##################
###########################################
YYY<-matrix(0,2000,2000)
for(i in 1:2000){
YYY[i,]<-survgev(XIG,THETAG[i],y)
}
IND<-sample(1:2000,10)
sy<-apply(YYY,2,mean)
#plot the survival curve based on the Weibull model with parameter alpha and lambda=1.
fnwei<-function(para){
 t<-y
 beta<-para[1:2]
 alpha<-para[3]
 theta<-exp(X%*%beta) 
  S<-exp(-(t^alpha))
  f<-exp(-(t^alpha))*alpha*(t^(alpha-1))
  return(sum(-theta*(1-S)+delta*log(theta*f)))
}

#########Plotting Survival Curve for model vs. Kaplan Meier###########
####MLE is calculated by Learn Bayes package###
library(LearnBayes)  
BETA<-laplace(fnwei,c(2,0.6,3))$mode[1:2]
ALPHA<-laplace(fnwei,c(2,0.6,3))$mode[3]
THETA<-exp(X%*%BETA)
survwei<-function(alpha,theta,y){
t<-y
S<-exp(-t^alpha)
SS<-exp(-theta*(1-S))
return(SS)}
###########################################
#We need to take average ##################
###########################################
ZZZ<-matrix(0,2000,2000)
for(i in 1:2000){
ZZZ[i,]<-survwei(ALPHA,THETA[i],y)
}
IND1<-sample(1:2000,10)
sz<-apply(ZZZ,2,mean)
pdf("C:/Users/Reseearcher Mode/Desktop/SCGEV-GEV-WEI.pdf")
plot(fit1,xlab="Survival Time",ylab="Survival function estimate",xlim=c(0.15,0.40))
points(y,sy,type="p",cex=0.4,pch=19,col="Blue")
points(y,sz,type="p",cex=0.4,pch=19,col="Red")
legend("topright",c("By Kaplan-Meier method","Fitted mGEV Model","Fitted Weibull Model"),lty=1,col=c(1,4,2))
dev.off()

