###### generate random variables(simulation.r)
library(survival)
library(msm)
####Install new package for reaing csv
library(gdata)
library(gplots)
############### READING DATA ##################################################
dat<-read.csv(file="C:/Users/Reseearcher Mode/Desktop/Application_Part/bcancer.csv",head=TRUE)
n<-nrow(dat)
t<-dat$time
delta<-dat$vital
X<-as.matrix(dat[c(5,2,3,7,6)])

fn<-function(para){	
	#use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
        #define a function,return the minus of the log likelihood
	beta<-para[-c(1,2)]
	lambda<-para[1]
	alpha<-para[2]
	

	
	s<-0

	S<-exp(-(t/lambda)^alpha)
  	f<-exp(-(t/lambda)^alpha)*(alpha/lambda)*((t/lambda)^(alpha-1))
	logf<-log(f)		
		  	
	for (i in 1:n){
		if(delta[i]==0)
		s<-s+exp(t(X[i,])%*%beta)*(1-S[i])
		else s<-s+exp(t(X[i,])%*%beta)*(1-S[i])-t(X[i,])%*%beta - logf[i]
		}
	return(s)
			
	}


para0<-c(3,0.1,0,-1,0,0,0)
mle<- nlm(fn,para0,hessian=T)




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
sig_bet<-rep(7,5)
a_lam<-.1 ; b_lam<-2
a_alp<-0.01 ; b_alp<-0.01
a<-10
prop.s<-c(.005,.01,0.005,0.005,rep(.02,3)) # sd of proposal normal for parameters.
para0<-c(5,1,-1.05,0.05,-0.5,-1.08,-0.35) #initial value for parameters.

m<-60000
##  MH algorithm within Gibbs.
burnin <- 5000     #burn-in time

#para: Bayesian estimate
para <-matrix(nrow=m, ncol=7) 
acc.prob <-rep(0,7) 
current.para<-para0


for (tt in 1:m){ 
	
	prop.para<- current.para  
  
            
        #for lambda, lambda^2 uses Inverse Gamma distribution     
        j<-1 
        prop.para[j]<- rnorm(1,current.para[j],prop.s[j])

        #random walk metropolis
        #fn...-log
        loga <- fn(current.para)-fn(prop.para)-fn_prior_lam(current.para)+fn_prior_lam(prop.para)
    
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
        loga <- fn(current.para)-fn(prop.para)-fn_prior_alp(current.para)+fn_prior_alp(prop.para)
    
        u<-runif(1)  
        u1<-log(u)  
        
        if(u1<loga) 
	    { 
            current.para[j]<-prop.para[j]   
            acc.prob[j] <- acc.prob[j]+1 
            }  
             
	
       #for beta's

       for (j in 3:7)
  {
       prop.para<- current.para   
       prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
       #random walk metropolis
       loga <- fn(current.para)-fn(prop.para)-fn_prior_bet(current.para)+fn_prior_bet(prop.para)
    
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
burnin <-10000
#MCMC samples after burn-in
PARA<-para[-(1:burnin),]
samp<-matrix(0,10000,8)
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


para.b;sd.b
#95% HPD interval.
library(lme4)
HPDinterval(samp)

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

M<-read.table(file="C:/Users/Reseearcher Mode/Desktop/CPOmGEV1.txt",head=TRUE)
CG<-M$CPO_GEV
lCPO.i<-log(CPO.i)
d<-numeric(1697)
d=CG-lCPO.i

pdf("C:/Users/Reseearcher Mode/Desktop/mCPOPLOT1.pdf")
plot(d,xlab="Observation",ylab="Difference of log CPO",ylim=c(-1,1),col = ifelse(d < 0,'red','blue'))
dev.off()
























