dat<-read.csv(file="C:/Users/Reseearcher Mode/Desktop/Application_Part/bcancer.csv",head=TRUE)
n<-nrow(dat)
t<-dat$time
delta<-dat$vital
X<-as.matrix(dat[c(5,2,3,7,6)])

fn<-function(para){	
	#use the sum that extends over all possible values of the vector N && (A.1) as the likelihood.
    #define a function,return the minus of the log likelihood
	beta<-para[-c(1,2,3)]
	mu<-para[1]
	sigma<-para[2]
	xi<-para[3]

	s<-0
	if(xi==0){  
		S<-exp(-exp((log(t)-mu)/sigma))
		logf<- log(S)+(log(t)-mu)/sigma-log(sigma*t)
		 }	
	else { 
		S<- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
		f<-(S)*(1/(sigma*t))*((1+xi*(log(t)-mu)/sigma)^(1/xi-1)) #density function
	      logf<- log(f)                                # the log density function		
		 }		
		  

	for (i in 1:n){
		if(delta[i]==0)
		s<-s+exp(t(X[i,])%*%beta)*(1-S[i])
		else s<-s+exp(t(X[i,])%*%beta)*(1-S[i])-t(X[i,])%*%beta - logf[i]
		}	
	return(s)
			
	}


para0<-c(3,1,.32,rep(0,4))
mle<- nlm(fn,para0,hessian=T)


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
sig_bet<-c(5,6,6,6,6)
a_sig<-.01 ; b_sig<-2
a<-1

#the prior of xi is uniform(-a,a)

prop.s<-c(.01,.01,0.02,0.01,0.01,0.01,0.01,0.01) # sd of proposal normal for parameters.
para0<-c(4.3,1.1,0.18,0.03,-0.08,-1.07,-1.06) #initial value for parameters.

m<-60000

##  MH algorithm within Gibbs.

burnin <- 10000     #burn-in time

#para: Bayesian estimate

para <-matrix(nrow=m, ncol=8) 
acc.prob <-rep(0,8) 
current.para<-para0


for (tt in 1:m){ 
	
	prop.para<- current.para  
	
	#for mu
	j<-1 
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
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
     prop.para[j]<- rnorm(1,current.para[j],prop.s[j])
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
      prop.para[j]<- rtnorm(1,current.para[j],prop.s[j],lower=-a,upper=a ) 
   
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
    for (j in 4:8){
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
#MLE posterior mean and sd for beta0, beta1, alpha.
para.mle<-laplace(fnwei,c(2,0.6,3))$mode
sd.mle<-sqrt(diag(laplace(fnwei,c(2,0.6,3))$var))
#################################################

#95% HPD interval.
library(lme4)
HPDinterval(samp)



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

lCPO.i<-log(CPO.i)
data <- data.frame(CPO_GEV=lCPO.i)
write.table(data, file = "C:/Users/Reseearcher Mode/Desktop/CPOmGEV1.txt")






































