###### generate random variables(simulation.r)
library(survival)
library(msm)
#setwd("~/Documents/cc")
d<-read.table("data.txt",header=TRUE)
#only keep the type "Sigmoid Colon"
dat<-subset(d,county==2 & time!=0)
dat$county<-NULL
dat$yeardiag<-NULL

n<-nrow(dat)
t<-dat$time
delta<-dat$vitalstatus
X<-cbind(1,dat[,-c(1,2)])





######################without covariates ########################
#plot the survival curve based on the model.
fn0<-function(para){
	
	theta<-para[4]
	mu<-para[1]
	sigma<-para[2]
	xi<-para[3]
	s<-0
	
	if(xi==0){  
		S<-1-exp(-exp(-((log(t)-mu)/sigma)))
		logf<- log(1-S)-(log(t)-mu)/sigma-log(sigma*t)
		 }	
	else{ 
		S<- 1-exp(-(1+xi*(log(t)-mu)/sigma)^(-1/xi)) #the survival function
	logf<- log(1-S)-(1/xi+1)*log(1+xi*(log(t)-mu)/sigma)-log(sigma*t) # the log density function
		
		 }		
		  
	for (i in 1:n){
		if(delta[i]==0)
		s<-s+theta*(1-S[i])
		else s<-s+theta*(1-S[i])-log(theta) - logf[i]
		}

	return(s)
			
	}


###### MLE methods(no covariate)
ss<- nlm(fn0,c(1,1,0.03,1),hessian=T)
ss$estimate

surv<-function(para,t){
	theta<-para[4]
	mu<-para[1]
	sigma<-para[2]
	xi<-para[3]
return(exp(-theta*exp(-(1+xi*(log(t)-mu)/sigma)^(-1/xi)) ))

}

surv_t<-surv(ss$estimate,t)
library(survival)
#use the Kaplan-Meier estimator.
fit <- survfit(Surv(t,delta)~1)
summary(fit)

pdf("plot_km_non.pdf")
plot(fit,xlab="Survival Time",ylab="Survival function estimate",ylim=c(0.7,1))
points(t,surv_t,type="p",cex=0.05,pch=19)
dev.off()






####### get Bayesian estimation: mean, sd
#use the uniform prior for xi, nonformative prior for beta
#how to choose???????
sig_mu<-10
sig_bet<-10
a_sig<-1 ; b_sig<-1


fn_prior_mu<-function(para){
    #define a function,return the minus of the log likelihood
	theta<-para[4]
	mu<-para[1]
	sigma<-para[2]
	xi<-para[3]
	return(mu^2/(2*sig_mu^2))
		
	}
	

fn_prior_sig<-function(para){
    #define a function,return the minus of the log likelihood
	theta<-para[4]
	mu<-para[1]
	sigma<-para[2]
	sigma2<-sigma^2
	xi<-para[3]
	return((a_sig+1)*log(sigma2)+b_sig/sigma2-log(sigma))
		
	}
		
fn_prior_bet<-function(para){
    #define a function,return the minus of the log likelihood
	theta<-para[4]
	mu<-para[1]
	sigma<-para[2]
	xi<-para[3]
	return(theta^2/(2*sig_bet^2))
		
	}



##  MH algorithm within Gibbs.
burnin <- 2000     #burn-in time

prop.s<-rep(0.01,4) # sd of proposal normal for parameters.
m<-10000
para0<-ss$estimate #initial value for parameters.
#quite sensative for initial value

#para: Bayesian estimate
para <-matrix(nrow=m, ncol=4) 
acc.prob <-rep(0,4) 
current.para<-para0


#set.seed(123)

for (tt in 1:m){ 
	
	prop.para<- current.para  
	
	#for mu
	j<-1 
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    #fn0...-log
    loga <- fn0(current.para)-fn0(prop.para)+fn_prior_mu(current.para)-fn_prior_mu(prop.para)
    
    u<-runif(1)  
    u<-log(u)  
    if(u<loga) { 
            current.para<-prop.para   
            acc.prob[j] <- acc.prob[j]+1 
            }  
            
                      
     #for sigma       
     j<-2 
     prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    #fn...-log
    loga <- fn0(current.para)-fn0(prop.para)+fn_prior_sig(current.para)-fn_prior_sig(prop.para)
    
    u<-runif(1)  
    u<-log(u)  
    if(u<loga) { 
            current.para<-prop.para   
            acc.prob[j] <- acc.prob[j]+1 
            }  
            
            
      #for xi            
      j<-3
    prop.para<- current.para   
    prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    loga <- fn0(current.para)-fn0(prop.para)
    
    u<-runif(1)  
    u<-log(u)  
    if(u<loga) { 
            current.para<-prop.para   
            acc.prob[j] <- acc.prob[j]+1 
            }    

	 #for theta
	 j<-4
      prop.para[j]<- rnorm(1,current.para[j],prop.s[j] )   
    #random walk metropolis
    #fn...-log
    loga <- fn0(current.para)-fn0(prop.para)+fn_prior_bet(current.para)-fn_prior_bet(prop.para)  
    u<-runif(1)  
    u<-log(u)  
    if(u<loga) { 
            current.para<-prop.para   
            acc.prob[j] <- acc.prob[j]+1 
            }  
                          
    para[tt,]<-current.para  
    
    
} 



#write.table(para,file="bay_estimation_nocov.txt")
#write.table(acc.prob/m,file="acc.prob_nocov.txt")
#para<-read.table("bay_estimation_nocov.txt",header=T)
print(acc.prob/m)
#the acceptance probability.

#posterior mean and sd for beta0, beta1, xi.
para.b<-apply(para[(burnin+1):m,],2,mean)
sd.b<-apply(para[(burnin+1):m,],2,sd) 

para.b;sd.b
#95% HPD interval.
library(lme4)
HPDinterval(para)

#posterior means for cure rates.
pi_est<-exp(-para[(burnin+1):m,4])


pdf("hist.pdf")
hist(pi_est,main="",xlab="cure rates")
abline(v=exp(-ss$estimate[4]),lty="dotdash")
#0.0067 estimate for cure rates using MLE.
dev.off()


mean(pi_est)
median(pi_est)











########convergence diagnostics test

# convergence diagnostics plot 

erg.mean<-function(x){ # compute ergodic mean 
        n<-length(x)
        result<-cumsum(x)/cumsum(rep(1,n))
  }
m<-nrow(para)
step<-5
idx<-seq(1,m,step)
idx2<-seq(burnin+1,m)

  #trace plot
 plot(idx,para[idx,1],type="l")
 plot(idx,para[idx,2],type="l")
 plot(idx,para[idx,3],type="l")


ergbeta0<-erg.mean(para[,1])
ergbeta02<-erg.mean(para[idx2,1])
ylims0<-range(c(ergbeta0,ergbeta02))   
    
ergbeta1<-erg.mean(para[,2])
ergbeta12<-erg.mean(para[idx2,2])
ylims1<-range(c(ergbeta1,ergbeta12))

ergbeta2<-erg.mean(para[,3])
ergbeta22<-erg.mean(para[idx2,3])
ylims2<-range(c(ergbeta2,ergbeta22))

ergbeta3<-erg.mean(para[,3])
ergbetae32<-erg.mean(para[idx2,3])
ylimse<-range(c(ergbeta2,ergbeta22))


#Ergodic Mean Plot.
plot(idx, ergbeta0[idx], type='l', ylab='Values of beta0', xlab='Iterations', main=' Ergodic Mean Plot of beta0', ylim=ylims0)
lines(idx2, ergbeta02[idx2-burnin], col=2, lty=2)

plot(idx, ergbeta1[idx], type='l', ylab='Values of beta1',xlab='Iterations', main=' Ergodic Mean Plot of beta1', ylim=ylims1)
lines(idx2, ergbeta12[idx2-burnin], col=2, lty=2)

plot(idx, ergbeta2[idx], type='l', ylab='Values of beta1',xlab='Iterations', main=' Ergodic Mean Plot of xi', ylim=ylims2)
lines(idx2, ergbeta22[idx2-burnin], col=2, lty=2)


#Autocorrelations Plot.
index3<-seq(1,m,20) 
index4<-seq(burnin+1,m,step)

par(mfrow=c(1,3))
acf(para[index4,1], main='Autocorrelations Plot for beta0', sub='(Thin = 20 iterations)')
acf(para[index4,2], main='Autocorrelations Plot for beta1', sub='(Thin = 20 iterations)')
acf(para[index4,3], main='Autocorrelations Plot for xi', sub='(Thin = 20 iterations)')





