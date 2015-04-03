setwd("~/Desktop/Research/RCodes/Chapter2")
rm(list=ls())
library(coda)
library(lme4)
source("DRS_priors.R")
source("DRS_pdf.R")
source("DRS_MH_functions.R")
source("DRS_estfns.R")
set.seed(100)
test <-read.csv("rcode/DRS.csv")

y1 <- test[,6][test$treat==1]
y2 <- test[,6][test$treat==0]
delta1 <- test[,7][test$treat==1]
delta2 <- test[,7][test$treat==0]
n <- length(y1)
# X <- test[,5][test$treat==1]
# X[X==2] <- 0
X <- test[,4][test$treat==1]

## testing MLE 
# theta = (sig1, xi1, sig2, xi2, b01, b11, b02, b12)
# Assuming initial value of phi = 3
theta0 <- c(1, -0.01, 1, 0, 1.5, 0, 1.5, 0)
res_mle <- optim(theta0, L_wrap2, , control = list(fnscale=-1,maxit=10000), x=X, y1=y1, 
                 y2=y2, phi=2.140715, delta1=delta1, delta2=delta2, n=n, hessian = TRUE)

## Getting MCMC samples of theta1 and theta2 for fixed phi1

# phi here is phi1
phi <- 2.140715
phi1 <- phi

## Choose hyperparameters
# variances of b01, b11, b02, b12
sig_beta <- list()
sig_beta$b01 <- 10
sig_beta$b11 <- 10
sig_beta$b02 <- 10
sig_beta$b12 <- 10

# for sigma1 and sigma2 inverse gamma distribution
a_sig <- .01 ; b_sig <- 2

#the prior of xi is uniform(-a,a)
a <- 1

# sd for proposed densities (truncated normal)
propsd <- list()
propsd$mu1 <- 0.20
propsd$mu2 <- 0.20
propsd$sigma1 <- 0.20
propsd$sigma2 <- 0.20
propsd$xi1 <- 0.14
propsd$xi2 <- 0.14
propsd$b01 <- 0.2
propsd$b11 <- 0.04
propsd$b02 <- 0.2
propsd$b12 <- 0.03

# Calculate the initial estimate of each of the parameters except phi and put those in current para
para0 <- list()
para0$sigma1 <- res_mle$par[1]
para0$sigma2 <- res_mle$par[3]
para0$xi1 <- res_mle$par[2]
para0$xi2 <- res_mle$par[4]
para0$b01 <- res_mle$par[5]
para0$b11 <- res_mle$par[6]
para0$b02 <- res_mle$par[7]
para0$b12 <- res_mle$par[8]

# Run the MH MCMC and get samples
m <- 10000
DRS <- list(X=X, y1=y1, y2=y2, phi=phi, delta1=delta1, delta2=delta2, n=n)
MH.list <- do_MCMC(m, para0, DRS, propsd, sig_beta, a_sig, b_sig)

#Checking Diagnostics##############
para <- MH.list$para
plot(mcmc(para))
##############################
burnin <- 4000
#MCMC samples after burn-in
PARA <- para[-(1:burnin),]
plot(mcmc(PARA))

##### Thinning to reduce autocorrelation ###
samp<-matrix(0,1200,8)
g<-5
samp<-PARA[1:g==g,]

autocorr.plot(para)
autocorr.plot(PARA)
autocorr.plot(samp)
acceptance.rate <- 1 - rejectionRate(mcmc(samp))
acceptance.rate <- 1 - rejectionRate(mcmc(para))
acceptance.rate

###########################################
# file <- sprintf("C:/Desktop stuff/Copula GEV/phi_8.RData")
file <- sprintf("RESULTS/phi_hatfin.RData")
save(samp,file=file)
###########################################
#posterior mean and sd for beta0, beta1, alpha.
para.b<-apply(samp,2,mean) #parameter estimated by Bayes
sd.b<-apply(samp,2,sd) #its SD
#95% HPD interval.
samp1<-mcmc(samp)
HPDinterval(samp1)
#################################################
# Load the MCMC samples for each of the 4 individual phi1 ###
load("RESULTS/phi_pt2.RData")
PARA1 <- samp
load("RESULTS/phi_half.RData")
PARA2 <- samp
load("RESULTS/phi_1.RData")
PARA3 <- samp
load("RESULTS/phi_3.RData")
PARA4 <- samp
load("RESULTS/phi_1pt8.RData")
PARA5 <- samp
load("RESULTS/phi_2pt4.RData")
PARA6 <- samp
### Estimating phi_hat by maximizing B_phi_phi1 ####
phi1 <- 3
phi_hat <- optimize(f = B ,phi1=phi1,PARA=PARA4,y1=y1,y2=y2,delta1=delta1,
                    delta2=delta2,n=n,lower = 0.2,upper = 4, maximum = TRUE)

### Plot B_phi_phi1 with different values of phi1. Show 4 curves and show that they achieve maximum at different points.###

phi_grid <- seq(0.2,3,0.1)
B_vect <- Vectorize(B,"phi")
y_vect1 <- B_vect(phi_grid,phi1=0.2,PARA=PARA1,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n)
y_vect2 <- B_vect(phi_grid,phi1=0.5,PARA=PARA2,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n)
y_vect3 <- B_vect(phi_grid,phi1=1,PARA=PARA3,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n)
y_vect4 <- B_vect(phi_grid,phi1=3,PARA=PARA4,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n)
y_vect5 <- B_vect(phi_grid,phi1=1.8,PARA=PARA5,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n)
y_vect6 <- B_vect(phi_grid,phi1=2.4,PARA=PARA6,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n)

pdf(file = "plot1")
plot(phi_grid, y_vect1, type="l",col="red", ylim = c(0,70), xlim =c(0,3),xlab = expression(phi),
     ylab = expression('B'[phi][','][phi[1]]))
# abline(h =1.38324, col="red")
# abline(h =1, col="blue")
# abline(h =2.366155, col="black")
lines(phi_grid, y_vect2, col="blue")
lines(phi_grid, y_vect3, col="black")
lines(phi_grid, y_vect5, col="orange")
lines(phi_grid, y_vect6, col="green")
#lines(phi_grid, y_vect4, col="green")
#legend("topright",c(expression(phi*' = '*'0.2'),expression(phi*' = '*'0.5'),
#                   expression(phi*' = '*'1'),expression(phi*' = '*'3')),col=c("red","blue","black","green"))
legend("topright",c(expression(phi[1]*' = '*'0.2'),expression(phi[1]*' = '*'0.5'),
expression(phi[1]*' = '*'1'),expression(phi[1]*' = '*'1.8'),expression(phi[1]*' = '*'2.4'))
,lty = c(1,1,1,1,1),col= c("red","blue","black","orange","green"))

dev.off()
