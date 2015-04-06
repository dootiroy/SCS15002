rm(list=ls())
library(coda)
library(lme4)
source("DRS_priors.R")
source("DRS_pdf.R")
source("DRS_MH_functions.R")
source("DRS_estfns.R")
# set.seed(100)

# Read the data and extract data
test <- read.csv("DRS.csv")
y1 <- test[,6][test$treat==1]
y2 <- test[,6][test$treat==0]
delta1 <- test[,7][test$treat==1]
delta2 <- test[,7][test$treat==0]
n <- length(y1)
# X <- test[,5][test$treat==1]
# X[X==2] <- 0
X <- test[,4][test$treat==1]
DRS <- list(X=X, y1=y1, y2=y2, phi=-999, delta1=delta1, delta2=delta2, n=n)

# Initialize skeletal points #####
skel.points <- c(0.2,0.5,1,3)
# Make a list of all MCMC samples ####
PARA <- list()

for(i in seq_along(skel.points)) {
  ## testing MLE 
  # theta = (sig1, xi1, sig2, xi2, b01, b11, b02, b12)
  # Assuming initial value of phi = 3
  theta0 <- c(1, -0.01, 1, 0, 1.5, 0, 1.5, 0)
  res_mle <- optim(theta0, L_wrap2, , control = list(fnscale=-1,maxit=10000), x=X, y1=y1, 
                   y2=y2, phi=skel.points[i], delta1=delta1, delta2=delta2, n=n, hessian = TRUE)
  
  ## Getting MCMC samples of theta1 and theta2 for fixed phi1
  
  # phi here is phi1
  phi <- skel.points[i]
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
  ##############################
  burnin <- 4000
  #MCMC samples after burn-in
  PARA <- para[-(1:burnin),]
  ##### Thinning to reduce autocorrelation ###
  samp <- matrix(0,1200,8)
  g <- 5
  samp <- PARA[1:g==g,]
  ###########################################
  PARA[[i]] <- samp
}

H.list <- list()

# Initialized the H matrix list which evaluates the log likelihood for every row of each MCMC chains ###
for(k in seq_along(skel.points)){
  H.list[[k]] <- matrix(nrow=nrow(PARA[[1]]), ncol = length(skel.points))
}

# populate the H list of matrices
for(i in 1:nrow(PARA[[1]])){
  for(j in seq_along(skel.points)){
    for(k in seq_along(skel.points)){
      H.list[[j]][i,k] <- L(PARA[[j]][i,3] + PARA[[j]][i,4]*X,PARA[[j]][i,1],
                            PARA[[j]][i,2],PARA[[j]][i,7] + PARA[[j]][i,8]*X,
                            PARA[[j]][i,5],PARA[[j]][i,6],y1,y2,skel.points[k],
                            delta1,delta2,n)
    }
  }
}

# estimate eta's / ri's with r1 = 1
r_init <- c(5,4,1)
r_est <- optim(r_init, etafn, control = list(fnscale = -1), H=H.list)
r <- c(1,r_est$par)

# estimate final phi
phi_est <- optimize(B_fin, lower = 0.2, upper = 4, maximum = TRUE,
                    H = H.list, eta = r, dlist = DRS, PARA = PARA, skel.points = skel.points)
