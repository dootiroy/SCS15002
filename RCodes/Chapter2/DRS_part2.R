setwd("~/Desktop/Research/RCodes/Chapter2")
rm(list=ls())
library(coda)
library(lme4)
source("DRS_priors.R")
source("DRS_pdf.R")
source("DRS_MH_functions.R")
source("DRS_estfns.R")
set.seed(100)

# Read the data and extract data
test <- read.csv("rcode/DRS.csv")
y1 <- test[,6][test$treat==1]
y2 <- test[,6][test$treat==0]
delta1 <- test[,7][test$treat==1]
delta2 <- test[,7][test$treat==0]
n <- length(y1)
# X <- test[,5][test$treat==1]
# X[X==2] <- 0
X <- test[,4][test$treat==1]
DRS <- list(X=X, y1=y1, y2=y2, phi=-999, delta1=delta1, delta2=delta2, n=n)

# Load the MCMC samples for each of the 4 individual phi1 ###
load("RESULTS/phi_pt2.RData")
PARA1 <- samp
load("RESULTS/phi_half.RData")
PARA2 <- samp
load("RESULTS/phi_1.RData")
PARA3 <- samp
load("RESULTS/phi_3.RData")
PARA4 <- samp

# Initialize skeletal points #####
skel.points <- c(0.2,0.5,1,3)

# Make a list of all MCMC samples ####
PARA <- list(PARA1, PARA2, PARA3, PARA4)
H.list <- list()

# Initialized the H matrix list which evaluates the log likelihood for every row of each MCMC chains ###
for(k in seq_along(skel.points)){
  H.list[[k]] <- matrix(nrow=nrow(PARA1), ncol = length(skel.points))
}

# populate the H list of matrices
for(i in 1:nrow(PARA1)){
  for(j in seq_along(skel.points)){
    for(k in seq_along(skel.points)){
      H.list[[j]][i,k] <- L(PARA[[j]][i,3] + PARA[[j]][i,4]*X,PARA[[j]][i,1],
                            PARA[[j]][i,2],PARA[[j]][i,7] + PARA[[j]][i,8]*X,
                            PARA[[j]][i,5],PARA[[j]][i,6],y1,y2,skel.points[k],
                            delta1,delta2,n)
    }
  }
}
save(list = c("H.list"), file = "Hlist.RData")
# When you have all PARA1-4, dont run the above section again and again. Using save(), 
# save H.list and then just load it

# estimate eta's / ri's with r1 = 1
r_init <- c(5,4,1)
r_est <- optim(r_init, etafn, control = list(fnscale = -1), H=H.list)
r <- c(1,r_est$par)

# estimate final phi
phi_est <- optimize(B_fin, lower = 0.2, upper = 4, maximum = TRUE,
                    H = H.list, eta = r, dlist = DRS, PARA = PARA, skel.points = skel.points)
