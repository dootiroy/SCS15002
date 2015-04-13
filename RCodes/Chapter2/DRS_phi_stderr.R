setwd("~/Desktop/Research/cluster/output")
rm(list=ls())

phi <- matrix(0,ncol = 1, nrow = 100)

for(i in 0:99){
  load(paste("phiest_",i,".RData",sep =""))
  phi[i+1,1] <- phi_est$maximum
}
hist(phi[,1], xlab = expression(phi), col = "grey", main = expression(paste('Histogram of ',phi)))
sd(phi[,1])
quantile(phi[,1],c(0.025,0.975))
summary(phi[,1])
tau <- matrix(0,ncol = 1, nrow = 100)
tau[,1] <- phi[,1]/(phi[,1] + 2)
sd(tau[,1])
quantile(tau[,1],c(0.025,0.975))

