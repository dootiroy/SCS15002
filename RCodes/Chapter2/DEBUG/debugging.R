### DEBUGGING scripts ##

load("DEBUG/check.RData")

# data
X <- data.list$X
y1 <- data.list$y1
y2 <- data.list$y2
phi <- data.list$phi
delta1 <- data.list$delta1
delta2 <- data.list$delta2
n <- data.list$n

# parameters
b01 <- proposed$b01
b11 <- proposed$b11
b02 <- proposed$b02
b12 <- proposed$b12
# mu1 <- proposed$mu1
# mu2 <- proposed$mu2
mu1 <- b01 + b11*X
mu2 <- b02 + b12*X
xi1 <- proposed$xi1
xi2 <- proposed$xi2
sigma1 <- proposed$sigma1
sigma2 <- proposed$sigma2
