#Define Clayton Copula structure
Cphi <- function(u1,u2,phi)
{
  return(((u1^(-phi))+(u2^(-phi))-1)^(-1/phi))
}
#Define partial derivative structure
dCphi_du1 <- function(u1,u2,phi)
{
  b <- (u1^(-phi)+u2^(-phi)-1)^(-1/phi)
  a <- u1^(-phi-1)*((u1^(-phi)+u2^(-phi)-1)^(-1-(1/phi)))
  return(a)
}
dCphi_du2 <- function(u1,u2,phi)
{
  b <- (u1^(-phi)+u2^(-phi)-1)^(-1/phi)
  a <- u2^(-phi-1)*((u1^(-phi)+u2^(-phi)-1)^(-1-(1/phi)))
  
  return(a)
}
dCphi_du1du2 <- function(u1,u2,phi)
{
  b <- (u1^(-phi)+u2^(-phi)-1)^(-1/phi)
  a <- (1+phi)*((u1*u2)^(-phi-1))*((u1^(-phi)+u2^(-phi)-1)^(-2-(1/phi)))
  return(a)
}
#Define pdf of MGEV(0,1,xi)
Mgev_pdf <- function(t,xi)
{
  if(xi==0)
    return((1/t^2)*exp(-1/t))
  else 
    #return(exp(-(1+xi*log(t))^(-1/xi))*((1/t)*((1+xi*log(t))^(-1 - (1/xi)))))
    denom = t * ((1+xi*log(t))^(1+1/xi))
  num = (1+xi*log(t))^(-1/xi)
  num = exp(-num)
  return(num/denom)
  
}

#Define survival function of MGEV(0,1,xi)
Mgev_surv <- function(t,xi)
{
  if(xi==0)
    return(1-exp(-1/t))
  else 
    return(1-exp(-(1+xi*log(t))^(-1/xi)))
  
}

L <- function(theta1,theta2,y1,y2,phi,delta1,delta2,n=length(y1))
{
  
  
  logl <- delta1*delta2*log(dCphi_du1du2(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi)*Mgev_pdf(y1,theta1)*Mgev_pdf(y2,theta2))
  + delta1*(1-delta2)*log(dCphi_du1(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi)*(Mgev_pdf(y1,theta1)))
  + delta2*(1-delta1)*log(dCphi_du2(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi)*(Mgev_pdf(y2,theta2)))   
  + (1-delta1)*(1-delta2)*log(Cphi(Mgev_surv(y1,theta1),Mgev_surv(y2,theta2),phi))
  
  
  return(sum(logl))
  
}

data <-load("C:/Users/Reseearcher Mode/Desktop/mydata.RData")
y1 <- tc[,1]
y2 <- tc[,2]
delta1 <- delta[,1]
delta2 <- delta[,2]
n <- length(y1)

#################################################
B <-function(phi,phi1,PARA,y1,y2,delta1,delta2,n)
{
  s <- 0
  N = nrow(PARA)
  for(i in 1:N)
  {
    num = exp(L(PARA[i,1],PARA[i,2],y1,y2,phi,delta1,delta2,n)-L(PARA[i,1],PARA[i,2],y1,y2,phi1,delta1,delta2,n))
    s = s + num
  }
  ratio = s/N
  return(ratio)
}

phi_hat <-optimize(f = B ,phi1=1,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n,lower = 1,upper = 5, maximum = TRUE)
phi_hat <-optimize(f = B ,phi1=0.5,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n,lower = 0.1,upper = 5, maximum = TRUE)
####################################################################################################
load("C:/Users/Reseearcher Mode/Desktop/phi1.RData")
data1 <- samp
load("C:/Users/Reseearcher Mode/Desktop/phi2.RData")
data2 <- samp
load("C:/Users/Reseearcher Mode/Desktop/phi3.RData")
data3 <- samp
load("C:/Users/Reseearcher Mode/Desktop/phi4.RData")
data4 <- samp

phi.skel <- c(0.5,2,4,7)
PARA <- list(data1,data2,data3,data4)
N <- length(data1[,1])
J <- length(phi.skel)
H <- matrix(nrow=N,ncol=J)
for (i in 1:N){
  for (j in 1:J){
    H[i,j] <- L(PARA[[j]][i,1],PARA[[j]][i,2],y1,y2,phi.skel[j],delta1,delta2)
  }
}
H <- exp(H)
bigf <- function(r1,H){ #r1 is a vector of length J-1
  r <- c(1,r1)
  J <- dim(H)[2]
  t1 <- sum(log(H*exp(r)))
  t2 <- sum(log(rowSums(H*exp(r))))*J
  return(t1-t2)
}

init <- c(1,1,1)
r.est <- optim(init,bigf,method="L-BFGS-B",lower=0.001,control=list(fnscale=-1),H=H)$par
r.est <- c(1,r.est)

mat1 <- matrix(nrow=N,ncol=J)
biglist <- list()
for (j in 1:J){
  biglist[[j]] <- mat1
}

for (l in 1:N){
  for (j in 1:J){
    for(i in 1:J){
      biglist[[i]][l,j] <- N*exp(L(PARA[[j]][l,1],PARA[[j]][l,2],y1,y2,phi.skel[i],delta1,delta2))/r.est[i]
    }
  }
}

biglistmat <- matrix(data=0,nrow=N,ncol=J)
for (i in 1:J){
  biglistmat <- biglistmat + biglist[[i]]
}

B.big <-function(phi,PARA,y1,y2,delta1,delta2,biglistmat)
{
  s = 0
  N <- dim(biglistmat)[1]
  J <- dim(biglistmat)[2]
  for (i in 1:N){
    for (j in 1:J){
      temp <- exp(L(PARA[[j]][i,1],PARA[[j]][i,2],y1,y2,phi,delta1,delta2))/biglistmat[i,j]
      s <- s + temp
    }
  }
  return(s)
}

phi.est.final <- optimize(f = B.big ,PARA=PARA,y1=y1,y2=y2,delta1=delta1,delta2=delta2,biglistmat=biglistmat,lower = 0.1,upper = 10, maximum = TRUE)

