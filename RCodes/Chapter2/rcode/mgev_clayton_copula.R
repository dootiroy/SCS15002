#Generate 1000 samples from clayton copula family

library(copula)

cop2 <- claytonCopula(5, dim = 2)

#cop2 <- gumbelCopula(2.5, dim=2)

#cop2 <- frankCopula(2.5, dim=2)


v <- rCopula(1000, cop2)

pairs(v)
v1 <- 1-v

#Take inverse survival transformation on v to get mgev(0,1,xi) samples.
xi <- 0.5
for(i in 1:length(v1))
{
  #set original value: xi=0.4
  if(xi != 0)
  {
    t =  exp((-1 + (-log(v1))^xi)/xi)    #the survival time
  }
  else
  {
    t = -log(v1)
  }
}
pairs(t)


#Get your censored y and delta
c <- 2 #change this to adjust the censoring rate.
tc <- t
tc[tc > c] <- c
delta <- (t < c)*1
sum(tc[,1]==c)/length(tc[,1]) #censoring rate of y1
sum(tc[,2]==c)/length(tc[,2]) #censoring rate of y2
pairs(tc)

save(tc,delta,file="C:/Desktop stuff/Copula GEV/mydata.RData")



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
#Define pdf of mgev(0,1,xi)
mgev_pdf <- function(t,xi)
{
  
  if(xi == 0)
  {
    f = exp(-t)
  }
  else
  {
    check <- rep(exp(-1/xi),length(t))
    cond <- ((t > check & xi > 0) | (t < check & xi < 0))
    r = (1/t)*((1 + xi*log(t))^(1/xi - 1))*(exp(-(1 + xi*log(t))^(1/xi)))
    f = ifelse(cond,r,0)
  }
  
  return(f)
}

#Define survival function of mgev(0,1,xi)
mgev_surv <- function(t,xi)
{
  if(xi != 0)
  {
    s =  exp(-(1 + xi*log(t))^(1/xi))    #the survival time
  }
  else
  {
    s = exp(-t)
  }
  return(s)
}

L <- function(xi1,xi2,y1,y2,phi,delta1,delta2,n)
{
  
  
  logl <- delta1*delta2*log(dCphi_du1du2(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi)*mgev_pdf(y1,xi1)*mgev_pdf(y2,xi2))
  + delta1*(1-delta2)*log(dCphi_du1(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi)*(mgev_pdf(y1,xi1)))
  + delta2*(1-delta1)*log(dCphi_du2(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi)*(mgev_pdf(y2,xi2)))   
  + (1-delta1)*(1-delta2)*log(Cphi(mgev_surv(y1,xi1),mgev_surv(y2,xi2),phi))
  
  
  return(sum(logl))
  
}

# #Load the simulated data set from Desktop
data <-load("C:/Desktop stuff/Copula GEV/mydata.RData")

y1 <- tc[,1]
y2 <- tc[,2]
delta1 <- delta[,1]
delta2 <- delta[,2]
n <- length(y1)


xi1<-seq(0.1,1,0.01)
Lvalues <- unlist(lapply(xi1, phi=5, L, xi2=0.5, y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n))

plot(xi1, Lvalues)



phi<-seq(0.1,20,0.1)
Lvalues <- unlist(lapply(phi, L, xi1=0.5, xi2=0.5,y1=y1,y2=y2,delta1=delta1,delta2=delta2,n=n))

plot(phi, Lvalues)
