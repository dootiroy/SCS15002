#Define Clayton Copula structure
# Note: Zero density regions put as a small +ve number for numerical convenience
Cphi <- function(u1,u2,phi)
{
  b <- exp((-1/phi)*log((u1^(-phi))+(u2^(-phi))-1))
  ifelse(b>1e-100, b, 1e-100)
}

#Define partial derivative structure
dCphi_du1 <- function(u1,u2,phi)
{
  b <- exp((-1/phi)*log((u1^(-phi))+(u2^(-phi))-1))
  a <- exp((-phi-1)*log(u1) - (1/phi + 1)*log(u1^(-phi)+u2^(-phi)-1))
  ifelse( (b>0 & a>1e-100), a, 1e-100)
}
dCphi_du2 <- function(u1,u2,phi)
{
  b <- exp((-1/phi)*log((u1^(-phi))+(u2^(-phi))-1))
  a <- exp((-phi-1)*log(u2) - (1/phi + 1)*log(u1^(-phi)+u2^(-phi)-1))
  ifelse( (b>0 & a>1e-100), a, 1e-100)
}
dCphi_du1du2 <- function(u1,u2,phi)
{
  b <- exp((-1/phi)*log((u1^(-phi))+(u2^(-phi))-1))
  a <- exp(log(1+phi) - (phi+1)*log(u1*u2) -(1/phi+2)*log(u1^(-phi)+u2^(-phi)-1))
  ifelse( (b>0 & a>1e-100), a, 1e-100)
}
#Define pdf of mgev(mu,sigma,xi)
mgev_pdf <- function(t,mu,sigma,xi)
{
  
  if(xi == 0)
  {
    s = exp(-exp((log(t)-mu)/sigma))
    f = s*exp((log(t) - mu)/sigma)*(1/(sigma*t))
  }
  else
  {
    # check <- rep(exp(mu -sigma/xi),length(t))
    cond <- (t > exp(mu - sigma/xi) & xi > 0) | (t < exp(mu - sigma/xi) & xi < 0)
    S <- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
    r <- (S)*(1/(sigma*t))*((1+xi*(log(t)-mu)/sigma)^(1/xi-1)) #density function
    f = ifelse((cond & r>1e-100),r,1e-100)
  }
  
  return(f)
}

#Define survival function of mgev(0,1,xi)
mgev_surv <- function(t,mu,sigma,xi) {
  if(xi != 0){
    cond <- (t > exp(mu - sigma/xi) & xi > 0) | (t < exp(mu - sigma/xi) & xi < 0)
    s1 <- exp(-(1+xi*(log(t)-mu)/sigma)^(1/xi)) #the survival function
    s <- ifelse(cond,s1,1)
  }
  else{
    s = exp(-exp((log(t)-mu)/sigma))
  }
  return(s)
}

L <- function(mu1,sigma1, xi1, mu2,sigma2,xi2,y1,y2,phi,delta1,delta2,n){
  
  
  logl <- delta1*delta2*log(dCphi_du1du2(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi)*mgev_pdf(y1,mu1,sigma1,xi1)*mgev_pdf(y2,mu2,sigma2,xi2))
  + delta1*(1-delta2)*log(dCphi_du1(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi)*(mgev_pdf(y1,mu1,sigma1,xi1)))
  + delta2*(1-delta1)*log(dCphi_du2(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi)*(mgev_pdf(y2,mu2,sigma2,xi2)))   
  + (1-delta1)*(1-delta2)*log(Cphi(mgev_surv(y1,mu1,sigma1,xi1),mgev_surv(y2,mu2,sigma2,xi2),phi))
  
  
  return(sum(logl))
  
}

L_wrap2 <- function(theta, x, y1,y2,phi,delta1,delta2,n) {
  # theta is a vector here
  sigma1 <- theta[1]
  xi1 <- theta[2]
  sigma2 <- theta[3]
  xi2 <- theta[4]
  b01 <- theta[5]
  b11 <- theta[6]
  b02 <- theta[7]
  b12 <- theta[8]
  # mu's are now vectors, changing with each subject
  mu1 <- b01 + b11*x
  mu2 <- b02 + b12*x
  L(mu1,sigma1, xi1, mu2,sigma2,xi2,y1,y2,phi,delta1,delta2,n)
}

L_wrap3 <- function(theta, x, y1,y2,phi,delta1,delta2,n) {
  # theta is a list here
  # mu's are now vectors, changing with each subject
  mu1 <- theta$b01 + theta$b11*x
  mu2 <- theta$b02 + theta$b12*x
  L(mu1, theta$sigma1, theta$xi1, mu2, theta$sigma2, theta$xi2,
    y1, y2, phi, delta1, delta2, n)
}
