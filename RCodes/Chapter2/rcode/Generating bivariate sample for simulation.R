#Generate 10,000 samples from clayton copula family
library(copula)
cop2 <- claytonCopula(2.5, dim = 2)
v <- rCopula(100, cop2)
pairs(v)
v1 <- 1-v
#Take inverse survival transformation on v to get GEV samples.
xi <- 0.55
for(i in 1:length(v1))
  {
    #set original value: u=0,sigma=1,xi=0.5 for MGEV distribution
    logt<- (-1+(-log(v1))^(-xi))/xi 
    t<-exp(logt)  #the survival time
    
    }
pairs(t)
#Get your censored y and delta
c <- 9 #change this to adjust the censoring rate.
tc <- t
tc[tc > c] <- c
delta <- (t < c)*1
sum(tc[,1]==c)/length(tc[,1]) #censoring rate of y1
sum(tc[,2]==c)/length(tc[,2]) #censoring rate of y2
pairs(tc)

save(tc,delta,file="C:/Users/Dooti/Desktop/Copula GEV/mydata.RData")
