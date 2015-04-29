library(gdata)
rm(list=ls())
setwd("~/Desktop/SCS_Spring_2015/SCS15002")

#Read Old Storm Data#
old <- read.csv("Old.csv")
old1 <- old
#Convert all the important columns imported as factors to numeric#
old1$storm_mag1 <- as.numeric(levels(old$storm_mag))[old$storm_mag]
old1$o_n_g1 <- as.numeric(levels(old$o_n_g))[old$o_n_g]
old1$rain_ph1 <- as.numeric(levels(old$rain_ph))[old$rain_ph]
old1$samp_ph1 <- as.numeric(levels(old$samp_ph))[old$samp_ph]
old1$cod1 <- as.numeric(levels(old$cod))[old$cod]
old1$tss1 <- as.numeric(levels(old$tss))[old$tss]
old1$tp1 <- as.numeric(levels(old$tp))[old$tp]
old1$tkn1 <- as.numeric(levels(old$tkn))[old$tkn]
old1$no31 <- as.numeric(levels(old$no3))[old$no3]
old1$tot_copper1 <- as.numeric(levels(old$tot_copper))[old$tot_copper]
old1$tot_zinc1 <- as.numeric(levels(old$tot_zinc))[old$tot_zinc]
old1$tot_lead1 <- as.numeric(levels(old$tot_lead))[old$tot_lead]

#Read New Storm Data for Industry1 #
ind <- list()
for(i in 1:10){
  dat <- read.csv(paste("ind",i,".csv",sep=''))
  dat <- dat[,1:30]
  #Convert all the important columns imported as factors to numeric#
  ind[[i]] <- dat
  ind[[i]]$o_n_g1 <- as.numeric(levels(dat$o_n_g))[dat$o_n_g]
  ind[[i]]$rain_ph1 <- as.numeric(levels(dat$rain_ph))[dat$rain_ph]
  ind[[i]]$samp_ph1 <- as.numeric(levels(dat$samp_ph))[dat$samp_ph]
  ind[[i]]$cod1 <- as.numeric(levels(dat$cod))[dat$cod]
  ind[[i]]$tss1 <- as.numeric(levels(dat$tss))[dat$tss]
  ind[[i]]$tp1 <- as.numeric(levels(dat$tp))[dat$tp]
  ind[[i]]$tkn1 <- as.numeric(levels(dat$tkn))[dat$tkn]
  ind[[i]]$no31 <- as.numeric(levels(dat$no3))[dat$no3]
  ind[[i]]$tot_copper1 <- as.numeric(levels(dat$tot_copper))[dat$tot_copper]
  ind[[i]]$tot_zinc1 <- as.numeric(levels(dat$tot_zinc))[dat$tot_zinc]
  ind[[i]]$tot_lead1 <- as.numeric(levels(dat$tot_lead))[dat$tot_lead] 
}



