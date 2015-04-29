library(gdata)
library(data.table)
rm(list=ls())
setwd("~/Desktop/SCS_Spring_2015/SCS15002")

#Name of all the variables in Old.csv

# NAME  rwater	TOWN	BASIN	GSI	SIC	Municipalities	samp_location	samp_date	
# storm_mag	storm_dur	prev_storm	rain_ph	days_last_storm	o_n_g	samp_ph	
# cod	tss	tp	tkn	no3	fec_coli	tot_copper	tot_zinc	tot_lead	24_hr_lc50	48_hr_lc50	hardness

# cc_old <- c("factor", "factor", "factor", "factor", "factor", "factor", "factor", "factor", "factor",
#             "numeric", "numeric", "factor",  "numeric", "numeric", "numeric", "numeric", "numeric",
#             "numeric", "numeric", "numeric", "numeric", "factor", "numeric", "numeric", "numeric", "factor", 
#             "factor", "factor")

cc_old <- c("factor", "factor", "factor", "factor", "factor", "factor", "factor", "factor", "factor",
            "numeric", "numeric", "factor",  "numeric", "numeric", "numeric", "numeric", "factor",
            "factor", "numeric", "numeric", "numeric", "factor", "numeric", "numeric", "numeric", "factor", 
            "factor", "factor")


old <- read.csv("Old.csv",colClasses = cc_old)
old$cod <- as.numeric(paste(old$cod))
old$tss <- as.numeric(paste(old$tss))
old$samp_date1 <- as.Date(old$samp_date, format = "%m/%d/%y")

ind <- list()
for(i in 1:12){
  dat <- read.csv(paste("ind",i,".csv",sep=''), colClasses = "factor")
  dat <- dat[,1:30]
  
  #Convert all the important columns imported as factors to numeric#
  ind[[i]] <- dat
  ind[[i]]$SIC1 <- dat$SIC
  ind[[i]]$samp_date1 <- as.Date(dat$samp_date, format = "%m/%d/%y")
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
  ind[[i]]$Sector1 <- factor(i)
}

# Combine new data into one frame. Relevant cols for goal 1
dat_c <- rbindlist(ind)
dat_c <- data.frame(dat_c)
dat_c <- dat_c[,c(4,31:44)]

# old1 is compatible with new
old1 <- data.frame(old$TOWN)
names(old1) <- "town"
old1$SIC1 <- old$SIC
old1$samp_date1 <- old$samp_date1
old1$o_n_g1 <- old$o_n_g
old1$rain_ph1 <- old$rain_ph
old1$samp_ph1 <- old$samp_ph
old1$cod1 <- old$cod
old1$tss1 <- old$tss
old1$tp1 <- old$tp
old1$tkn1 <- old$tkn
old1$no31 <- old$no3
old1$tot_copper1 <- old$tot_copper
old1$tot_zinc1 <- old$tot_zinc
old1$tot_lead1 <- old$tot_lead


# Get Sectors for Old using SIC from Old and SIC1 and Sector1 in new data
SICs <- dat_c$SIC1
Secs <- dat_c$Sector1
SIC_bad <- character()
old1$Sector1 <- rep(NA, nrow(old1))

for(i in 1:nrow(old1)) {
  if(!is.na(which( SICs %in% old1$SIC[i] )[1] )) {
    old1$Sector1[i] <- Secs[which( SICs %in% old1$SIC[i] )[1]]
  } else {
    old1$Sector1[i] <- NA
    SIC_bad <- append(SIC_bad, as.character(old1$SIC[i]))
  }
}
write.csv(unique(SIC_bad), file = "badSIC.csv")

old1$Sector1 <- factor(old1$Sector1)
#names(old1)[1] <- "SIC1"

# Concatenate the two dataframes
dat <- rbind(old1, dat_c)
save(dat, file = "combined.RData")

#Relationship between discharge_quality with location and sector

tsslm <- lm(dat$tss1 ~ dat$town + dat$Sector1 + dat$town*dat$Sector1)
tss.res = resid(tsslm)
fitted <- tsslm$fitted.values
qqnorm(tss.res)
plot(tsslm)
plot(dat$tss1)
