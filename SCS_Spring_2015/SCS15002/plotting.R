library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
rm(list=ls())
setwd("~/Desktop/SCS_Spring_2015/SCS15002")
# read in the data
load("combined.RData")
source("multiplot_code.R")

# Goal1: Find time trends for each industries
# Loess curves
vec <- c("tss1", "tp1", "tkn1", "no31", "tot_copper1", "tot_zinc1", "tot_lead1")
inter <- c(90, 0.4, 2.3, 1.1, 0.059, 0.16, 0.076)
tplot <- list()
fin_plot <- list()
fin_plot1 <- list()


# compare sectors
splot <- ggplot(dat, aes(x = samp_date1, y = tss1)) +
         stat_smooth(se=FALSE, method = "auto", span=0.9) + 
         geom_hline(yintercept=90)  +
         theme_bw()
splot

# compare sectors
splot <- ggplot(dat, aes(x = samp_date1, y = tss1, col = Sector1)) +
  ggtitle("Time Trend for Total Suspended Solids by Sectors") +
  xlab("Year") +
  ylab("TSS mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=90)  +
  theme_bw()
splot

splot <- ggplot(dat, aes(x = samp_date1, y = tp1, col = Sector1)) +
  ggtitle("Time Trend for Total Phosphorus by Sectors") +
  xlab("Year") +
  ylab("TP mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=0.4)  +
  theme_bw()
splot

splot <- ggplot(dat, aes(x = samp_date1, y = tkn1, col = Sector1)) +
  ggtitle("Time Trend for Total Potassium Nitrate by Sectors") +
  xlab("Year") +
  ylab("TKN mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=2.3)  +
  theme_bw()
splot

splot <- ggplot(dat, aes(x = samp_date1, y = no31, col = Sector1)) +
  ggtitle("Time Trend for Total Nitrate by Sectors") +
  xlab("Year") +
  ylab("NO3 mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=1.1)  +
  theme_bw()
splot

splot <- ggplot(dat, aes(x = samp_date1, y = tot_copper1, col = Sector1)) +
  ggtitle("Time Trend for Total Copper by Sectors") +
  xlab("Year") +
  ylab("Cu mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=0.059)  +
  theme_bw()
splot

splot <- ggplot(dat, aes(x = samp_date1, y = tot_zinc1, col = Sector1)) +
  ggtitle("Time Trend for Total Zinc by Sectors") +
  xlab("Year") +
  ylab("Zn mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=0.16)  +
  theme_bw()
splot

splot <- ggplot(dat, aes(x = samp_date1, y = tot_lead1, col = Sector1)) +
  ggtitle("Time Trend for Total Lead by Sectors") +
  xlab("Year") +
  ylab("Pb mg/L") +
  stat_smooth(se=FALSE, method = "loess", span=0.9) + 
  geom_hline(yintercept=0.076)  +
  theme_bw()
splot
################################################################
tss <- ggplot(dat, aes(x = samp_date1, y = tss1)) +
       ggtitle("Time Trend for Total Suspended Solids") +
       xlab("Year") +
       ylab("TSS mg/L") +
       stat_smooth(method = "loess", span=0.9) + 
       geom_hline(yintercept=90)

tp <- ggplot(dat, aes(x = samp_date1, y = tp1)) +
      ggtitle("Time Trend for Total Phosphorus") +
      xlab("Year") +
      ylab("TP mg/L") +
      stat_smooth(method = "loess", span=0.9) + 
      geom_hline(yintercept=0.4)

tkn <- ggplot(dat, aes(x = samp_date1, y = tkn1)) +
       ggtitle("Time Trend for Total Potassium Nitrate") +
       xlab("Year") +
       ylab("TKN mg/L") +
       stat_smooth(method = "loess", span=0.9) + 
       geom_hline(yintercept=2.3)

no3 <- ggplot(dat, aes(x = samp_date1, y = no31)) +
       ggtitle("Time Trend for Total Nitrate") +
       xlab("Year") +
       ylab("NO3 mg/L") +
       stat_smooth(method = "loess", span=0.9) + 
       geom_hline(yintercept= 1.1)

plot1 <- grid.arrange(tss, tp, tkn, no3, ncol = 1)



