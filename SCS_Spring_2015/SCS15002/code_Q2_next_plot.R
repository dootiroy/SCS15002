for(i in 1:1){
  
  dat_sub <- filter(dat, Sector1 == i)
  
  tss <- ggplot(dat_sub, aes(x = samp_date1, y = tss1)) + 
    stat_smooth(SE = TRUE, method = "loess", span =0.9) + 
    ggtitle("Time Trend for Total Suspended Solids") +
    xlab("Year") +
    ylab("TSS mg/L") +
    #geom_point(alpha=1/5) +
    geom_hline(yintercept= 90) 
    tss <- tss + coord_cartesian(ylim=c(0, 200))
  
  tp <- ggplot(dat_sub, aes(x = samp_date1, y = tp1)) +
    ggtitle("Time Trend for Total Phosphorus") +
    xlab("Year") +
    ylab("TP mg/L") +
    stat_smooth(SE = TRUE, method = "loess", span=0.9) + 
    #geom_point(alpha=1/5) +
    geom_hline(yintercept=0.4)
    tp <- tp + coord_cartesian(ylim=c(0, 0.75))
  
  tkn <- ggplot(dat_sub, aes(x = samp_date1, y = tkn1)) +
    ggtitle("Time Trend for Total Potassium Nitrate") +
    xlab("Year") +
    ylab("TKN mg/L") +
    stat_smooth(SE = TRUE, method = "loess", span=0.9) + 
    #geom_point(alpha=1/5) +
    geom_hline(yintercept=2.3)
    tkn <- tkn + coord_cartesian(ylim=c(1.5, 3.5))
  
  no3 <- ggplot(dat_sub, aes(x = samp_date1, y = no31)) +
    ggtitle("Time Trend for Total Nitrate") +
    xlab("Year") +
    ylab("NO3 mg/L") +
    stat_smooth(SE = TRUE, method = "loess", span=0.9) + 
    #geom_point(alpha=1/5) +
    geom_hline(yintercept= 1.1)
    no3 <- no3 + coord_cartesian(ylim=c(0.5, 1.5))
  
  cop <- ggplot(dat_sub, aes(x = samp_date1, y = tot_copper1)) +
    ggtitle("Time Trend for Total Copper Content") +
    xlab("Year") +
    ylab("Cu mg/L") +
    stat_smooth(SE = TRUE, method = "loess", span=0.9) + 
    #geom_point(alpha=1/5) +
    geom_hline(yintercept= 0.059)
    cop <- cop + coord_cartesian(ylim=c(0, 0.25))
  
  zin <- ggplot(dat_sub, aes(x = samp_date1, y = tot_zinc1)) +
    ggtitle("Time Trend for Total Zinc Content") +
    xlab("Year") +
    ylab("Zn mg/L") +
    stat_smooth(SE = TRUE, method = "loess", span=0.9) + 
    #geom_point(alpha=1/5) +
    geom_hline(yintercept= 0.16)
    zin <- zin + coord_cartesian(ylim=c(0, 0.75))
  
  lead <- ggplot(dat_sub, aes(x = samp_date1, y = tot_lead1)) +
    ggtitle("Time Trend for Total Lead Content") +
    xlab("Year") +
    ylab("Pb mg/L") +
    stat_smooth(SE = TRUE, method = "loess", span=0.9) + 
    #geom_point(alpha=1/5) +
    geom_hline(yintercept= 0.076)
    lead <- lead + coord_cartesian(ylim=c(0, 0.1))
  
  fin_plot1[[i]] <- grid.arrange(arrangeGrob(tss,tp,tkn,no3,cop,zin,lead, ncol = 2, nrow = 4),
                                main=textGrob(paste("Sector:",i,sep =""), gp=gpar(fontsize=11,font=3)))
}