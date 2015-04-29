rm(list=ls())
gpclibPermit()
setwd("~/Desktop/SCS_Spring_2015/SCS15002/")
save(dat_twn_yr, file = "dat_twn_yr.RData")

load("dat_twn_yr.RData")
read.csv("plot.csv")

dat_twn_yr$group_yr <- cut(dat_twn_yr$yr, breaks = c(1995, 2000, 2005, 2010, 2014), 
                    labels = c("1995 - 2000", "2001 - 2005", "2006 - 2010", "2011 - 2014"), include.lowest = TRUE)

dat_hot_twn <- dat_twn_yr %>% group_by(group_yr,NAME10, GEOID_AFF1)%>% summarise(
  ong_grp = mean(ong, na.rm = TRUE),
  cod_grp = mean(cod, na.rm = TRUE),
  tss_grp = mean(tss, na.rm = TRUE),
  tp_grp = mean(tp, na.rm = TRUE),
  tkn_grp = mean(tkn, na.rm = TRUE),
  no3_grp = mean(no3, na.rm = TRUE),
  cop_grp = mean(cop, na.rm = TRUE),
  zinc_grp = mean(zinc, na.rm = TRUE),
  lead_grp = mean(lead, na.rm = TRUE)
)
hot_towns <- filter(dat_hot_twn, sum((tss_grp > 90), (tp_grp > 0.4), (tkn_grp > 2.3),
            (cop_grp > 0.059), (zinc_grp > 0.16), (lead_grp > 0.076), (no3_grp > 1.1)) > 2)

name_hot <- table(hot_towns$group_yr, hot_towns$NAME10)
write.table(name_hot, file = "hot_town_list.csv")
