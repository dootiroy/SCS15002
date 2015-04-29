library(ggplot2)
library(dplyr)
library(Cairo)
library(rgdal)
library(ggmap)
library(scales)
library(maptools)
library(animation)
rm(list=ls())
gpclibPermit()
setwd("~/Desktop/SCS_Spring_2015/SCS15002/Pollutants/")

# twn contains shapefile info for CT towns
twn <- readOGR(dsn = "townct_37800_0000_2010_s100_census_1_shp/wgs84/", 
               layer = "townct_37800_0000_2010_s100_census_1_shp_wgs84")
twn <- fortify(twn, region = "GEOID_AFF1")
# twn2 <- readOGR(dsn = "townct_37800_0000_2010_s100_census_1_shp/nad83/", 
#                 layer = "townct_37800_0000_2010_s100_census_1_shp_nad83_meters")
#twn2 <- fortify(twn2)

twngeo <- read.csv("towns_geo.csv", stringsAsFactors = FALSE)
twngeo <- select(twngeo, GEOID_AFF1, NAME10, INTPTLAT10, INTPTLON10)
twngeo <- filter(twngeo, NAME10 != "County subdivisions not defined")
twnid <- read.csv("towns.csv", colClasses = c("integer","character"))
twn_fin <- read.csv("plot.csv")
twn_geoid <- inner_join(twngeo, twnid, by = c("NAME10" = "town_name"))
twn_plot <- inner_join(twn_geoid, twn_fin, by = c("NAME10" = "Town"))

load("combined.RData")
dat$townid <- as.integer(paste(dat$town))

# Add GEOID10 to dat
dat <- left_join(dat, twn_geoid, by = "townid")

# For mapping, need to aggregate by towns and year
# Currently using simple mean for aggregation
# Remove rows with no year information
dat$yr <- as.numeric(format(dat$samp_date1, "%Y"))
dat1 <- dat %>% select(o_n_g1:tot_lead1, GEOID_AFF1, NAME10, yr)
dat_twn_yr <- dat1 %>% group_by(yr, NAME10, GEOID_AFF1) %>%
              summarise(
                ong = mean(o_n_g1, na.rm = TRUE),
                cod = mean(cod1, na.rm = TRUE),
                tss = mean(tss1, na.rm = TRUE),
                tp = mean(tp1, na.rm = TRUE),
                tkn = mean(tkn1, na.rm = TRUE),
                no3 = mean(no31, na.rm = TRUE),
                cop = mean(tot_copper1, na.rm = TRUE),
                zinc = mean(tot_zinc1, na.rm = TRUE),
                lead = mean(tot_lead1, na.rm = TRUE)
              ) %>%
              filter(!is.na(yr))
dat_twn_yr$ong_dis <- cut(dat_twn_yr$ong, breaks = c(0, 5, 25, 420.4), labels = c("< 5", "5 - 25", "> 25"), include.lowest = TRUE)
dat_twn_yr$cod_dis <- cut(dat_twn_yr$cod, breaks = c(0, 75, 100, 2873), labels = c("< 75", "75 - 100", "> 100"), include.lowest = TRUE)
dat_twn_yr$tss_dis <- cut(dat_twn_yr$tss, breaks = c(0, 90, 200, 6417), labels = c("< 90", "90 - 200", "> 200"), include.lowest = TRUE)
dat_twn_yr$tp_dis <- cut(dat_twn_yr$tp, breaks = c(0, 0.4, 3, 26.03), labels = c("< 0.4", "0.4 - 3", "> 3"), include.lowest = TRUE)
dat_twn_yr$tkn_dis <- cut(dat_twn_yr$tkn, breaks = c(0, 2.3, 10, 105.6), labels = c("< 2.3", "2.3 - 10", "> 10"), include.lowest = TRUE)
dat_twn_yr$no3_dis <- cut(dat_twn_yr$no3, breaks = c(0, 1.1, 10, 170), labels = c("< 1.1", "1.1 - 10", "> 10"), include.lowest = TRUE)
dat_twn_yr$cop_dis <- cut(dat_twn_yr$cop, breaks = c(0, 0.059, 2, 8.729), labels = c("< 0.059", "0.059 - 2", "> 2"), include.lowest = TRUE)
dat_twn_yr$zinc_dis <- cut(dat_twn_yr$zinc, breaks = c(0, 0.16, 3, 38), labels = c("< 0.16", "0.16 - 3", "> 3"), include.lowest = TRUE)
dat_twn_yr$lead_dis <- cut(dat_twn_yr$lead, breaks = c(0, 0.076, 1.5, 3.319), labels = c("< 0.076", "0.076 - 1.5", "> 1.5"), include.lowest = TRUE)

# Testing for year = 2000
# for ong
dat_2010 <- dat_twn_yr %>% filter(yr == 2008)
names(dat_2010)[3] <- "id"
mapdat_2010 <- left_join(twn, dat_2010)
# 06000US0900100000, 06000US0900740710
mapdat_2010 <- filter(mapdat_2010, id!= "06000US0900100000")
# put NA ong as zero, change for other pollutants
mapdat_2010$ong_dis[which(is.na(mapdat_2010$ong_dis))] <- "< 75"
# mapdat_2010 <- filter(mapdat_2010, !is.na(NAME10))
mycols = c("#FFCCCC", "#FF3333", "#CC0000")
p <- ggplot() + geom_polygon(data = mapdat_2010, aes(x = long, y = lat, 
    group = group, fill = ong_dis), color = "white", size = 0.5) +
    coord_map() + 
    scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
    xlab("") + ylab("") +
    ggtitle("Chemical Oxygen Demand Pollutant Loadings in CT Towns in 2010 ") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_blank(), 
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), legend.position = "top",
        legend.title = element_text(size=14) , legend.text = element_text(size=12),
        plot.title=element_text(family="Times",face="bold", size=20)) +
    with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
p


# Animation
# ONG
# replace NAs by zero
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$ong_dis[which(is.na(mapdat_yr$ong_dis))] <- "< 5"

    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = ong_dis), 
                                 color = "white", size = 0.5) +
        coord_map() + 
       #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
        scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
        xlab("") + ylab("") + 
        ggtitle(paste("Oil and Grease Pollutant Loadings in CT Towns in ", years[i], sep="")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
        with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plotong", imgdir = "autoplots", htmlfile = "yearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "ONG Pollutant Loadings in CT")
dev.off()

# COD
# replace NAs by zero
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$cod_dis[which(is.na(mapdat_yr$cod_dis))] <- "< 75"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = cod_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Chemical Oxygen Demand Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plotcod", imgdir = "autoplots", htmlfile = "codyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "COD Pollutant Loadings in CT")
dev.off()

#TSS
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$tss_dis[which(is.na(mapdat_yr$tss_dis))] <- "< 90"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = tss_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
        ggtitle(paste("Total Suspended Solids Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plottss", imgdir = "autoplots", htmlfile = "tssyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "Tss Pollutant Loadings in CT")
dev.off()

#TP
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$tp_dis[which(is.na(mapdat_yr$tp_dis))] <- "< 0.4"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = tp_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Total Phosphates Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plottp", imgdir = "autoplots", htmlfile = "tpyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "Tp Pollutant Loadings in CT")
dev.off()

#TKN
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$tkn_dis[which(is.na(mapdat_yr$tkn_dis))] <- "< 2.3"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = tkn_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Total Potassium Nitrate Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plottkn", imgdir = "autoplots", htmlfile = "tknyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "Tkn Pollutant Loadings in CT")
dev.off()

#No3
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$no3_dis[which(is.na(mapdat_yr$no3_dis))] <- "< 1.1"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = no3_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Total Nitrate Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plotno3", imgdir = "autoplots", htmlfile = "no3yearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "No3 Pollutant Loadings in CT")
dev.off()

#Cu
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$cop_dis[which(is.na(mapdat_yr$cop_dis))] <- "< 0.059"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = cop_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Total Copper Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plotcop", imgdir = "autoplots", htmlfile = "copyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "Cu Pollutant Loadings in CT")
dev.off()

#Zinc
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$zinc_dis[which(is.na(mapdat_yr$zinc_dis))] <- "< 0.16"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = zinc_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Total Zinc Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plotzn", imgdir = "autoplots", htmlfile = "znyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "Zn Pollutant Loadings in CT")
dev.off()

#Lead
years <- unique(dat_twn_yr$yr)
saveHTML({
  for (i in seq_along(years)) {
    # create a choropleth for each year
    dat_yr <- dat_twn_yr %>% filter(yr == years[i])
    names(dat_yr)[3] <- "id"
    mapdat_yr <- left_join(twn, dat_yr)
    mapdat_yr <- filter(mapdat_yr, id!= "06000US0900100000")
    # put NA ong as zero, change for other pollutants
    mapdat_yr$lead_dis[which(is.na(mapdat_yr$lead_dis))] <- "< 0.076"
    
    mycols = c("#FFCCCC", "#FF3333", "#990000")
    
    p <- ggplot() + geom_polygon(data = mapdat_yr, aes(x = long, y = lat, group = group, fill = lead_dis), 
                                 color = "white", size = 0.5) +
      coord_map() + 
      #scale_fill_gradient2(high = "red", "Avg. Oil and Grease mg/l", limits = c(0,425)) +
      scale_fill_manual(values = mycols, name = "Pollutant Level in mg/l") +
      xlab("") + ylab("") + 
      ggtitle(paste("Total Lead Pollutant Loadings in CT Towns in ", years[i], sep="")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), legend.position = "top",
            legend.title = element_text(size=14) , legend.text = element_text(size=12),
            plot.title=element_text(family="Times",face="bold", size=20)) +
      with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
    print(p)
  }
}, img.name = "plotpb", imgdir = "autoplots", htmlfile = "pbyearplots.html", 
outdir = getwd(), autobrowse = FALSE, ani.height = 800, ani.width = 1200, 
verbose = FALSE, autoplay = TRUE, title = "Zn Pollutant Loadings in CT")
dev.off()



# Testing

# 06000US0900100000, 06000US0900740710

names(twn_plot)[1] <- "id"
twn_plot1 <- left_join(twn, twn_plot)
twn_plot2 <- filter(twn_plot1, id!= "06000US0900100000")

# put NA ong as zero, change for other pollutants
twn_plot2$temp[which(is.na(twn_plot2$temp))] <- "5"
#mapdat_2010 <- filter(mapdat_2010, !is.na(NAME10))
mycols = c("red", "orange", "gold4","yellow", "grey92")
p <- ggplot() + geom_polygon(data = twn_plot2, aes(x = long, y = lat, 
                                                     group = group, fill = temp), color = "white", size = 0.5) +
  coord_map() + 
  scale_fill_manual(values = mycols, name = "Benchmark Exceeded",
                    labels=c("For Past 20 Years", "For Past 15 Years", "For Past 10 Years", "For Past 5 Years", "Not Consistently")) +
  xlab("") + ylab("") +
  ggtitle("Connecticut Towns Exceeding Benchmark for 3/7 Safety Parameters between 1995 - 2014") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_blank(), 
        axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), legend.position = "top",
        legend.title = element_text(size=14) , legend.text = element_text(size=12),
        plot.title=element_text(family="Times",face="bold", size=20)) +
  with(twngeo, annotate(geom="text", x = INTPTLON10, y=INTPTLAT10, label = NAME10, size = 2.5))
p