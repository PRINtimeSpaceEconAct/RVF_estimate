rm(list = ls())

library(dplyr)
library(sm)
library(maptools)
library(tigris)
library(tmap)
library(ggplot2)
library(plm)
library(spectralGP)
library(readxl)
library(foreach)
library(RColorBrewer)
library(spdep)
library(wesanderson)
library(rgeos)
library(plyr)
library(classInt)
library(raster)
library(rgdal)
library(xtable)
library(funrar)

library(doMC)
registerDoMC(cores=6)
library(expm)
library(data.table)

source("lib/pictureMoranScatterplot.R")
source("lib/dfForPicture.R")

load("../datasets_it/data_x_wx_allYears.RData")
load("../datasets_it/W_SLL_2001_2011_2019.RData")

#Read the polygon shapefile 
shp_SLL_2001 <- readOGR("../datasets_it/SLL2001_shapefile/SLL2001_wgs84_EURO_rev_20150205.shp")
shp_SLL_2011 <- readOGR("../datasets_it/SLL2011_shapefile/SLL2011_wgs84_EURO.shp")
FUA <- read_excel("../datasets_it/FUA_eurostat.xlsx")

#-------------------------------------------------#
# Pooled observations for all years in subperiods #
#-------------------------------------------------#
dfPooled <- dfForPicture(df=df_allYears, fua=FUA, years=c(1984:1991, 1994:2001, 2002, 2004:2011, 2012:2019))
df.1984 <- dfPooled %>% filter(time==1984)
df.2019 <- dfPooled %>% filter(time==2019)

# Figure 3 (a) ----

#1984
pictureMoranScatterplot(dfPooled=dfPooled, df=df.1984, year="1984", shp_comm=shp_comm_2001,  shp_sll=shp_SLL_2001, scatterPlot=TRUE, map=TRUE, typeCities="highPop")
# dev.copy2eps (file=paste("moranScatterplot_",1984,".eps",sep=""))

# Figure 3 (b) ----

#2019
pictureMoranScatterplot(dfPooled=dfPooled, df=df.2019, year="2019", shp_comm=shp_comm_2019,  shp_sll=shp_SLL_2011, scatterPlot=TRUE, map=TRUE, typeCities="highPop")
# dev.copy2eps (file=paste("moranScatterplot_",2019,".eps",sep=""))
