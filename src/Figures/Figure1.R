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

source("lib/chang")
source("lib/pictureMoranScatterplot.R")
source("lib/changeIstatCodes_adjustOverTime.R")

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

# Figure 1 (a) ----

#1984
pictureMoranScatterplot(dfPooled=dfPooled, df=df.1984, year="1984", shp_comm=shp_comm_2001,  shp_sll=shp_SLL_2001, scatterPlot=TRUE, map=TRUE, typeCities="fua")

# Figure 1 (b) ----

#2019
pictureMoranScatterplot(dfPooled=dfPooled, df=df.2019, year="2019", shp_comm=shp_comm_2019,  shp_sll=shp_SLL_2011, scatterPlot=TRUE, map=TRUE, typeCities="fua")


# Figure 1 (c) ----
load("../datasets_it/data_x_wx_allYears.RData")
load("../datasets_it/W_SLL_2001_2011_2019.RData")
FUA <- read_excel("../datasets_it/FUA_eurostat.xlsx")

df_allYears<- changeIstatCodes_adjustOverTime(df_allYears)
df_allYears <- df_allYears %>% left_join(FUA, by=c("Comune"="NAME"))
df_allYears <- df_allYears %>% mutate(citiesFUA=if_else(!is.na(CODE), 1, 0))
df_allYears <- df_allYears %>% dplyr::select(-c(CODE))
df.1984 <- df_allYears %>% filter(time==1984)
df.2019 <- df_allYears %>% filter(time==2019)
df.1984 <- df.1984 %>% dplyr::rename(popDensAct=popDensity)
df.2019 <- df.2019 %>% dplyr::rename(popDensFut=popDensity)

df <- df.1984 %>% inner_join(df.2019, by="geo")  %>% dplyr::select(c(geo, "Comune"="Comune.x",popDensAct, popDensFut,pop.x, citiesFUA=citiesFUA.x))
df <- df %>% mutate(popDiffLog = (log(popDensFut)-log(popDensAct)))

shp_sll <- shp_SLL_2011
shp_comm <- shp_comm_2019
shp_comm@data <- shp_comm@data %>% left_join(df, by=c("PRO_COM_T"="geo"))
shp_comm_sf <- st_as_sf(shp_comm)
cities <- shp_comm_sf %>% filter(citiesFUA==1)
cities <- st_intersection(cities, st_union(shp_comm_sf))
cities <- st_cast(cities,"POLYGON")

qq <- quantile(df$popDiffLog, seq(0, 1, by=0.1))
tm_obj <- tm_shape(shp_comm) + tm_fill("popDiffLog", title="Population density (log) diff", breaks=qq, palette ="RdBu", midpoint = 0) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(cities) + 
    tm_dots(size = 0.1, col="red") + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5)

dev.new()
tm_obj
# tmap_save(tm = tm_obj,  filename = "mapDiffLog1984_2019.pdf")
