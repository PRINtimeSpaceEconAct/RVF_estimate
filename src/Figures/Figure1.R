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

# load("../datasets_it/data_x_wx_allYears.RData")
load("../datasets_it/municipalPopulation_1982_2019.RData")
load("../datasets_it/W_SLL_2001_2011_2019.RData")

#Read the polygon shapefile 
shp_SLL_2001 <- readOGR("../datasets_it/SLL2001_shapefile/SLL2001_wgs84_EURO_rev_20150205.shp")
shp_SLL_2011 <- readOGR("../datasets_it/SLL2011_shapefile/SLL2011_wgs84_EURO.shp")
FUA <- read_excel("../datasets_it/FUA_eurostat.xlsx")

#-------------------------------------------------#
# Pooled observations for all years in subperiods #
#-------------------------------------------------#
df.1984 <- data.pop.all %>% filter(time==1984) %>% mutate(popDensAct = as.numeric(x))
df.2019 <- data.pop.all %>% filter(time==2019) %>% mutate(popDensFut = as.numeric(y))
df <- df.1984 %>% inner_join(df.2019, by="geo")  %>% dplyr::select(c(geo,Comune=Comune.x,popDensAct, popDensFut, pop=pop.x)) %>% na.omit()
#Cities bigger that 100000 inhabitants
df <- df %>% mutate(citiesHighPop=if_else(pop>=20000, 1, 0))

# df <- df %>% mutate(popDiffLog = ifelse((popDensFut-popDensAct)>0, log(popDensFut-popDensAct+1), 
# ifelse((popDensFut-popDensAct)<0,-log(popDensAct-popDensFut+1),0)))
df <- df %>% mutate(popDensActLog=log(popDensAct), popDensFutLog=log(popDensFut), popDiffLog = (log(popDensFut)-log(popDensAct)))

# #FUA as eurostat
# df <- df %>% left_join(FUA, by=c("Comune"="NAME"))
# df <- df %>% mutate(citiesFUA=if_else(!is.na(CODE), 1, 0))
# df <- df %>% dplyr::select(-c(CODE))

shp_sll <- shp_SLL_2011
shp_comm <- shp_comm_2019
shp_comm@data <- shp_comm@data %>% left_join(df, by=c("PRO_COM"="geo"))
shp_comm_sf <- st_as_sf(shp_comm)
cities <- shp_comm_sf %>% filter(cities==1)
cities <- st_intersection(cities, st_union(shp_comm_sf))
cities <- st_cast(cities,"POLYGON")

#Map
qq <- quantile(dfPooled$x, seq(0, 1, by=0.1))

# tm_obj <- tm_shape(shp_comm) + tm_fill("x", title="Population density (log)", breaks=qq, palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5)

tm_obj <- tm_shape(shp_comm) + tm_fill("popDensAct", title="Population density (log)", breaks=qq, palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(cities) +  tm_dots(size = 0.05, col="red")   + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5)

dev.new()
tm_obj
tmap_save(tm = tm_obj,  filename = paste("map_popDensSLL_", year, ".pdf", sep=""))










qq <- quantile(df$popDiffLog, seq(0, 1, by=0.1),na.rm=T)

#qq <- sort(c(qq,0.0))
tm_obj <- tm_shape(shp_comm) + tm_fill("popDiffLog", title="Population density (log) diff", breaks=qq, palette ="RdBu", midpoint = 0) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(cities) + 
    tm_dots(size = 0.1, col="yellow") + tm_shape(shp_sll) + tm_borders("black", alpha=1, lwd=0.5)

dev.new()
tm_obj
# tmap_save(tm = tm_obj,  filename = "mapDiffLog1984_2019.pdf")

