# Forecast using estimated RVF 

rm(list = ls())

library(dplyr)
library(wordspace)
library(sm)
library(pracma)
library(doMC)
registerDoMC(cores=10)


source("lib/forward_estimation_RSS.R")
source("lib/forward_estimation_RSS_bootstrap.R")

load("../datasets_it/data_x_wx_allYears.RData")



#-----------------#
#   1984-2019     #
#-----------------#
f_name = "1984-2019"
load(paste("../datasets_it/randomVectorField_estimation_viaRSS_alphaOpt_TRUE",f_name,".RData",sep=""))
load(paste("../datasets_it/bootstrapAnalysis",f_name,".RData",sep=""))
dfTuttoForward = df_allYears %>% filter(time == 2019)



### Forward mean arrows 
# number of discrete periods. Each periods is of length 2019-1984 = 35 years
nP=200

xMat = matrix(0, nrow = nrow(dfTuttoForward), ncol=nP)
WxMat = matrix(0, nrow = nrow(dfTuttoForward), ncol=nP)
xMat[,1] = dfTuttoForward$logPopDensity
WxMat[,1] = dfTuttoForward$WlogPopDensity


for (i in 2:nP){
    print(i)
    randomVectorField_estimation_viaRSS$estimationRVF$y_Rel = matrix(rep(xMat[,i-1],2),ncol=2)
    randomVectorField_estimation_viaRSS$estimationRVF$Wy_Rel = matrix(rep(WxMat[,i-1],2),ncol=2)
    
    stepForward = forward_estimation_RSS(randomVectorField_estimation_viaRSS$estimationRVF,continuousForward=TRUE,distributionFutDensity=NA)
    
    xMat[,i] = stepForward$obs_pred[,1]
    WxMat[,i] = stepForward$obs_pred[,2]
}
save(xMat,WxMat,file="xMatWxMat.RData")

# precomputed 
# load(file="../datasets_it/xMatWxMat.RData")


### Forward for each bootstrap sample

nP=200
nBoot = ncol(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot)

xMatBoot = matrix(NA, nrow=length(dfTuttoForward$logPopDensity), ncol = nBoot)
WxMatBoot = matrix(NA, nrow=length(dfTuttoForward$WlogPopDensity), ncol = nBoot)

combinedStepForwardBoot <- rbind(xMatBoot, WxMatBoot)

print("Bootstrap")
combinedStepForwardBoot <- foreach(b=1:nBoot, .combine=cbind, .verbose=TRUE) %dopar% {
  
  print(paste("Bootstrap number equal to ",b))
    
    gc()
  
    xMat = matrix(0, nrow = nrow(dfTuttoForward), ncol=nP)
    WxMat = matrix(0, nrow = nrow(dfTuttoForward), ncol=nP)
    
    xMat[,1] = dfTuttoForward$logPopDensity
    WxMat[,1] = dfTuttoForward$WlogPopDensity
    
    for (i in 2:nP){
        print(i)
        randomVectorField_estimation_viaRSS$estimationRVF$y_Rel = matrix(rep(xMat[,i-1],2),ncol=2)
        randomVectorField_estimation_viaRSS$estimationRVF$Wy_Rel = matrix(rep(WxMat[,i-1],2),ncol=2)
        
        stepForward = forward_estimation_RSS_bootstrap(randomVectorField_estimation_viaRSS$estimationRVF,bootstrapAnalysis, continuousForward=TRUE,distributionFutDensity=NA,b)
        
        xMat[,i] = stepForward$obs_pred[,1]
        WxMat[,i] = stepForward$obs_pred[,2]
    }
    
    c(as.vector(xMat[,nP]), as.vector(WxMat[,nP]))
}    

xMatBoot <- combinedStepForwardBoot[1:length(dfTuttoForward$logPopDensity),]
WxMatBoot <- combinedStepForwardBoot[(length(dfTuttoForward$logPopDensity)+1):nrow(combinedStepForwardBoot),]

save(xMatBoot,WxMatBoot,file="xMatWxMatBoot.RData")

# precomputed
# load("../datasets_it/xMatWxMatBoot.RData")


randomVectorField_estimation_viaRSS$estimationRVF$y_Rel = cbind(xMat[,1],xMat[,nP])
randomVectorField_estimation_viaRSS$estimationRVF$Wy_Rel = cbind(WxMat[,1],WxMat[,nP])
df = data.frame(xAct = xMat[,1],wxAct = WxMat[,1], xFut = xMat[,nP], wxFut = WxMat[,nP])


geoRVF = rownames(randomVectorField_estimation_viaRSS$xObs)
load("../datasets_it/data_x_wx_allYears.RData")
df.1 <- df_allYears[df_allYears$time == 1984,]
df.2 <- df_allYears[df_allYears$time == 2019,]
df.1 <- df.1 %>% dplyr::rename(xAct=logPopDensity, wxAct=WlogPopDensity)
df.2 <- df.2 %>% dplyr::rename(xFut=logPopDensity, wxFut=WlogPopDensity)
df <- df.2 %>% dplyr::left_join(df.1, by="geo") %>% dplyr::select(c(geo, xAct, wxAct, xFut, wxFut,cities.x,Comune.x)) 
df <- df %>% dplyr::rename(cities = cities.x,name.com = Comune.x)

xFF = randomVectorField_estimation_viaRSS$estimationRVF$y_Rel[,2]
WxFF = randomVectorField_estimation_viaRSS$estimationRVF$Wy_Rel[,2]

my = mean(randomVectorField_estimation_viaRSS$estimationRVF$actObs_y_Rel)
mWy = mean(randomVectorField_estimation_viaRSS$estimationRVF$actObs_Wy_Rel)

df = df %>% bind_cols(xFF = xFF,WxFF = WxFF)
df = df %>% mutate(quadrant = ifelse(xFF >= .GlobalEnv$my & WxFF >= .GlobalEnv$mWy, "HH", 
                                     ifelse(xFF >= .GlobalEnv$my & WxFF <= .GlobalEnv$mWy, "HL",
                                            ifelse(xFF <= .GlobalEnv$my & WxFF <= .GlobalEnv$mWy, "LL","LH"))))
df = df %>% mutate(colPlot = ifelse(quadrant %in% c("HH","HL"),"UA", 
                                    ifelse(quadrant == "LL","RA","SA")))
df = df %>% mutate(colPlotPlot = ifelse(quadrant %in% c("HH","HL"),"red", 
                                        ifelse(quadrant == "LL","forestgreen","yellow")))  
dfAttractors =  df
save(dfAttractors, file="dfAttractors.RData")

# precomputed
# load(file="../datasets_it/dfAttractors.RData")



# For each municipalities probability of being in a specific attractor
df = dfAttractors 
dfBoot = df
centerUx = mean(filter(df,colPlot=="UA")$xFF)
centerUWx = mean(filter(df,colPlot=="UA")$WxFF)
centerSAx = mean(filter(df,colPlot=="SA")$xFF)
centerSAWx = mean(filter(df,colPlot=="SA")$WxFF)
centerRx = mean(filter(df,colPlot=="RA")$xFF)
centerRWx = mean(filter(df,colPlot=="RA")$WxFF)
my = mean(randomVectorField_estimation_viaRSS$estimationRVF$actObs_y_Rel)
mWy = mean(randomVectorField_estimation_viaRSS$estimationRVF$actObs_Wy_Rel)

dfBootAttr = matrix(NA, nrow=nrow(dfBoot), ncol=nBoot)
ii = (xMatBoot >= my & WxMatBoot >= mWy) 
dfBootAttr[ii] = "UA"

ii = (xMatBoot >= my & WxMatBoot <= mWy) 
dfBootAttr[ii] = "UA"

ii = (xMatBoot <= my & WxMatBoot <= mWy) 
dfBootAttr[ii] = "RA"

ii = (xMatBoot <= my & WxMatBoot >= mWy) 
dfBootAttr[ii] = "SA"

# percentage in every boostrap
percAttractors = matrix(NA, nrow=nrow(dfBootAttr), ncol=3)
colnames(percAttractors) = c("UAPerc", "SAPerc", "RAPerc")

for (i in 1:nrow(dfBootAttr)){
    percAttractors[i,colnames(percAttractors)=="UAPerc"] = sum(dfBootAttr[i,]=="UA")
    percAttractors[i,colnames(percAttractors)=="SAPerc"] = sum(dfBootAttr[i,]=="SA")
    percAttractors[i,colnames(percAttractors)=="RAPerc"] = sum(dfBootAttr[i,]=="RA")
}
percAttractors = percAttractors/100
df = df %>% bind_cols(percAttractors)

redRGB = col2rgb("red", alpha = FALSE)
yellowRGB = col2rgb("yellow", alpha = FALSE)
forestGreenRGB = col2rgb("forestgreen", alpha = FALSE)

for (i in 1:nrow(df)){
    df$rgbCol[i] = rgbCol = do.call("rgb",c(as.list(df$UAPerc[i]*redRGB + df$SAPerc[i]*yellowRGB + df$RAPerc[i]*forestGreenRGB), maxColorValue = 255))
}

# save(df,file="dfAttractorsFuzzy.RData")
# precomputed 
load("../datasets_it/dfAttractorsFuzzy.RData")

# Figure 5 (a) ----
dev.new()
randomVectorField_estimation <- randomVectorField_estimation_viaRSS$estimationRVF

# moranAct <- lm(df$wxAct ~ df$xAct)
moranAct.sm <- sm.regression(df$xAct,df$wxAct,se=TRUE,display="none")
moranFut.sm <- sm.regression(df$xFut,df$wxFut,se=TRUE,display="none")

#Where we can calculate a directional vectors
where.est.directions <- randomVectorField_estimation$stacked_expected_Delta_y!= 0 | randomVectorField_estimation$stacked_expected_Delta_Wy!= 0
where.est.directions[where.est.directions==0] <- NA
where.est.directions[where.est.directions==TRUE] <- 1

arrowLength = sqrt(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot^2 + bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot^2)
arrowAngles = atan(bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot/arrowLength + bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot)
arrowAngles[arrowAngles == "NaN"] <- NA
sdArrowAngles = apply(arrowAngles, MARGIN=1, FUN=sd, na.rm=T)

# standard deviation of angles
image(x=randomVectorField_estimation$grid_y_Rel, y=randomVectorField_estimation$grid_Wy_Rel, z=matrix(sdArrowAngles, ncol=randomVectorField_estimation$numGrid, nrow=randomVectorField_estimation$numGrid,byrow=TRUE),xlab="Population density (log)",ylab="Spatially lagged (LMA) population density (log)",col=hcl.colors(1000,"Blues",rev=TRUE), ylim=c(min(randomVectorField_estimation$grid_Wy_Rel),max(randomVectorField_estimation$grid_Wy_Rel)), xlim=c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)))

#Observations in the moran space
points(df$xFut,df$wxFut,cex=0.25,pch=19,col=df$rgbCol)
abline(h = mean(randomVectorField_estimation$actObs_y_Rel), lty = 2)
abline(v = mean(randomVectorField_estimation$actObs_Wy_Rel), lty = 2)
lines(c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)),c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)))
grid()

lenghtArrows = 0.2

#Estimated RVF at any significance
for (j in 1: nrow(randomVectorField_estimation$evalPoints_Stack)){
    if ((randomVectorField_estimation$stacked_expected_Delta_y[j]!=0) |  (randomVectorField_estimation$stacked_expected_Delta_Wy[j]!=0) ){
        
        arrows(randomVectorField_estimation$evalPoints_Stack[j,1],
               randomVectorField_estimation$evalPoints_Stack[j,2],
               x1=randomVectorField_estimation$evalPoints_Stack[j,1]+ randomVectorField_estimation$stacked_expected_Delta_y[j]*(lenghtArrows),
               y1=randomVectorField_estimation$evalPoints_Stack[j,2] + randomVectorField_estimation$stacked_expected_Delta_Wy[j]*(lenghtArrows),
               code=2,length=0.06,lwd=1.2, col="#CC4600")
    }
}

points(centerUx,centerUWx,cex=5,col="red",lwd=5)
points(centerSAx,centerSAWx,cex=5,col="yellow",lwd=5)
points(centerRx,centerRWx,cex=5,col="forestgreen",lwd=5)

# dev.copy2eps(file="RVFAttractorsFuzzy.eps")

# Figure 5 (b) ----
library(dplyr)
library(tmap)
library(rgdal)
load("../datasets_it/W_SLL_2001_2011_2019.RData")
shp_SLL_2011 <- readOGR("../datasets_it/SLL2011_shapefile/SLL2011_wgs84_EURO.shp")
shp_com = shp_comm_2019
shp_com@data = shp_com@data %>% dplyr::left_join(df,by=c("PRO_COM"="geo"))


tm_obj <- tm_shape(shp_com) + tm_fill("rgbCol", title="", labels=c("Rural","Suburban","Urban"), palette=c("forestgreen","yellow","red")) + tm_shape(shp_com) + tm_fill("rgbCol", title="") + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tm_shape(shp_SLL_2011) + tm_borders("black", alpha=1, lwd=0.5)
dev.new()
tm_obj
tmap_save(tm = tm_obj,  filename = "MapAttractorsFuzzy.eps")

