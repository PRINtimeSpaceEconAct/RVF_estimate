rm(list = ls())

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
library(dplyr)
library(classInt)
library(raster)
library(rgdal)
library(funrar)
library(doMC)
registerDoMC(cores=6)
library(expm)
library(data.table)

load("../datasets_it/W_SLL_2001_2011_2019.RData")
load("../datasets_it/municipalPopulation_1982_2019.RData")
source("lib/decompositionOfVariance.R")


data.pop.all <- data.pop.all %>% mutate(density=(pop/km2)) %>% mutate(log.density=log(pop/km2)) 
data.pop.all <- data.pop.all %>% filter(pop!=0) %>% mutate(time=as.numeric(as.character(time))) %>% filter(time>1983)
data.pop.all.p <- pdata.frame(data.pop.all, index = c("geo", "time"))
data.pop.all.p <- make.pbalanced(data.pop.all.p, balance.type = "shared.individuals")
data.pop.ts <- data.pop.all.p %>% dplyr::group_by(time) %>% dplyr::summarise(meanPop=mean(pop), meanDensity=mean(density), meanLogDensity=mean(log.density), sumPop=sum(pop), sdDensity=sd(density), CV=sdDensity/meanDensity) 
data.popAll.ts <- data.pop.all %>% dplyr::group_by(time) %>% dplyr::summarise( sumPop=sum(pop)/1000000) 




#Figure 2 (a) ----
df <- data.pop.all.p
shp1 <- shp_comm_2001@data  %>% dplyr::select(PRO_COM_T, SLL)
df1 <- df %>% filter(as.numeric(as.character(time))<=2001) %>% left_join(shp1, by=c("geo"="PRO_COM_T"))
shp2 <- shp_comm_2011@data  %>% dplyr::select(PRO_COM_T, SLL)
df2 <- df %>% filter(as.numeric(as.character(time))>2001) %>% left_join(shp2, by=c("geo"="PRO_COM_T"))
df <- rbind(df1, df2)

timeSeriesVarDec <- matrix(NA, nrow=3, ncol=length(years))
rownames(timeSeriesVarDec) <- c("within", "between", "varTot")
colnames(timeSeriesVarDec) <- years
timeSeriesMean <- matrix(NA, nrow=1, ncol=length(years))
colnames(timeSeriesMean) <- years
for (t in 1:length(years)){
    data <- df %>% filter(time==years[t]) 
    varDec <- decompositionOfVariance(data=data, var.name="density", group.name="SLL")
    timeSeriesVarDec[1,t] <- varDec$varWithinPerc
    timeSeriesVarDec[2,t] <- varDec$varBetweenPerc
    timeSeriesVarDec[3,t] <- varDec$varTot
    
    timeSeriesMean[,t] <- mean(data$density)
}

dfVarDec <- matrix_to_stack(
    timeSeriesVarDec[1:2,],
    value_col = "value",
    row_to_col = "Variance decomposition:",
    col_to_col = "time"
)

CV <- sqrt(timeSeriesVarDec[3,])/timeSeriesMean


dev.new()
year2001Index = 18
nyears = length(colnames(timeSeriesVarDec))
plot(colnames(timeSeriesVarDec)[1:year2001Index], timeSeriesVarDec[1,1:year2001Index], type="l", col="black", xlim=as.numeric(c(colnames(timeSeriesVarDec)[1],colnames(timeSeriesVarDec)[nyears])),ylim=c(40,60), xaxt="n", xlab="", ylab="")
lines(colnames(timeSeriesVarDec)[(year2001Index+1):nyears], timeSeriesVarDec[1,(year2001Index+1):nyears], type="l", col="black")
points(colnames(timeSeriesVarDec), timeSeriesVarDec[1,], pch=8)
lines(colnames(timeSeriesVarDec)[1:year2001Index], timeSeriesVarDec[2,1:year2001Index],  col="black")
lines(colnames(timeSeriesVarDec)[(1+year2001Index):nyears], timeSeriesVarDec[2,(1+year2001Index):nyears],  col="black")
points(colnames(timeSeriesVarDec), timeSeriesVarDec[2,], pch=21)
grid()
axis(side=1, at= years, labels= years,las=2,cex.axis=0.7)
par(new=T)
plot(colnames(timeSeriesVarDec), CV, type="l", col="black", axes=F, xlab="", ylab="")
points(colnames(timeSeriesVarDec), CV, pch=19)
lab4 <- round(seq(min(CV), max(CV), length=8), digits=2)
axis(4,col="black",lwd=1, at= lab4, labels= lab4)
mtext(side=4,text="Lcjhdsklhjfkldsjfds",line=2)
abline(v=2013,lty=2)
abline(v=1991, lty=3)
abline(v=2001, lty=3)
abline(v=2011, lty=3)
l <- legend("top", cex = 1, legend = c("Variance within (%)","Variance between (%)", "CV"), lty=c(1,1,1), pch=c(8,21,19), col = c("black", "black", "black"), bg="white")
# dev.copy2pdf(file="varDecomposition_CV.pdf")

# Figure 2 (b) ----
load("../datasets_it/data_x_wx_allYears.RData")

iMoran <- vector("numeric", length = length(years))
seiMoran <- vector("numeric", length = length(years))
for (t in 1:length(years)){
    df <- df_allYears %>% filter(time==years[t])
    iMoran[t] <- coef(lm(WpopDensity ~ popDensity, data=df))[2]
    seiMoran[t] <- sqrt(diag(vcov(lm(WpopDensity ~ popDensity, data=df))))[2]
}

dev.new()
years <- as.numeric(as.character(data.pop.ts$time))
year2001Index = 18
nyears = length(years)

plot(years, data.pop.ts$meanDensity, type="l", ylab="", lwd=2, xaxt="n", xlab="")
grid()
axis(side=1, at= years, labels= years,las=2,cex.axis=0.7)
par(new=T)
plot(years[1:year2001Index], iMoran[1:year2001Index], type="l", ylab="", lwd=2, xaxt="n",col="black",axes=F, lty="dashed", xlab="", ylim=c(min(iMoran-1.96*seiMoran), max(iMoran+1.96*seiMoran)),xlim=c(years[1],years[nyears]))
lines(years[(year2001Index+1):nyears], iMoran[(year2001Index+1):nyears], lwd=2,col="black", lty="dashed")

lines(years[1:year2001Index], iMoran[1:year2001Index]-1.96*seiMoran[1:year2001Index], col="black", lty="dotted")
lines(years[(year2001Index+1):nyears], iMoran[(year2001Index+1):nyears]-1.96*seiMoran[(year2001Index+1):nyears], col="black", lty="dotted")

lines(years[1:year2001Index], iMoran[1:year2001Index]+1.96*seiMoran[1:year2001Index], col="black", lty="dotted")
lines(years[(year2001Index+1):nyears], iMoran[(year2001Index+1):nyears]+1.96*seiMoran[(year2001Index+1):nyears], col="black", lty="dotted")

lab4 <- round(seq(min(iMoran), max(iMoran), length=8), digits=2)
axis(4,col="black",lwd=1, at= lab4, labels= lab4)
mtext(side=4,text="Lcjhdsklhjfkldsjfds",line=2)
abline(v=2013,lty=2)
abline(v=1991, lty=3)
abline(v=2001, lty=3)
abline(v=2011, lty=3)
l <- legend("topleft", cex = 1, legend = c("Mean density","Moran's I"), lty=c(1,2), lwd=c(2,2), col = c("black", "black"), bg="white")
# dev.copy2pdf(file="timeSeriesDensityMeanMoranI.pdf")

#Figure 2 (c) ----
years <- as.numeric(as.character(data.popAll.ts$time))
dev.new()
plot(years, data.popAll.ts$sumPop, type="l", lwd=2, xaxt="n", ylab="Population (million)", xlab="")
axis(side=1, at= years, labels= years,las=2,cex.axis=0.7)
grid()
abline(v=2013,lty=2)
abline(v=1991, lty=3)
abline(v=2001, lty=3)
abline(v=2011, lty=3)
# dev.copy2pdf(file="timeSeriesTotPop.pdf")
