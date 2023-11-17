rm(list = ls())


library(dplyr)
library(sm)
library(RColorBrewer)
library(classInt)
library(funrar)
library(wordspace)
library(pracma)
library(doMC)
registerDoMC(cores=10)



#Functions
source("lib/RVF_besth_adaptive_bootstrap.R")
source("lib/estimationVectorField_best_h.R")
source("lib/estimationVectorField_best_h_best_alpha.R")
source("lib/bootstrapStandardErrorsAdaptive_RVF_parallel.R")
source("lib/estimationVectorFieldAdaptive.R")
source("lib/epaKernelBivAdaptive.R")
source("lib/normKernelBiv_vec.R")
source("lib/epaKernelBiv_vec.R")
source("lib/forward_estimation_RSS.R")

source("lib/plotEstimatedRVF.R")



#-----------------#
#   1984-2019     #
#-----------------#
load("../datasets_it/data_x_wx_allYears.RData")
df.1 <- df_allYears[df_allYears$time == 1984,]
df.2 <- df_allYears[df_allYears$time == 2019,]
df.1 <- df.1 %>% dplyr::rename(xAct=logPopDensity, wxAct=WlogPopDensity)
df.2 <- df.2 %>% dplyr::rename(xFut=logPopDensity, wxFut=WlogPopDensity)
df <- df.1 %>% dplyr::inner_join(df.2, by="geo") %>% dplyr::select(c(geo, xAct, wxAct, xFut, wxFut,cities.x,Comune.x)) 
df <- df %>% dplyr::rename(cities = cities.x,name.com = Comune.x)
df <- df %>% dplyr::filter((wxAct != 0) & (wxFut != 0))

# estimate RVF
f_name = "1984-2019"
RVF <- RVF_besth_adaptive_bootstrap(df=df, numGrid=50, numGrid.h=10,
                                    numGrid.alpha=10, title=f_name,
                                    nboot=100, pValue=0.05, figure=FALSE,
                                    bestAlpha=TRUE,continuousForward = TRUE)

# load precomputed estimates
load(paste("../datasets_it/randomVectorField_estimation_viaRSS_alphaOpt_TRUE",f_name,".RData",sep=""))
load(paste("../datasets_it/bootstrapAnalysis",f_name,".RData",sep=""))

# Figure 4 ----
# plot estimated RVF
dev.new()
plotEstimatedRVF(estimationRVF=randomVectorField_estimation_viaRSS$estimationRVF,df=df,lenghtArrows = 0.2,bootstrapAnalysis=bootstrapAnalysis,df_all = df_allYears)
