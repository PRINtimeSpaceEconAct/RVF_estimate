###########################
# RVF ESTIMATION AND PLOT #
###########################
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
#   1985-2020     #
#-----------------#
load("../datasets_it/dataset_x_wx_cells_GHSL.RData")
df = df.cells.sll.notZero %>% dplyr::rename(xAct=logPopDen1985, wxAct=WlogPopDen1985, xFut=logPopDen2020, wxFut=WlogPopDen2020) %>% 
    dplyr::select(xAct,wxAct,xFut,wxFut)


# estimate RVF
f_name = "1985-2020_GHSL"
RVF <- RVF_besth_adaptive_bootstrap(df=df, numGrid=50, numGrid.h=10,
                                    numGrid.alpha=10, title=f_name,
                                    nboot=100, pValue=0.05, figure=FALSE,
                                    bestAlpha=TRUE,continuousForward = TRUE,adaptive=FALSE)
save(df,RVF,file=paste("RVF_",f_name,".RData",sep=""))

# load precomputed estimates
load(paste("../datasets_it/randomVectorField_estimation_viaRSS_alphaOpt_TRUE",f_name,".RData",sep=""))
load(paste("../datasets_it/bootstrapAnalysis",f_name,".RData",sep=""))


# plot estimated RVF
dev.new()
plotEstimatedRVF(estimationRVF=randomVectorField_estimation_viaRSS$estimationRVF,df=df,lenghtArrows = 0.2,bootstrapAnalysis=bootstrapAnalysis,df_all = df_allYears)

# save
dev.copy2pdf(file=paste("RVF_",f_name,".pdf",sep=""),width = 6, height = 6)





