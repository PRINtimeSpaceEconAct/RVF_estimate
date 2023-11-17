require(sm)
require(RColorBrewer)

RVF_besth_adaptive_bootstrap <- function(df, numGrid=50, numGrid.h=10, 
                                         numGrid.alpha=10, title="",
                                         nboot=100, pValue=0.05, figure=FALSE, 
                                         weiNorm=rep(1/nrow(df),nrow(df)),
                                         bestAlpha=FALSE,continuousForward=FALSE){
    
    
    
    y <- cbind(df$xAct, df$xFut)
    Wy <- cbind(df$wxAct, df$wxFut)
    
    #Grid for relative y
    grid_y <- seq((min(y)-1),(max(y)+1),length=numGrid)
    #Grid for relative Wy
    grid_Wy <- seq((min(Wy)-1),(max(Wy)+1),length=numGrid)
    
    #*******************#
    # Estimation best h #
    #*******************#
    if (bestAlpha==TRUE){
      print("Choosing best h and alpha")
      randomVectorField_estimation_viaRSS <- estimationVectorField_best_h_best_alpha(
                                            numGrid=numGrid, 
                                            grid_y_Rel= grid_y, 
                                            grid_Wy_Rel= grid_Wy, 
                                            y_Rel=y, Wy_Rel=Wy, 
                                            numLag=1, 
                                            numGrid.h=numGrid.h,
                                            numGrid.alpha=numGrid.alpha,
                                            nrepsimPar=100, 
                                            min.size.drawnPar=1,
                                            weiNorm=weiNorm,
                                            continuousForward=continuousForward)
    }else{
      print("Choosing best h")
      randomVectorField_estimation_viaRSS <- estimationVectorField_best_h(
                                            numGrid=numGrid, 
                                            grid_y_Rel= grid_y, 
                                            grid_Wy_Rel= grid_Wy, 
                                            y_Rel=y, Wy_Rel=Wy, 
                                            numLag=1, 
                                            numGrid.h=numGrid.h,
                                            nrepsimPar=100, 
                                            min.size.drawnPar=1,
                                            weiNorm=weiNorm,
                                            continuousForward=continuousForward)
    }
    
    save(randomVectorField_estimation_viaRSS,
         file=paste("randomVectorField_estimation_viaRSS_alphaOpt_",bestAlpha, 
                    title, ".RData",sep=""))
    
    #***********#
    # BOOTSTRAP #
    #***********#
    estimationRVFList <- vector("list", length=1)
    estimationRVFList[[1]]$estimationRVF <- randomVectorField_estimation_viaRSS$estimationRVF
    randomVectorField_estimation <- randomVectorField_estimation_viaRSS$estimationRVF
    
    print("Bootstrap")
    bootstrapAnalysis <- bootstrapStandardErrorsAdaptive_RVF_parallel(
      stacked_expected_Delta_y=randomVectorField_estimation$stacked_expected_Delta_y,
      stacked_expected_Delta_Wy=randomVectorField_estimation$stacked_expected_Delta_Wy,
      estimationRVFList,optimalBandwith=randomVectorField_estimation_viaRSS$hOpt,
        nboot=nboot, weiNorm=weiNorm, alpha=randomVectorField_estimation_viaRSS$alphaOpt)
    
    save(bootstrapAnalysis, 
         file=paste("bootstrapAnalysis", title, ".RData", sep=""))
    

    #********************************#
    # Inference on RVF via bootstrap #
    #********************************#
    stacked_expected_Delta_y_boot_ConfInt <- vector(mode="numeric",length=nrow(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot))	
    stacked_expected_Delta_Wy_boot_ConfInt <- vector(mode="numeric",length=nrow(bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot))
    for (j in 1:nrow(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot)){	  
        stacked_expected_Delta_y_boot_ConfInt[j] <-  max(sum(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot[j,]<0),sum(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot[j,]>0))/nboot
        
        stacked_expected_Delta_Wy_boot_ConfInt[j] <-  max(sum(bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot[j,]<0),sum(bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot[j,]>0))/nboot
    }
    
    stacked_expected_Delta_y_stat_sign <- randomVectorField_estimation$stacked_expected_Delta_y
    stacked_expected_Delta_y_stat_sign[stacked_expected_Delta_y_boot_ConfInt<(1-pValue)] <-0
    
    stacked_expected_Delta_Wy_stat_sign <- randomVectorField_estimation$stacked_expected_Delta_Wy
    stacked_expected_Delta_Wy_stat_sign[stacked_expected_Delta_Wy_boot_ConfInt<(1-pValue)] <-0
    
    # variance of angle and length 
    arrowLength = sqrt(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot^2 + 
                bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot^2)
    arrowAngles = atan(bootstrapAnalysis$stacked_expected_Delta_Wy_Averaging_boot/
                           arrowLength + bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot)
    arrowAngles[arrowAngles == "NaN"] <- NA
    sdArrowAngles = apply(arrowAngles, MARGIN=1, FUN=sd, na.rm=T)
    
    #*******#
    # Moran #
    #*******#
    moranAct.sm <- sm.regression(df$xAct,df$wxAct,se=TRUE,display="none")
    moranFut.sm <- sm.regression(df$xFut,df$wxFut,se=TRUE,display="none")
    
    
    #*******#
    # Plot  #
    #*******#
    if (figure==FALSE){
        return(list(randomVectorField_estimation_viaRSS=randomVectorField_estimation_viaRSS, 
             bootstrapAnalysis=bootstrapAnalysis, 
             stacked_expected_Delta_y_stat_sign=stacked_expected_Delta_y_stat_sign, 
             stacked_expected_Delta_Wy_stat_sign=stacked_expected_Delta_Wy_stat_sign, 
             moranAct.sm=moranAct.sm, moranFut.sm=moranFut.sm))
    }
    
    

    
    dev.new()
    
    #Area of estimation
    image(x=randomVectorField_estimation$grid_y_Rel, y=randomVectorField_estimation$grid_Wy_Rel, z=matrix(sdArrowAngles, ncol=randomVectorField_estimation$numGrid, nrow=randomVectorField_estimation$numGrid,byrow=TRUE),xlab="",ylab="",col=hcl.colors(1000,"Blues",rev=TRUE), ylim=c(0,10), xlim=c(0,10))
    
    
    #Observations in the moran space
    abline(h = mean(randomVectorField_estimation$actObs_y_Rel), lty = 2)
    abline(v = mean(randomVectorField_estimation$actObs_Wy_Rel), lty = 2)
    lines(c(-2,13),c(-2,13))
    grid()
    
    #Estimated RVF significant at 5% level  
    for (j in 1: nrow(randomVectorField_estimation$evalPoints_Stack)){
        if ((stacked_expected_Delta_y_stat_sign[j]!=0) |  (stacked_expected_Delta_Wy_stat_sign[j]!=0) ){
            
            arrows(randomVectorField_estimation$evalPoints_Stack[j,1], randomVectorField_estimation$evalPoints_Stack[j,2], x1=randomVectorField_estimation$evalPoints_Stack[j,1]+ randomVectorField_estimation$stacked_expected_Delta_y[j]*(0.1), y1=randomVectorField_estimation$evalPoints_Stack[j,2] +  randomVectorField_estimation$stacked_expected_Delta_Wy[j]*(0.1),code=2,length=0.06,lwd=1.2, col="red")
        }
    }
    
    lines(moranAct.sm$eval.points,moranAct.sm$estimate, col="blue")
    lines(moranAct.sm$eval.points,moranAct.sm$estimate+2*moranAct.sm$se, lty = 2, col="blue")
    lines(moranAct.sm$eval.points,moranAct.sm$estimate-2*moranAct.sm$se, lty = 2, col="blue")
    
    lines(moranFut.sm$eval.points,moranFut.sm$estimate, col="purple")
    lines(moranFut.sm$eval.points,moranFut.sm$estimate+2*moranFut.sm$se, lty = 2, col="purple")
    lines(moranFut.sm$eval.points,moranFut.sm$estimate-2*moranFut.sm$se, lty = 2, col="purple")
    
    points(df$xAct[df$cities==1], df$wxAct[df$cities==1], col="orange", pch=19, cex=1)
    points(df$xAct[df$cities==1], df$wxAct[df$cities==1], col="black", pch=21, cex=1)
    
    points(df$xFut[df$cities==1], df$wxFut[df$cities==1], col="brown", pch=19, cex=1)
    points(df$xFut[df$cities==1], df$wxFut[df$cities==1], col="black", pch=21, cex=1)
    arrows(df$xAct[df$cities==1], df$wxAct[df$cities==1], x1=df$xFut[df$cities==1], y1=df$wxFut[df$cities==1],code=2,length=0,lwd=1.2, col="black")
    text(df$xAct[df$cities==1], df$wxAct[df$cities==1], labels=df$name.com[df$cities==1], pos=3, cex=0.8, col="black")
    
    
    dev.copy2eps(file=paste("pictureSignificantDynamicsRVF_", title, ".eps", sep=""))
    
    
    list(randomVectorField_estimation_viaRSS=randomVectorField_estimation_viaRSS,          
         bootstrapAnalysis=bootstrapAnalysis, 
         stacked_expected_Delta_y_stat_sign=stacked_expected_Delta_y_stat_sign, 
         stacked_expected_Delta_Wy_stat_sign=stacked_expected_Delta_Wy_stat_sign, 
         moranAct.sm=moranAct.sm, moranFut.sm=moranFut.sm,
         weiNorm=weiNorm)
    
    
    
    
    
    
    
    
}