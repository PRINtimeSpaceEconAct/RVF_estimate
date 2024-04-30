plotEstimatedRVF <- function(df, estimationRVF, bootstrapAnalysis=NA, lenghtArrows=1, df_all = NA,munics = NA,drawAll=F){
    
    randomVectorField_estimation <- estimationRVF
   
    # moranAct <- lm(df$wxAct ~ df$xAct)
    moranAct.sm <- sm.regression(df$xAct,df$wxAct,se=TRUE,display="none")
    moranFut.sm <- sm.regression(df$xFut,df$wxFut,se=TRUE,display="none")
    
    #Where we can calculate a directional vectors
    where.est.directions <- randomVectorField_estimation$stacked_expected_Delta_y!= 0 | randomVectorField_estimation$stacked_expected_Delta_Wy!= 0
    where.est.directions[where.est.directions==0] <- NA
    where.est.directions[where.est.directions==TRUE] <- 1
    
    dev.new()
    
    if (!any(is.na(bootstrapAnalysis))){
        #Test of statistical significance of Delta y
        nboot <- ncol(bootstrapAnalysis$stacked_expected_Delta_y_Averaging_boot)
        #Only the directions belonging to (1-pValue)% confidence interval
        pValue <- 0.01
        
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
        
        #dev.new()
        # standard deviation of angles
        image(x=randomVectorField_estimation$grid_y_Rel, y=randomVectorField_estimation$grid_Wy_Rel, z=matrix(sdArrowAngles, ncol=randomVectorField_estimation$numGrid, nrow=randomVectorField_estimation$numGrid,byrow=TRUE),xlab="Population density (log)",ylab="Spatially lagged (LMA) population density (log)",col=hcl.colors(1000,"Blues",rev=TRUE), ylim=c(min(randomVectorField_estimation$grid_Wy_Rel),max(randomVectorField_estimation$grid_Wy_Rel)), xlim=c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)))
        
        #Observations in the moran space
        points(randomVectorField_estimation$y_Rel[,1],randomVectorField_estimation$Wy_Rel[,1],cex=0.25,pch=19,col="yellow")
        points(randomVectorField_estimation$y_Rel[,ncol(randomVectorField_estimation$y_Rel)],randomVectorField_estimation$Wy_Rel[,ncol(randomVectorField_estimation$Wy_Rel)],cex=0.25,pch=19,col="green")
        abline(h = mean(randomVectorField_estimation$actObs_y_Rel), lty = 2)
        abline(v = mean(randomVectorField_estimation$actObs_Wy_Rel), lty = 2)
        lines(c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)),c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)))
        grid()
        
        #Estimated RVF significant at 5% level
        for (j in 1: nrow(randomVectorField_estimation$evalPoints_Stack)){
            if (drawAll == T){
                arrows(randomVectorField_estimation$evalPoints_Stack[j,1], randomVectorField_estimation$evalPoints_Stack[j,2], x1=randomVectorField_estimation$evalPoints_Stack[j,1]+ randomVectorField_estimation$stacked_expected_Delta_y[j]*(lenghtArrows), y1=randomVectorField_estimation$evalPoints_Stack[j,2] +  randomVectorField_estimation$stacked_expected_Delta_Wy[j]*(lenghtArrows),code=2,length=0.06,lwd=1.2, col="red")
            }
            else{
                if ((stacked_expected_Delta_y_stat_sign[j]!=0) |  (stacked_expected_Delta_Wy_stat_sign[j]!=0) ){
                    
                    arrows(randomVectorField_estimation$evalPoints_Stack[j,1], randomVectorField_estimation$evalPoints_Stack[j,2], x1=randomVectorField_estimation$evalPoints_Stack[j,1]+ randomVectorField_estimation$stacked_expected_Delta_y[j]*(lenghtArrows), y1=randomVectorField_estimation$evalPoints_Stack[j,2] +  randomVectorField_estimation$stacked_expected_Delta_Wy[j]*(lenghtArrows),code=2,length=0.06,lwd=1.2, col="red")
                }
            }
        }
        
    }else{
        image(x=randomVectorField_estimation$grid_y_Rel, y=randomVectorField_estimation$grid_Wy_Rel, z=matrix(where.est.directions, ncol=randomVectorField_estimation$numGrid, nrow=randomVectorField_estimation$numGrid,byrow=TRUE),xlab="",ylab="",col="gray93", ylim=c(min(randomVectorField_estimation$grid_Wy_Rel),max(randomVectorField_estimation$grid_Wy_Rel)), xlim=c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)))
        
        #Observations in the moran space
        points(randomVectorField_estimation$y_Rel[,1],randomVectorField_estimation$Wy_Rel[,1],cex=0.25,pch=19,col="yellow")
        points(randomVectorField_estimation$y_Rel[,ncol(randomVectorField_estimation$y_Rel)],randomVectorField_estimation$Wy_Rel[,ncol(randomVectorField_estimation$Wy_Rel)],cex=0.25,pch=19,col="green")
        abline(h = mean(randomVectorField_estimation$actObs_y_Rel), lty = 2)
        abline(v = mean(randomVectorField_estimation$actObs_Wy_Rel), lty = 2)
        lines(c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)),c(min(randomVectorField_estimation$grid_y_Rel),max(randomVectorField_estimation$grid_y_Rel)))
        grid()
        
        #Estimated RVF at any significance
        for (j in 1: nrow(randomVectorField_estimation$evalPoints_Stack)){
            if ((randomVectorField_estimation$stacked_expected_Delta_y[j]!=0) |  (randomVectorField_estimation$stacked_expected_Delta_Wy[j]!=0) ){
                
                arrows(randomVectorField_estimation$evalPoints_Stack[j,1], randomVectorField_estimation$evalPoints_Stack[j,2], x1=randomVectorField_estimation$evalPoints_Stack[j,1]+ randomVectorField_estimation$stacked_expected_Delta_y[j]*(lenghtArrows), y1=randomVectorField_estimation$evalPoints_Stack[j,2] + randomVectorField_estimation$stacked_expected_Delta_Wy[j]*(lenghtArrows),code=2,length=0.06,lwd=1.2, col="red")
            }
        }
    }
    
    lines(moranAct.sm$eval.points,moranAct.sm$estimate, col="blue")
    lines(moranAct.sm$eval.points,moranAct.sm$estimate+2*moranAct.sm$se, lty = 2, col="blue")
    lines(moranAct.sm$eval.points,moranAct.sm$estimate-2*moranAct.sm$se, lty = 2, col="blue")
    
    lines(moranFut.sm$eval.points,moranFut.sm$estimate, col="purple")
    lines(moranFut.sm$eval.points,moranFut.sm$estimate+2*moranFut.sm$se, lty = 2, col="purple")
    lines(moranFut.sm$eval.points,moranFut.sm$estimate-2*moranFut.sm$se, lty = 2, col="purple")
    
    df <- df %>% mutate(name.com=ifelse(name.com=="Milano", "Milan", 
                                        ifelse(name.com=="Firenze", "Florence",
                                               ifelse(name.com=="Torino", "Turin",
                                                      ifelse(name.com=="Venezia", "Venice",
                                                             ifelse(name.com=="Napoli", "Naples",
                                                                    ifelse(name.com=="Genova", "Genoa",
                                                                           ifelse(name.com=="Roma", "Rome",name.com))))))))
    # to draw metropolitan cities
    if (!is.null(df$cities)){
        points(df$xAct[df$cities==1], df$wxAct[df$cities==1], col="orange", pch=19, cex=1)
        points(df$xAct[df$cities==1], df$wxAct[df$cities==1], col="black", pch=21, cex=1)

        points(df$xFut[df$cities==1], df$wxFut[df$cities==1], col="brown", pch=19, cex=1)
        points(df$xFut[df$cities==1], df$wxFut[df$cities==1], col="black", pch=21, cex=1)
        arrows(df$xAct[df$cities==1], df$wxAct[df$cities==1], x1=df$xFut[df$cities==1], y1=df$wxFut[df$cities==1],code=2,length=0,lwd=1.2, col="black")
        text(df$xAct[df$cities==1], df$wxAct[df$cities==1], labels=df$name.com[df$cities==1], pos=3, cex=0.8, col="black")
    }


}

