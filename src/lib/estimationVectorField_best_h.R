estimationVectorField_best_h <- function(numGrid, grid_y_Rel, grid_Wy_Rel, y_Rel,
                                         Wy_Rel, W,numLag,numGrid.h, nrepsimPar=100,
                                         min.size.drawnPar=1,
                                         weiNorm=rep(1/length(y_Rel),length(y_Rel)),
                                         continuousForward=FALSE,adaptive=TRUE){
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    #Selection of the best bandwidth h using the in-sample RSS #
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#  
    
    #INPUTS	
    #numGrid: number of grid points for the estimate of vector field (numGrid x numGrid is the total number)
    #grid_y_Rel: grid for relative y
    #grid_Wy_Rel: grid for relative Wy
    #y_Rel: y in relative terms
    #Wy_Rel: W*y_Rel
    #numLag: number of lag to calculate Delta y and Delta Wy
    #numGrid.h: number of grid for searching the optimal bandwidth
    #nrepsim: number of replication of the simulation of the last-year distribution
    #min.size.draw: minimum size of the sampling of the transitions for the last-year distribution
    
    #OUTPUTS
    #estimationRVF: estimatin results for model with the optimal h.bandwidth via min RSS
    #RSS_for_best_h: residual sum of squares of all model estimated
    #minRSS: min RSS
    #hOpt: optimal bandwidth associated with the minimum RSS
    
    
    #Evaluation points for the stochastic kernel in 2-column format
    evalPoints_Stack <- matrix(NA,nrow=numGrid*numGrid,ncol=2)
    for (i in 1:numGrid){
        evalPoints_Stack[(1+(i-1)*numGrid):(i*numGrid),] <- cbind(rep(grid_y_Rel[i],numGrid),grid_Wy_Rel)
    }
    
    
    #------------------------------------# 
    # Number of lags for the transitions #
    #------------------------------------#          
    N <- nrow(y_Rel)
    numTrans <- (dim(y_Rel)[2]-numLag)
    NT <- N*numTrans
    
    #----------------------------------------# 
    # Stacked actual and future observations #
    #----------------------------------------# 
    futObs_y_Rel <- c(y_Rel[,(1+numLag):dim(y_Rel)[2]])
    actObs_y_Rel <- c(y_Rel[,1:(dim(y_Rel)[2]-numLag)])
    
    futObs_Wy_Rel <- c(Wy_Rel[,(1+numLag):dim(Wy_Rel)[2]])
    actObs_Wy_Rel <- c(Wy_Rel[,1:(dim(Wy_Rel)[2]-numLag)])
    
    #Delta y and Wy
    delta_y_Rel <- futObs_y_Rel-actObs_y_Rel
    delta_Wy_Rel <- futObs_Wy_Rel-actObs_Wy_Rel
    
    #Matrix of observations
    x <- cbind(actObs_y_Rel,actObs_Wy_Rel)
    
    #Inverse of covariance matrix
    invS <- solve(cov(x))
    #Determinant of covariance matrix
    detS <- det(cov(x))
    
    #Initilization of the bandwidth (use the Epanechnikov bandwidth as in Silverman p.87)
    h.initial <- c(1.77*nrow(x)^(-1/6),1.77*nrow(x)^(-1/6))
    # #First estimation with h.initial (the one of the normal kernel)
    # estimateDensity_Stack <- foreach(j=1:(numGrid* numGrid), .combine= cbind, .verbose=FALSE) %dopar%{   
    # epaKernelBiv(y=evalPoints_Stack[j,],x=x,invS=invS,detS=detS,h.bandwidth=h.initial)  
    # }
    
    #Grid for h.bandwidth as Bowmanm and Azzalini (function "hcv" for sm)
    hstart <- h.initial/10 
    hend <- h.initial*2
    hgrid <- log(hstart) + (log(hend) - log(hstart)) * (0:(numGrid.h-1))/(numGrid.h-1)
    hgrid <- exp(hgrid)
    
    #Use the same h.bandwidth for both dimensions of x
    hGrid <- cbind(hgrid,hgrid)
    
    #RSS_for_best_h <- vector("numeric",length=nrow(numGrid.h))
    RSS_for_best_h <- vector("numeric",length=(numGrid.h))
    
    
    # compute density estimation for future observation
    distributionFutDensity <- rowSums(epaKernelBivAdaptive(evalPoints_Stack,cbind(futObs_y_Rel,futObs_Wy_Rel)))
    
    RSS_for_best_h <- foreach(i=1:(numGrid.h), .combine= cbind, .verbose=FALSE) %dopar%{
        
        print(paste("testing h = ",hGrid[i,1],", i = ",i,"/",numGrid.h,sep=""))
        
        #Estimation of density
        estimationRVF <-  estimationVectorFieldAdaptive(numGrid=numGrid, 
                                                        grid_y_Rel= grid_y_Rel,
                                                        grid_Wy_Rel= grid_Wy_Rel, 
                                                        y_Rel=y_Rel, 
                                                        Wy_Rel=Wy_Rel, 
                                                        numLag=numLag, 
                                                        figure =FALSE, 
                                                        h.bandwidth=hGrid[i,], 
                                                        weiNorm=weiNorm,
                                                        adaptive=adaptive)
        
        # new method based directly on mean of arrows
        simulation_RSS <- forward_estimation_RSS(estimationRVF,
                                                 continuousForward=continuousForward,
                                                 distributionFutDensity=distributionFutDensity)
        
        # RSS_for_best_h[i] <- simulation_RSS$RSS 		
        return(simulation_RSS$RSS)
    }
    
    minRSS <- min(RSS_for_best_h) 
    positionMin <- which.min(RSS_for_best_h)
    hOpt <- hGrid[positionMin,]
    
    print(paste("finished choice of h. best h = (",hOpt[1],",",hOpt[2],")",sep = ""))
    
    ################################################
    # Estimation for the hOpt selected via min RSS #
    ################################################
    print("starting final adaptive estimation using best h")
    estimationRVF <-  estimationVectorFieldAdaptive(numGrid=numGrid, 
                                                    grid_y_Rel= grid_y_Rel,
                                                    grid_Wy_Rel= grid_Wy_Rel, y_Rel=y_Rel, 
                                                    Wy_Rel=Wy_Rel, numLag=numLag, 
                                                    figure =FALSE, 
                                                    h.bandwidth=hOpt,
                                                    weiNorm=weiNorm,
                                                    adaptive=adaptive)

    
    # redo prediction with optimal h and adaptive estimation
    simulation_RSS <- forward_estimation_RSS(estimationRVF,continuousForward=continuousForward,distributionFutDensity=distributionFutDensity)
    estimatedLastYearDistr <- simulation_RSS$obs_pred
    list(estimationRVF= estimationRVF, RSS_for_best_h= RSS_for_best_h, minRSS= minRSS, hOpt= hOpt, xObs=x, estimatedLastYearDistr= estimatedLastYearDistr, alphaOpt=0.5)
}