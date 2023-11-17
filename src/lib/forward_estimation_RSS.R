forward_estimation_RSS <- function(estimationRVF,continuousForward=TRUE,
                                   distributionFutDensity=NA){
    if (continuousForward == TRUE){
        return(forwardEstimationContinuousInterpolation(estimationRVF,distributionFutDensity))
    } else{
        return(forwardEstimationOneStep(estimationRVF))
    }
}

forwardEstimationOneStep <- function(estimationRVF){
    #--------------------------------------------------------#
    #               Angela Parenti, Cristiano Ricci          #
    #  forward estimation single step without interpolation  # 
    #--------------------------------------------------------#
    
    y_Rel_act <- estimationRVF$y_Rel[,1]
    y_Rel_fut <- estimationRVF$y_Rel[,2]
    
    Wy_Rel_act <- estimationRVF$Wy_Rel[,1]
    Wy_Rel_fut <- estimationRVF$Wy_Rel[,2]
    
    weiNorm <- estimationRVF$weiNorm
    
    stacked_expected_Delta_y  <- estimationRVF$stacked_expected_Delta_y 
    stacked_expected_Delta_Wy <- estimationRVF$stacked_expected_Delta_Wy
    
    nObs = length(y_Rel_act)
    nEval = length(stacked_expected_Delta_y)
    
    evalPoints <- estimationRVF$evalPoints_Stack     # (ngrid*ngrid) x 2
    
    obs_act = cbind(y_Rel_act,Wy_Rel_act)           # nObs x 2
    obs_fut = cbind(y_Rel_fut,Wy_Rel_fut)           # nObs x 2
    obs_pred = matrix(NA, nrow = nObs, ncol = 2)    # nObs x 2
    
    arrows = cbind(stacked_expected_Delta_y,stacked_expected_Delta_Wy)  
    
    for (i in 1:nObs){
        closestEval_j = which.min( rowNorms( t(obs_act[i,] - t(evalPoints)) ))
        obs_pred[i,] = obs_act[i,] + arrows[closestEval_j,]
    }
    
    RSSall = rowNorms(obs_fut - obs_pred)
    RSS = sum(weiNorm*rowNorms(obs_fut - obs_pred))
    
    list(RSS=RSS,RSSall=RSSall,obs_pred=obs_pred)
}

forwardEstimationContinuousInterpolation <- function(estimationRVF,distributionFutDensity=NA){
    #--------------------------------------------------------#
    #               Angela Parenti, Cristiano Ricci          #
    #           Forward estimation continuous time step      #
    #   with linear interpolation between evaluation Points  #
    #--------------------------------------------------------#
    
    # observations
    y_Rel_act <- estimationRVF$y_Rel[,1]
    y_Rel_fut <- estimationRVF$y_Rel[,2]
    Wy_Rel_act <- estimationRVF$Wy_Rel[,1]
    Wy_Rel_fut <- estimationRVF$Wy_Rel[,2]
    weiNorm <- estimationRVF$weiNorm
    
    nObs = length(y_Rel_act)
    obs_act = cbind(y_Rel_act,Wy_Rel_act)           # nObs x 2
    obs_fut = cbind(y_Rel_fut,Wy_Rel_fut)           # nObs x 2
    
    # evaluation Points and arrows 
    evalPoints <- estimationRVF$evalPoints_Stack     # (nGrid*nGrid) x 2
    stacked_expected_Delta_y  <- estimationRVF$stacked_expected_Delta_y 
    stacked_expected_Delta_Wy <- estimationRVF$stacked_expected_Delta_Wy
    arrows = cbind(stacked_expected_Delta_y,stacked_expected_Delta_Wy)  
    nGrid = sqrt(length(stacked_expected_Delta_y))
    
    # choosing dt depending on longest arrow and grid resolution
    # assuming cartesian grid of eval Points (dx possibly != dy)
    maxArrLen = max(rowNorms(arrows))       # maximum length of arrows
    evalPoint_y <- unique(evalPoints[,1])
    evalPoint_Wy <- unique(evalPoints[,2])
    dyGrid = sort(evalPoint_y)[2]-sort(evalPoint_y)[1]      # grid dist on y
    dWyGrid = sort(evalPoint_Wy)[2]-sort(evalPoint_Wy)[1]   # grid dist on Wy
    dMin = min(dyGrid,dWyGrid)  
    
    # maxArrLen/dMin = max number of cells done with one 'unitary' step
    # dMin/maxArrLen = fraction of a step to follow an arrow of length maxArrLen
    #                   without skipping one cell
    # dMin/maxArrLen / 10 = for extra smoothness, 10 steps per cells 
    # nSteps = 1/ ( dMin/maxArrLen / 10 ) # number of steps, rounded from above
    nSteps = ceiling(1/( dMin/maxArrLen / 10 ))
    dt = 1/nSteps
    
    # time iteration     
    obs_pred = obs_act          # used as temporary storage for time iteration
    for (t in 1:nSteps){
        Interp_y = interp2(evalPoint_y,evalPoint_Wy,
                           Reshape(stacked_expected_Delta_y,nGrid,nGrid),
                           obs_pred[,1],obs_pred[,2])
        Interp_Wy = interp2(evalPoint_y,evalPoint_Wy,
                            Reshape(stacked_expected_Delta_Wy,nGrid,nGrid),
                            obs_pred[,1],obs_pred[,2])
        
        obs_pred = obs_pred + dt * cbind(Interp_y,Interp_Wy)
    }
    
    RSSall = rowNorms(obs_fut - obs_pred)
    RSS = sum(weiNorm*rowNorms(obs_fut - obs_pred))
    return(list(RSS=RSS,RSSall=RSSall,obs_pred=obs_pred))
    
    # distributionPredDensity <- rowSums(epaKernelBivAdaptive(evalPoints,obs_pred))
    # Area = dyGrid*dWyGrid
    # RSS = 1 - (1/Area)*sum(sqrt(distributionPredDensity*distributionFutDensity))
    # return(list(RSS=RSS,obs_pred=obs_pred))
}
