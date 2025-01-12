bootstrapStandardErrorsAdaptive_RVF_parallel <- function(stacked_expected_Delta_y, stacked_expected_Delta_Wy, estimationRVFList, optimalBandwith=c(1.77*nrow(estimationRVFList[[1]]$estimationRVF$y_Rel)^(-1/6),1.77*nrow(estimationRVFList[[1]]$estimationRVF$y_Rel)^(-1/6)), nboot=100,
 weiNorm=rep(1/length(estimationRVFList[[1]]$estimationRVF$y_Rel),length(estimationRVFList[[1]]$estimationRVF$y_Rel)), alpha=0.5,adaptive=TRUE){

  ###############################################
  # Code for "Spatial Clubs in European Regions #
  #                                             #
  #   D. Fiaschi, L. Gianmoena and A. Parenti   #
  #                                             #
  #       Update: October 2016                  #
  ###############################################
  
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
	#Bootstrap Standard Errors of the estimated random vector field #
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
	
	#INPUTS	
	#stacked_expected_Delta_y: estimated direction both for y averaged over the W
	#stacked_expected_Delta_Wy: estimated direction both for Wy averaged over the W
  #estimationRVF: RVF estimation
	#nboot: number of bootstrap
	
	#OUTPUTS
	#estimateDensity_StackAveraging_boot: bootstrap estimated densities for each point in the grid

    #Evaluation points of RVF (equal for each W)
    evalPoints_Stack <- estimationRVFList[[1]]$estimationRVF$evalPoints_Stack
    #Number of grid of RVF (equal for each W)
    numGrid <- estimationRVFList[[1]]$estimationRVF$numGrid
    #Lag of Markovian process (equal for each W)
    numLag <- estimationRVFList[[1]]$estimationRVF$numLag
    grid_y_Rel <- estimationRVFList[[1]]$estimationRVF$grid_y_Rel
    grid_Wy_Rel <- estimationRVFList[[1]]$estimationRVF$grid_Wy_Rel
    

	#Matrix of expected Delta y and Wy
	stacked_expected_Delta_y_boot <- matrix(NA, nrow= numGrid* numGrid, ncol=length(estimationRVFList))
	stacked_expected_Delta_Wy_boot <- matrix(NA, nrow= numGrid* numGrid, ncol=length(estimationRVFList))

	#Matrix of expected Delta y and Wy
	stacked_expected_Delta_y_Averaging_boot <- matrix(NA,nrow=length(stacked_expected_Delta_y),ncol=nboot)
	stacked_expected_Delta_Wy_Averaging_boot <- matrix(NA,nrow=length(stacked_expected_Delta_Wy),ncol=nboot)
  
	combined.Delta_y.Delta_Wy_Averaging_boot <- rbind(stacked_expected_Delta_y_Averaging_boot,
	                                                  stacked_expected_Delta_Wy_Averaging_boot)
  #Bootstrap 
	combined.Delta_y.Delta_Wy_Averaging_boot <- foreach(b=1:nboot, .combine= cbind, .verbose=FALSE) %dopar%{
	  
	  print(paste("Bootstrap number equal to ",b))
	  
	  gc()
	  
	  #Bootstrapped matrix of observations
	  y_Rel <- estimationRVFList[[1]]$estimationRVF$y_Rel	
	  Wy_Rel <- estimationRVFList[[1]]$estimationRVF$Wy_Rel	
	  
	  #Stacked actual and future observations
	  futObs_y_Rel <- c(y_Rel[,(1+numLag):ncol(y_Rel)])
	  actObs_y_Rel <- c(y_Rel[,1:(ncol(y_Rel)-numLag)])	
	  futObs_Wy_Rel <- c(Wy_Rel[,(1+numLag):ncol(Wy_Rel)])
	  actObs_Wy_Rel <- c(Wy_Rel[,1:(ncol(Wy_Rel)-numLag)])
	  
	  #Matrix of observations
	  x.1 <- cbind(actObs_y_Rel,actObs_Wy_Rel)
	  uuu <- sample(1:nrow(x.1),replace=TRUE)
	  
	  for (i in 1:length(estimationRVFList)){
	    print(i)
	    y_Rel <- estimationRVFList[[i]]$estimationRVF$y_Rel	
	    Wy_Rel <- estimationRVFList[[i]]$estimationRVF$Wy_Rel	
	    
	    #Bootstrapped sample
	    y_Rel_boot <- y_Rel[uuu,]
	    Wy_Rel_boot <- Wy_Rel[uuu,]
	    
	    estimationRVF_boot <-  estimationVectorFieldAdaptive(numGrid=numGrid, 
	                                                    grid_y_Rel= grid_y_Rel,
	                                                    grid_Wy_Rel= grid_Wy_Rel, 
	                                                    y_Rel=y_Rel_boot, 
	                                                    Wy_Rel=Wy_Rel_boot, 
	                                                    numLag=numLag, 
	                                                    figure =FALSE, 
	                                                    h.bandwidth=optimalBandwith,
	                                                    alpha=alpha,
	                                                    weiNorm=weiNorm,
	                                                    adaptive=adaptive)
	    
	    # #Stacked actual and future observations
	    # futObs_y_Rel <- c(y_Rel[,(1+numLag):ncol(y_Rel)])
	    # actObs_y_Rel <- c(y_Rel[,1:(ncol(y_Rel)-numLag)])	
	    # futObs_Wy_Rel <- c(Wy_Rel[,(1+numLag):ncol(Wy_Rel)])
	    # actObs_Wy_Rel <- c(Wy_Rel[,1:(ncol(Wy_Rel)-numLag)])
	    # 
	    # #Matrix of observations
	    # x <- cbind(actObs_y_Rel,actObs_Wy_Rel)
	    # #Bootstrap sample		  
	    # x_boot <-  x[uuu,]
	    # weiNorm_boot <- weiNorm[uuu]
	    # weiNorm_boot <- weiNorm_boot/sum(weiNorm_boot)
	    # 
	    # #Bootstrap of Delta y and Delta Wy
	    # delta_y_Rel_boot <- estimationRVFList[[i]]$estimationRVF$delta_y_Rel[uuu]
	    # delta_Wy_Rel_boot <- estimationRVFList[[i]]$estimationRVF$delta_Wy_Rel[uuu]
	    # 
	    # #Determinant of covariance matrix
	    # detS_boot <- det(cov(x_boot))	  
	    # if (abs(detS_boot) < 1e-8){ break }
	    # 
	    # #Inverse of covariance matrix
	    # invS_boot <- solve(cov(x_boot))
	    # #h.bandwidth
	    # h.bandwidth.par <- optimalBandwith
	    # 
	    # estimateDensity_Stack_boot <- epaKernelBivAdaptive(evalPoints_Stack=evalPoints_Stack,x= x_boot,invS= invS_boot,detS= detS_boot, ngrid=numGrid, h.bandwidth=h.bandwidth.par, weiNorm=weiNorm_boot, alpha=alpha)
	    # # estimateDensity_Stack_boot <- matrix(0,nrow=nrow(x_boot),ncol=(numGrid*numGrid))
	    # 
	    # 
	    # #Matrix of variation of y and Wy
	    # delta_y_Rel_matrix_boot <- matrix(delta_y_Rel_boot,nrow=length(estimationRVFList[[i]]$estimationRVF$delta_y_Rel),ncol=ncol(estimationRVFList[[i]]$estimationRVF$estimateDensity_Stack))
	    # 
	    # delta_Wy_Rel_matrix_boot <- matrix(delta_Wy_Rel_boot,nrow=length(estimationRVFList[[i]]$estimationRVF$delta_y_Rel),ncol=ncol(estimationRVFList[[i]]$estimationRVF$estimateDensity_Stack))
	    # 
	    # 
	    # #Normalizing weights for local mean estimator (see Bowman and Azzalini, p. 49, 1997)
	    # normWeights <- colSums(estimateDensity_Stack_boot)
	    # 
	    # #Check for zero weights (0/1=0, instead of 0/0 that is indeterminate but it should be zero)
	    # normWeights[normWeights==0] <-1
	    # stacked_expected_Delta_y_boot[,i] <- colSums(delta_y_Rel_matrix_boot*estimateDensity_Stack_boot)/normWeights
	    # stacked_expected_Delta_Wy_boot[,i] <- colSums(delta_Wy_Rel_matrix_boot* estimateDensity_Stack_boot)/normWeights
	    	  	  
	  }
	  c(as.vector(estimationRVF_boot$stacked_expected_Delta_y),as.vector(estimationRVF_boot$stacked_expected_Delta_Wy))
	}
	
	stacked_expected_Delta_y_Averaging_boot <- combined.Delta_y.Delta_Wy_Averaging_boot[1:(numGrid*numGrid),]
	stacked_expected_Delta_Wy_Averaging_boot <- combined.Delta_y.Delta_Wy_Averaging_boot[((numGrid*numGrid)+1):nrow(combined.Delta_y.Delta_Wy_Averaging_boot),]

	# #Check if bootstrap is well-done
	# plot(stacked_expected_Delta_y,rowMeans(stacked_expected_Delta_y_boot), pch=19, cex=0.25)
	# lines(c(-1,1),c(-1,1))
	# plot(stacked_expected_Delta_Wy,rowMeans(stacked_expected_Delta_Wy_boot), pch=19, cex=0.25)
	# lines(c(-1,1),c(-1,1))
	
	list(stacked_expected_Delta_y_Averaging_boot=stacked_expected_Delta_y_Averaging_boot,
	     stacked_expected_Delta_Wy_Averaging_boot=stacked_expected_Delta_Wy_Averaging_boot)
	
}