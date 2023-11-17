estimationVectorFieldAdaptive <- function(numGrid, 
                                grid_y_Rel, grid_Wy_Rel, y_Rel, Wy_Rel= Wy_Rel, 
                                numLag, figure=TRUE,alpha=0.5,
                                h.bandwidth=c(1.77*length(y_Rel)^(-1/6),1.77*length(y_Rel)^(-1/6)),
                                weiNorm=rep(1/length(y_Rel),length(y_Rel))){

  
  
  	#++++++++++++++++++++++++++++++++++++++++++++#
	#Estimation of the random vector field (RVF) #
	#++++++++++++++++++++++++++++++++++++++++++++#
	
	#INPUTS	
	#numGrid: number of grid points for the estimate of vector field (numGrid x numGrid is the total number)
  #grid_y_Rel: grid for relative y
	#grid_Wy_Rel: grid for relative y
	#y_Rel: y in relative terms
	#Wy_Rel: W*y_Rel
	#numLag: number of lag to calculate Delta y and Delta Wy 
	#h.bandwidth: width of the bandwidth to be used in the estimation of the kernel
	
	#OUTPUTS
	#estimateDensity_Stack: estimated density for each point in the grid
	#stacked_expected_Delta_y: expected Delta y
	#stacked_expected_Delta_Wy: expected Delta Wy
	
		
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
	
	estimateDensity_Stack <- epaKernelBivAdaptive(evalPoints_Stack=evalPoints_Stack,x=x,invS=invS,detS=detS,ngrid= numGrid, h.bandwidth=h.bandwidth, alpha=alpha, weiNorm=weiNorm)  
	
	#################################################
	#Estimate directions for each point in the grid #
	#################################################
	
	#------------------------------------------------------------#
	# Expected delta of the observed y and Wy from t to t+numLag #
	#------------------------------------------------------------#
	
	delta_y_Rel_matrix <- matrix(delta_y_Rel,nrow=length(delta_y_Rel),ncol=ncol(estimateDensity_Stack))
	
	delta_Wy_Rel_matrix <- matrix(delta_Wy_Rel,nrow=length(delta_y_Rel),ncol=ncol(estimateDensity_Stack))
	
	#stacked_expected_Delta_y <- colSums(delta_y_Rel_matrix*estimateDensity_Stack)/numTrans
	#stacked_expected_Delta_Wy <- colSums(delta_Wy_Rel_matrix*estimateDensity_Stack)/numTrans
	
	#Normalizing weights for local mean estimator (see Bowman and Azzalini, p. 49, 1997)
	normWeights <- colSums(estimateDensity_Stack)

	#Check for zero weights (0/1=0, instead of 0/0 that is indeterminate but it shoudl be zero)
	normWeights[normWeights==0] <-1
	stacked_expected_Delta_y <- colSums(delta_y_Rel_matrix*estimateDensity_Stack)/normWeights
	stacked_expected_Delta_Wy <- colSums(delta_Wy_Rel_matrix*estimateDensity_Stack)/normWeights

	
	if (figure==TRUE){
			plot(y_Rel[,1],Wy_Rel[,1],cex=0.25,pch=19,xlab="y",
		     ylab="Wy",xlim=c(grid_y_Rel[1],grid_y_Rel[numGrid]),
		     ylim=c(grid_Wy_Rel[1],grid_Wy_Rel[numGrid]))
		points(y_Rel[,ncol(y_Rel)],Wy_Rel[,ncol(y_Rel)],cex=0.25,pch=19,col="blue")
		
		for (j in 1: nrow(evalPoints_Stack)){
		  
		  if ( (stacked_expected_Delta_y[j]!=0) |  (stacked_expected_Delta_Wy[j]!=0) ){
		  arrows(evalPoints_Stack[j,1], evalPoints_Stack[j,2],
		         evalPoints_Stack[j,1]+ stacked_expected_Delta_y[j]*(1/numLag),evalPoints_Stack[j,2] + stacked_expected_Delta_Wy[j]*(1/numLag),code=2,length=0.1,lwd=1.2, col="red")
		  }
		}
		lines(c(0,2),c(0,2))
		abline(v=1)
		abline(h=1)
		legend(0,2,c("First year","Last year"),col=c("black","blue"),pch=19,bg="lightgray")		
	}

	
	
	list(estimateDensity_Stack = estimateDensity_Stack,evalPoints_Stack=evalPoints_Stack, stacked_expected_Delta_y=stacked_expected_Delta_y, stacked_expected_Delta_Wy=stacked_expected_Delta_Wy , y_Rel=y_Rel, Wy_Rel=Wy_Rel, grid_y_Rel=grid_y_Rel, grid_Wy_Rel= grid_Wy_Rel, numGrid=numGrid, numLag= numLag, actObs_y_Rel= actObs_y_Rel, futObs_y_Rel= futObs_y_Rel, actObs_Wy_Rel= actObs_Wy_Rel, futObs_Wy_Rel= futObs_Wy_Rel, delta_y_Rel= delta_y_Rel, delta_Wy_Rel= delta_Wy_Rel, weiNorm=weiNorm)
	
	
	
}