epaKernelBivAdaptive <- function(evalPoints_Stack,x,
                                 invS=solve(cov(x)),
                                 detS=det(cov(x)),
                                 weiNorm=rep(1/nrow(x),nrow(x)),
                                 alpha=0.5, ngrid,
                                 h.bandwidth=c(1.77*nrow(x)^(-1/6),1.77*nrow(x)^(-1/6)),
                                 adaptive=TRUE){
    
    ###################################
    ## Calculation of lambdas from the density pilot
    
    ##Pilot estimate, Silverman (1986), p. 101
    # pb = txtProgressBar(min = 1, max = nrow(x), initial = 0) 
    # densityPilot <- foreach(i=1:nrow(x), .combine=c, .verbose=FALSE) %do%{  
    #     setTxtProgressBar(pb,i)
    #     normKernelBiv(y=x[i,],x=x,invS=invS,detS=detS,weiNorm=weiNorm,lambdas=rep(1,dim(x)[1]), 
    #                   h.bandwidth=h.bandwidth)
    #     # print(paste("i = ",i,"/",nrow(x)))
    # }
    # close(pb)
    
    if (adaptive==TRUE){
        print("Start pilot estimate")
        densityPilot <- normKernelBiv_vec(y=x,x=x,invS=invS,detS=detS,
                                          weiNorm=weiNorm,lambdas=rep(1,dim(x)[1]))
        print("end estimate pilot")
        
        # correct 
        lambdas <- (densityPilot/exp(sum(weiNorm*log(densityPilot))))^(-alpha)
    }else{
        lambdas <- rep(1,dim(x)[1])
    }
    
    print("Start estimate of the stochastic kernel")
    estimateDensity_Stack <- epaKernelBiv_vec(y=evalPoints_Stack,x=x,
                                          invS=invS,detS=detS,
                                          weiNorm=weiNorm,lambdas=lambdas,
                                          h.bandwidth=h.bandwidth)
    print("End estimate of the stochastic kernel")
    return(estimateDensity_Stack) # nobs x (ngrid*ngrid) 
}                    
