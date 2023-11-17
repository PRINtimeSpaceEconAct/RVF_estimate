epaKernelBiv_vec <- function(y,x,
		  lambdas=rep(1,dim(x)[1]),
		  invS=solve(cov(x)),
		  detS=det(cov(x)),
		  weiNorm=rep(1/dim(x)[1],dim(x)[1]),
		  h.bandwidth=c(1.77*nrow(x)^(-1/6),1.77*nrow(x)^(-1/6))){
		  
    #y: vector of NEval x 2
    #x: matrix with observation of dimensions NObs x 2
    # weiNorm: weights of observations
    
    
    epaPDF <- function(z){
        (2/pi)*(1-z)
    }
    
    NObs = nrow(x)
    NEval = nrow(y)
    MEval1 = matrix(rep(y[,1],NObs),NEval,NObs)
    MEval2 = matrix(rep(y[,2],NObs),NEval,NObs)
    MObs1 = matrix(rep(x[,1],NEval),NEval,NObs,byrow=TRUE)
    MObs2 = matrix(rep(x[,2],NEval),NEval,NObs,byrow=TRUE)
    MLambdas = matrix(rep(lambdas,NEval),NEval,NObs,byrow=TRUE)
    
    M1 = MEval1 - MObs1
    M2 = MEval2 - MObs2 
    
    # quadratic form 
    QF = invS[1,1]*M1^2 + invS[1,2]*M1*M2 + invS[2,1]*M1*M2 + invS[2,2]*M2^2
    
    # argument of standard gaussian kernel 
    MArg = (1/prod(h.bandwidth))*QF/MLambdas^2
    MK =  epaPDF(MArg)
    MK[MArg >= 1] = 0
    MWei = matrix(rep(weiNorm,NEval),NEval,NObs,byrow=TRUE)
    fx = (detS)^(-1/2)/(prod(h.bandwidth))*(MWei*MK/MLambdas^2)
    
    return(t(fx))   # NObs x NEval
    
}                    


