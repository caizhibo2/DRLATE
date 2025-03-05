library(dcov)
library(glmnet)

HOPG <- function(x, y, xp=x, reducedDim = 5, isScreen = TRUE, screenDim = 0, isFast = TRUE, isPCA = TRUE) 
{   
  
  # data information
  n <- nrow(x)
  p <- ncol(x)
  y = as.matrix(y)
  q = ncol(y)
  
  if (reducedDim >= p) 
    return(list(x = x, xp = xp, dir = diag(1, p), allDir=diag(1, p)))

  # variable screening before reduction
  if (isScreen) 
  {
    if (screenDim > 0)
      sn = min(screenDim, p)
    else 
      sn= min(floor(n/log(n)), p)
    
    if (p > sn)
    {
      dcr = rep(0, p)
      for (i in 1:p)
      {
        for (j in 1:q)
          dcr[i] =  dcr[i] + dcor2d(y[,j], x[,i])
      }
      
      I = order(dcr, decreasing = TRUE)
      Iscreen = I[1:sn]
      x = x[,Iscreen]
      xp = xp[,Iscreen]
    }
  } 
  else {
    Iscreen = 1:p 
  }
  pD = length(Iscreen)
  
  # PCA process eliminating the unimportant directions
  if (isPCA) {
    s <- eigen( t(x)%*%x/n + diag(1/n, pD))    
    value0 = abs(s$values)
    value = value0
    Ieig = which( value > 1.0e-8)
    value[Ieig] = 1/value[Ieig];
    value[-Ieig] = 0;
    ss <- s$vectors %*% diag( value^0.5) %*% t(s$vectors)
    
    x <- x %*% ss
    xp <- xp %*% ss
  } 
  
  # Calculate parameter alpha in kernel function
  V = as.matrix(x)
  vdcor = rep(0, pD)  # parameter alpha - dcor based on V
  for (i in 1:pD)
  {
    for (j in 1:q)
      # vdcor[i] =  vdcor[i] + dcor2d(y[,j]-sum(V[,i]*y[,j])/sum(V[,i]^2)*V[,i], V[,i])
      vdcor[i] =  vdcor[i] + dcor2d(y[,j]-cov(V[,i], y[,j])/var(V[,i])*V[,i], V[,i])
  }
  vdcor = vdcor/max(vdcor) # modify to 0-1
  #  vdcor = vdcor^0.25
  vdcorSum = sum(vdcor) # for bandwidth calculation
  
  # Calculate bandwidth for gradient estimation
  h <- 1.2*mean(apply(V, 2, sd))/n^(1/(vdcorSum+4))

  # Find a small number of points for the calculation of OPG
  if (isFast) {
    I <- Ifast(V)
  } else {
    I <- 1:n
  }

  # compute sum of the OPG matrix
  M0 <- 0 # sum of OPG matrix
  onex = cbind(x, rep(1, n))  
  for (i in I) {
    dx <- x - matrix(x[i, ], n, pD, byrow = TRUE)
    ker <- AlphaKernel(dx, h, vdcor)
    
    xk <- t(onex*matrix(rep(ker, pD+1), ncol = pD+1))
    abi <- solve(xk %*% onex + diag(log(p)/n * 1e-9, pD+1), xk%*%y) 
    
    M0 <- M0 + abi[1:pD,] %*%  t(abi[1:pD,]) # sum of outer product
  }
  
  # eigenvalue decomposition
  M0.eigen <- eigen(M0)
  Ii <- order(M0.eigen$values,decreasing = T)
  val = M0.eigen$values[Ii]
  B <- M0.eigen$vectors[,Ii]

  # compute outputs
  x = x %*% B[,1:reducedDim] # reduced x
  xp = xp %*% B[,1:reducedDim] # reduced xp
  dirB = ss %*% B # directions after screen 
  dirB = dirB/matrix(sqrt(colSums(dirB^2)), pD, pD, byrow=TRUE) # standardize
  
  allDir = matrix(0, p, p)
  allDir[Iscreen, 1:pD] = dirB # all directions in the order of importance
  dir = allDir[,1:reducedDim]
  
  return(list(x = x, xp = xp, dir = dir, allDir = allDir))
  
}

AlphaKernel <- function(U, h, alpha, type = "Gaussian") {
  # Calculate weights based on kernel functions
  # INPUT:
  # U - n x p matrix
  # h - bandwidth
  # alpha - parameter
  # type - type of kernel function
  # OUTPUT:
  # res - n vector of kernel weights
  
  # data information
  U = as.matrix(U)
  n = nrow(U)
  p = ncol(U)
  
  # data processing for kernel calculation
  S = matrix( h^alpha, n, p, byrow = TRUE)
  h2 = h^sum(alpha)
  censorU = U / S * (-S <= U & U <= S)

  # compute for weights
  if (tolower(type) == "uniform") {
    res <- apply(censorU, 1, prod) / h2 / 2^p
  }
  else if (tolower(type) == "triangular") {
    res <- apply((1 - abs(censorU)), 1, prod) / h2
  }
  else if (tolower(type) == "parabolic") {
    res <- apply((1 - censorU^2), 1, prod) / h2 * 0.75^p
  }
  else if (tolower(type) == "cosine") {
    res <- apply(cos(censorU*pi/2), 1, prod) / h2 * (pi/4)^p
  } 
  else if (tolower(type) == "gaussian") {
    res <- exp( - apply( (U/S)^2, 1, sum ) / 2 ) / h2 / (2*pi)^(p/2)
  }
  else {
    stop("Kernel function not included.")
  }
  
  return(res)
}

Ifast <- function(V)
{
  # Select points for faster computation
  # Input: 
  # V - data matrix for selection
  # Output:
  # Ifast - indices for selected points
  
  n = nrow(V)
  pD = ncol(V)
  n0 <- min(100, n)
  npoints <- n0 + floor( (n-n0)^0.7 )
  
  Ifast = 1:n
  if (ceiling(n/npoints) > 1)
  {
    k = floor(n/npoints);
    q = n - npoints*k
    K = rep(k, npoints)
    if (q > 0)  K[1:q] = K[1:q]+1
    
    Vremain = V
    Iremain = 1:n;
    Ifast = c();
    for (jk in 1:npoints)
    {
      mi = length(Iremain);
      if (mi > 1)
      {
        xi = Vremain[1,];    
        di = 0;
        for (ip in 1:pD) di = di + (Vremain[,ip] - xi[ip])^2
        
        Ii = order(di);    
        Ii = Ii[1:min(K[jk], mi)]; # select k nearest points to xi
        
        if (length(Ii)>1)
        {
          vi = colMeans(Vremain[Ii,]) 
          di = 0 
          for (ip in 1:pD) di = di + (Vremain[Ii,ip] - vi[ip])^2
          
          J = order(di)
          Ifast = c(Ifast, Iremain[Ii[J[1]]]); # select points with min variance
        }else
        {
          Ifast = c(Ifast, Iremain[Ii])
        }
        
        Iremain = setdiff(Iremain, Iremain[Ii]);
        Vremain = V[Iremain,]
      }
      else
      {
        Ifast = c(Ifast, Iremain[1])
      }
    }
  }        
  
  return(Ifast)
  
}



