library(dcov)
library(Rdimtools)
source("HOPG.R")

DRLATE <- function(X, y, T, alpha, ndim, method = "HOPG", screen = TRUE) {
  # DRLATE (Dimension-Reduced Linearly Adjusted Treatment Effect Estimator) method
  # Input: 
  # X - covariates matrix
  # y - response
  # T - treatment
  # alpha - significance level
  # ndim - number of reduced dimension
  # method - method for dimension reduction ("HOPG", "SIR", "PCA")
  # screen - whether screening is used
  # Output:
  # a list of ATE, SE and CI
  
  n <- nrow(X); 
  X <- as.matrix(X)
  
  ####### screening #######
  if (screen & floor(n/log(n)) < p & floor(n/log(n)) > ndim) {
    cor <- c()
    for (i in 1:ncol(sample_X)) {
      if (sd(X[,i]) == 0) {
        cbind(cor, c(i, 0, 0))
      } else {
        cor <- cbind(cor, c(i, cor(y, X[,i]), dcor2d(y, X[,i])))
      }
    }
    # order.cor <- cor[,order(abs(cor[2,]), decreasing = TRUE)]
    order.dcor <- cor[,order(abs(cor[3,]), decreasing = TRUE)]
    # sample_X_cor <- sample_X[, order.cor]
    X <- X[, order.dcor[1:floor(n/log(n))]]
  }
  
  ####### Dimension Reduction #######
  if (toupper(method) == "HOPG") {
    print("Doing DRLATE with HOPG.")
    HOPG <- HOPG(x = X, y = y, reducedDim = ndim, isScreen = FALSE)
    B <- HOPG$dir
  }
  else if (toupper(method) == "PCA") {
    print("Doing DRLATE with PCA.")
    PCA <- do.pca(X = X, ndim = ndim)
    B <- PCA$projection
  }
  else if (toupper(method) == "SIR") {
    print("Doing DRLATE with SIR.")
    SIR <- do.sir(X = X, response = y, ndim = ndim)
    B <- SIR$projection
  }
  else {
    stop("Dimension reduction method not included.")
  }
  
  ######## linear adjusted ATE ########
  Z <- X %*% B
  centeredZ <- Z - matrix(apply(Z, 2, mean), nrow(Z), ncol(Z), byrow = TRUE)
  
  data <- data.frame(cbind(Z, centeredZ, T, y))
  LinearModel <- lm(y~T+Z+T*centeredZ, data)
  summaryLinearModel <- summary(LinearModel)
  
  ATE <- as.numeric(LinearModel$coefficients["T"])
  ATE.SE <- summaryLinearModel$coefficients["T", "Std. Error"]
  CI = ATE + ATE.SE*qt(c(alpha/2, 1-alpha/2), df = LinearModel$df)
  
  return(list(ATE = ATE, ATE.SE = ATE.SE, CI = CI))
  
}