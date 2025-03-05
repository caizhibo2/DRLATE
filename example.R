source("DRLATE.R")

### Example 1 ########################

X <- matrix(runif(200*100), 200, 100)
T <- (runif(200) > 0.5) + 0

y = 2*T + 10*rnorm(200) +
  10*sin(pi*X[ ,1]*X[ ,2]) + 20*(X[ ,3]-0.5)^2 + 10*X[ ,4] + 5*X[ ,5]

DRLATE(X, y, T, 0.05, ndim = 10, method = "HOPG")


### Example 2 ########################

X <- matrix(runif(200*100), 200, 100)
T <- (runif(200) > 0.5) + 0

B = matrix(0, 100, 20)
B[sample(1:100, 20), 1:20] = matrix(rnorm(20*20), ncol=20)
B = svd(B)$u
Z = X %*% B
y = 2*T + 10*rnorm(200) + 
  10*sin(pi*Z[ ,1]*Z[ ,2]) + 20*(Z[ ,3]-0.5)^2 + 10*Z[ ,4] + 5*Z[ ,5]

DRLATE(X, y, T, 0.05, ndim = 10, method = "HOPG")
