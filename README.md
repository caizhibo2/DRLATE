# DRLATE
The implementation for Dimension-Reduced Linearly Adjusted Treatment Effect Estimator (DRLATE) in R. 
DRLATE first applies the High-dimensional Outer-Product-of-Gradient (HOPG) method for sufficient dimension reduction, preserving the relationship between covariates and the outcome. 
It then estimates the ATE using linear regression on the reduced-dimensional data. 
DRLATE is particularly suited for ultrahigh-dimensional data (where the number of features exceeds the sample size) and offers interpretability as a linear model.

## Codes
**DRLATE.R** - the R code for DRLATE.  
**HOPG.R** - the R code for HOPG.  
**example.R** - the small examples .  

## Usage
Two small examples for the usage of DRLATE. 

- Example 1   
``` r
source("DRLATE.R")

X <- matrix(runif(200*100), 200, 100)
T <- (runif(200) > 0.5) + 0

y = 2*T + 10*rnorm(200) +
  10*sin(pi*X[ ,1]*X[ ,2]) + 20*(X[ ,3]-0.5)^2 + 10*X[ ,4] + 5*X[ ,5]

DRLATE(X, y, T, 0.05, ndim = 10, method = "HOPG")
```

- Example 2  
``` r
source("DRLATE.R")

X <- matrix(runif(200*100), 200, 100)
T <- (runif(200) > 0.5) + 0

B = matrix(0, 100, 20)
B[sample(1:100, 20), 1:20] = matrix(rnorm(20*20), ncol=20)
B = svd(B)$u
Z = X %*% B
y = 2*T + 10*rnorm(200) + 
  10*sin(pi*Z[ ,1]*Z[ ,2]) + 20*(Z[ ,3]-0.5)^2 + 10*Z[ ,4] + 5*Z[ ,5]

DRLATE(X, y, T, 0.05, ndim = 10, method = "HOPG")
```
