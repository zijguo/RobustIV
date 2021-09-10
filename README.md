# RobustIV
This package provides functions for robust inference in the presence of potentially invalid instrumental variables, including the Two-Stage Hard Thresholding(TSHT) method, the Endogeneity Testing method, Searching-Sampling method, the Control Function method, and the SpotIV method.


# Installation
The package can be installed from Github using the following code:
```
# install.packages("devtools")
library(devtools)
devtools::install_github("https://github.com/zijguo/RobustIV")
```

# Examples

First, we consider the linear model.

## TSHT
### Low-dimensional setting.
When dimension $p$ is low, we implement TSHT with OLS reduced-form estimator. 
```
library(RobustIV)
library(MASS)
n = 500; L = 10; s = 3
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = 0.5 + Z %*% gamma + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]
TSHT(Y,D,Z)
TSHT(Y,D,Z,max_clique=TRUE)
```


