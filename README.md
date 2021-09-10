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
### Low-dimensional setting
When dimension is low, we implement TSHT with OLS reduced-form estimator. 
```
library(RobustIV)
library(MASS)
### generate low dimensional data ###
n = 500; L = 10; s = 3
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = 0.5 + Z %*% gamma + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]

### basic usage ###
TSHT(Y,D,Z)

### usage with maximum clique
TSHT(Y,D,Z,max_clique=TRUE)
```

### High-dimensional setting
When dimension is high, we use Debiased lasso estimator as reduced-form estimator.
```
### generate high dimensional data ###
n = 500; L = 600; s = 3; nRelevant = 10
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,nRelevant),rep(0,L-nRelevant))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D =  0.5 + Z %*% gamma + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]

### basic usage with debiased lasso ###
TSHT(Y,D,Z,method="DeLasso")

### usage with debiased lasso and maximum clique option ###
TSHT(Y,D,Z,method="DeLasso",max_clique=TRUE)
```

## Endogeneity test
It uses same reduced form estimator as TSHT in each setting.

### Low-dimensional setting
```
### Generate low-dimensional data ###
n <- 1000; pz <- 9; px <- 5;
p <- pz+px
s = 10; nRelevant = 7
beta <- 1
phi <- seq(0.6,1.0,length.out=px)
psi <- seq(1.1,1.5,length.out=px)
gamma = c(rep(1,nRelevant),rep(0,pz-nRelevant))
epsilonSigma = matrix(c(1.5,0.5,0.5,1.5),2,2)
W <- matrix(rnorm(n*p),n,p)
Z <- W[,1:pz]; X <- W[,(pz+1):p]
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D <- 1 + Z%*%gamma + X%*%psi + epsilon[,1]
Y <- -1 + D*beta + X%*%phi + epsilon[,2]

### basic usage ###
endo.test(Y,D,Z,X)
```

### High-dimensional setting
```
### Define covariance structure of Z and X ###
ar1_cor <- function(p, rho) {
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}
Lambda <- ar1_cor(p,0.5)

### generate high-dimensional data ### 
n = 200
rho1 = 0.3
Const = 0.5
s = 10; nRelevant = 7
pz <- 100; px <- 150; p <- pz+px
S <- seq(1,nRelevant)
I <- seq(1,pz) # candidates of instrument variables
beta <- 1
phi <- c(seq(0.6,1.5,length.out=s),rep(0,px-s))
psi <- c(seq(1.1,2.0,length.out=s),rep(0,px-s))
sig1 <- 1.5; sig2 <- 1.5 
epsilonSigma = matrix(c(sig1,0.5,0.5,sig2),2,2) 
gamma = Const*c(rep(1,nRelevant-1),rho1,rep(0,pz-nRelevant))
W <- mvrnorm(n,mu=rep(0,p),Sigma = Lambda)
Z <- W[,1:pz]; X <- W[,(pz+1):p]
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D <- 1 + Z%*%gamma + X%*%psi + epsilon[,1]
Y <- -1 + D*beta + X%*%phi + epsilon[,2]

### usage with debaised lasso ###
endo.test(Y,D,Z,X,method = "DeLasso")
```
