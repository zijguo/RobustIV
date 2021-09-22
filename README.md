# RobustIV
This package provides functions for robust inference in the presence of potentially invalid instrumental variables, including the Two-Stage Hard Thresholding(TSHT) method, the Endogeneity Testing method, Searching-Sampling method, the Control Function method, and the SpotIV method.


# Installation
The package can be installed from Github using the following code:
```R
# install.packages("devtools")
library(devtools)
devtools::install_github("https://github.com/zijguo/RobustIV")
```

# Examples

First, we consider the linear model.

## TSHT
### Low-dimensional setting
When dimension is low, we implement TSHT with OLS reduced-form estimator. 
```R
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
```R
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
```R
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
```R
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

## Searching-Sampling method
We use Searching-Sampling method in low-dimensional setting to resolve post-selection problem.
```R
### Define covariance matrix ###
A1gen<-function(rho,p){
    A1=matrix(0,p,p)
    for(i in 1:p){
     for(j in 1:p){
        A1[i,j]<-rho^(abs(i-j))
      }
    }
    A1
  }
n = 1000;
IV.str=0.5;
VIO.str=0.4;
pi.value<-IV.str*VIO.str
beta = 1;
px <- 10;
L = 10;
p=L+px
phi<-rep(0,px)
psi<-rep(0,px)
phi[1:px]<-(1/px)*seq(1,px)+0.5
psi[1:px]<-(1/px)*seq(1,px)+1
rho=0.5
Cov<-(A1gen(rho,p))

### Generate a setting ###
L = 10;
s1 = 2;
s2 = 2;
s=s1+s2
alpha = c(rep(0,L-s1-s2),rep(pi.value,s1),-seq(1,s2)/2);
gamma=rep(IV.str,L)

epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
W<-mvrnorm(n, rep(0, p), Cov)
Z=W[,1:L]
X=W[,(L+1):p]
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = 0.5 + Z %*% gamma+ X%*% psi + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon[,2]
### Using Searching method ###
Searching.Sampling(Y,D,Z,X,Sampling=FALSE)
### Using Searching method with maximum clique option ###
Searching.Sampling(Y,D,Z,X,Sampling=FALSE,max_clique=TRUE)
### Using Sampling method with maximum clique and specifying sampling threshold ### 
Searching.Sampling(Y,D,Z,X,max_clique=TRUE,alpha0 = 0.01)
```

## Control function method
We use the control function method in additive model.
```R
### control function method ###

### Generate a setting ###
n <- 10000
mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
X <- rnorm(n); Z <- rnorm(n)
err <- mvrnorm(n,mu=mu,Sigma = V)
u1 <- err[,1]; v2 <- err[,2]
D <- 1+X/8+Z/3+Z^2/8+v2
Y <- 1+X+10*D +10*D^2+u1

### Implement the control function method ###
cf(Y~X+D+I(D^2),D~X+Z+I(Z^2))

### pretest using control function method ###
pretest(Y~X+D+I(D^2),D~X+Z+I(Z^2))
```

## SpotIV method
We use the SpotIV method in semiparametric model.
```R
### Generate a setting ###
n = 500; J = 5; s = 3; d1=-1; d2=1; z0=c(rep(0, J-1),0.1); x0 = c(0.1,0.2)
Z <- matrix(rnorm(n * J, 0, 1) , ncol = J, nrow = n)
gam <- c(rep(0.8, floor(J / 2)), rep(-0.8, J - floor(J / 2)))
cov.noise<-matrix(c(1,0.25, 0.25, 1),ncol=2)
noise.vec<-mvrnorm(n, rep(0,2), cov.noise)
v.vec<-noise.vec[,1]
X<-matrix(runif(n*2), ncol=2)
D = 0.5+Z %*% gam + v.vec
pi0 <- c(rep(0, s), 0.8, 0.4)
beta0 <- 0.25
u.vec<- noise.vec[,2]
Y = (-0.5 + Z %*% pi0 + D * beta0 + u.vec>=0)

### Implement SpotIV method assuiming all IVs are valid ###
SpotIV(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, V= 1:J, w0 = c(z0,x0), parallel=FALSE)

### Implement SpotIV method without parallel computing option ###
SpotIV(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, w0 = c(z0,x0), parallel=FALSE)

### Implement SpotIV method without parallel computing option ###
SpotIV(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, w0 = c(z0,x0), parallel=TRUE)

### Implement SpotIV method when we assume probit model with valid IV assumption ###
ProbitControl(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, w0 = c(z0,x0), method='valid', intercept=TRUE)

### Implement SpotIV method when we assume probit model with possibly invalid IV assumption and majority rule ###
ProbitControl(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, w0 = c(z0,x0), method='majority', intercept=TRUE)
```

