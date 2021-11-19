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
### generate low dimensional data ###
n = 500; L = 10; s = 3
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)
epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
D = 0.5 + Z %*% gamma + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]

### basic usage ###
RobustIV::TSHT(Y,D,Z)
RobustIV::TSHT(Y,D,Z,voting = 'Conservative')
RobustIV::TSHT(Y,D,Z,boot.SHat = TRUE,voting = 'Conservative')

```

### High-dimensional setting
When dimension is high, we use Debiased lasso estimator as reduced-form estimator.
```R
### generate high dimensional data ###
n = 500; L = 600; s = 3; nRelevant = 10
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,nRelevant),rep(0,L-nRelevant))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)
epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
D =  0.5 + Z %*% gamma + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]

### basic usage with debiased lasso ###
RobustIV::TSHT(Y,D,Z,method="DeLasso")

### usage with debiased lasso and voting option ###
RobustIV::TSHT(Y,D,Z,method="DeLasso",voting = 'MP')
RobustIV::TSHT(Y,D,Z,method="DeLasso",voting = 'Conservative')

### usage with debiased lasso and bootstrap threshold option ###
RobustIV::TSHT(Y,D,Z,method="DeLasso",boot.SHat = TRUE)
```

### Reference
Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), [Confidence intervals for causal effects with invalid instruments by using two-stage hard thresholding with voting](https://doi.org/10.1111/rssb.12275), J. R. Stat. Soc. B, 80: 793-815. 

## Searching-Sampling method
We use Searching-Sampling method to resolve post-selection problem.
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
n = 1000; IV.str=0.5; VIO.str=0.4; pi.value<-IV.str*VIO.str
beta = 1;
px <- 10; L = 10; p=L+px
phi<-rep(0,px); psi<-rep(0,px); phi[1:px]<-(1/px)*seq(1,px)+0.5; psi[1:px]<-(1/px)*seq(1,px)+1
rho=0.5; Cov<-(A1gen(rho,p))

### Generate a setting ###
L = 10; s1 = 2; s2 = 2; s=s1+s2
alpha = c(rep(0,L-s1-s2),rep(pi.value,s1),-seq(1,s2)/2); gamma=rep(IV.str,L)

W<-MASS::mvrnorm(n, rep(0, p), Cov)
Z=W[,1:L]
X=W[,(L+1):p]
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
D = 0.5 + Z %*% gamma+ X%*% psi + epsilon[,1]
Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon[,2]
### Using Searching method ###
RobustIV::Searching.Sampling(Y,D,Z,X,Sampling=FALSE)
### Using Searching method with voting option ###
RobustIV::Searching.Sampling(Y,D,Z,X,Sampling=FALSE,voting = 'Conservative')
### Using Sampling method and specifying sampling threshold ### 
RobustIV::Searching.Sampling(Y,D,Z,X,alpha0 = 0.01,boot.SHat = TRUE,voting = 'Conservative')
```

### Reference
Guo, Z. (2021), [Causal  Inference  with  Invalid  Instruments: Post-selection Problems and A Solution Using Searching and Sampling](https://arxiv.org/abs/2104.06911), Preprint arXiv:2104.06911.

## Endogeneity test in high dimension
It uses same reduced form estimator as TSHT in each setting.

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


### generate high-dimensional data ### 
n = 200
rho1 = 0.3
Const = 0.5
s = 10; nRelevant = 7
pz <- 100; px <- 150; p <- pz+px
Lambda <- ar1_cor(p,0.5)
S <- seq(1,nRelevant)
I <- seq(1,pz) # candidates of instrument variables
beta <- 1
phi <- c(seq(0.6,1.5,length.out=s),rep(0,px-s))
psi <- c(seq(1.1,2.0,length.out=s),rep(0,px-s))
sig1 <- 1.5; sig2 <- 1.5 
epsilonSigma = matrix(c(sig1,0.5,0.5,sig2),2,2) 
gamma = Const*c(rep(1,nRelevant-1),rho1,rep(0,pz-nRelevant))
W <- MASS::mvrnorm(n,mu=rep(0,p),Sigma = Lambda)
Z <- W[,1:pz]; X <- W[,(pz+1):p]
epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
D <- 1 + Z%*%gamma + X%*%psi + epsilon[,1]
Y <- -1 + D*beta + X%*%phi + epsilon[,2]

### usage with debaised lasso ###
RobustIV::endo.test(Y,D,Z,X,method = "DeLasso",voting = 'Conservative')

```

### Low-dimensional setting
If you need, you can also implement the method in low dimensional setting.
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
epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
D <- 1 + Z%*%gamma + X%*%psi + epsilon[,1]
Y <- -1 + D*beta + X%*%phi + epsilon[,2]

### basic usage ###
RobustIV::endo.test(Y,D,Z,X,method ="OLS",voting = 'MP')
```

### Reference
Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), [Testing endogeneity with high dimensional covariates](https://www.sciencedirect.com/science/article/pii/S0304407618301325), Journal of Econometrics, Elsevier, vol. 207(1), pages 175-187.

## Control function method
We use the control function method in additive model with continuous outcome and valid IVs.
```R
### control function method ###

### Generate a setting ###
n <- 10000
mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
X <- rnorm(n); Z <- rnorm(n)
err <- MASS::mvrnorm(n,mu=mu,Sigma = V)
u1 <- err[,1]; v2 <- err[,2]
D <- 1+X/8+Z/3+Z^2/8+v2
Y <- 1+X+10*D +10*D^2+u1

### Implement the control function method ###
RobustIV::cf(Y~X+D+I(D^2),D~X+Z+I(Z^2))

### pretest using control function method ###
RobustIV::pretest(Y~X+D+I(D^2),D~X+Z+I(Z^2))

```

### References
Guo, Z. and D. S. Small (2016), [Control function instrumental variable estimation of nonlinear
causal effect models](https://www.jmlr.org/papers/volume17/14-379/14-379.pdf), The Journal of Machine Learning Research 17(1), 3448–3482.

## SpotIV
We use the SpotIV method in semiparametric model with valid IV or possibly invalid IV assumption. This method can take a lot of time when the dimensions are a little bigger. 

```R
### Generate a setting ###
m = 500; J = 5; d1=-1; d2=1;
z0 = rep(0,J); x0 = c(0.1,0.2);
Z <- matrix(rnorm(m * J, 0, 1) , ncol = J, nrow = m)
gam <- c(rep(0.8, floor(J / 2)), rep(-0.8, J - floor(J / 2)))
cov.noise<-matrix(c(1,0.25, 0.25, 1),ncol=2)
noise.vec<-MASS::mvrnorm(m, rep(0,2), cov.noise)
v.vec<-noise.vec[,1]
X<-matrix(runif(m*2), ncol=2)
psi0 <- c(0.1,0.3)
D = 0.5+Z %*% gam + X%*%psi0+v.vec
phi0 <- c(0.2,0.1)
beta0 <- 0.25
u.vec<- noise.vec[,2]
Y = (-0.5 + X%*%phi0 + D * beta0 + u.vec>=0)

### Approximate CATE ###

u1.r<-rnorm(20000,0,sd=1)
cace0 <- mean((as.numeric(-0.5+d1 * beta0 + x0 %*% psi0)+ u1.r )>=0) -
  mean((as.numeric(-0.5+d2 * beta0 + x0 %*% psi0) + u1.r)>=0)
cace0 

### Implement SpotIV method assuming all IVs are valid ###
RobustIV::SpotIV(Y=Y, D=D, Z=Z, X=X, d1 = d1, d2 = d2, w0 = c(z0,x0),invalid = FALSE ,parallel=FALSE)
### Implement SpotIV method assuming there is an invalid IV ###
RobustIV::SpotIV(Y=Y, D=D, Z=Z, X=X, d1 = d1, d2 = d2, w0 = c(z0,x0),invalid = TRUE ,parallel=FALSE)
```
## Control function method for probit model
Especially, we can use the method for probit model assumption.
```R
### Generate a setting ###
   n = 500; J = 5; s = 3; d1=-1; d2=1; z0=c(rep(0, J-1),0.1); 
   Z <- matrix(rnorm(n * J, 0, 1) , ncol = J, nrow = n)
   gam <- c(rep(0.8, floor(J / 2)), rep(-0.8, J - floor(J / 2)))
   cov.noise<-matrix(c(1,0.25, 0.25, 1),ncol=2)
   noise.vec<-MASS::mvrnorm(n, rep(0,2), cov.noise)
   v.vec<-noise.vec[,1]
   D = 0.5+Z %*% gam + v.vec
   pi0 <- c(rep(0, s), 0.8, 0.4)
   beta0 <- 0.25
   u.vec<- noise.vec[,2]
   Y = (-0.5 + Z %*% pi0 + D * beta0 + u.vec>=0)

### Implement SpotIV method when we assume probit model with valid IV assumption ###
RobustIV::ProbitControl(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, w0 = c(z0,x0), method='valid', intercept=TRUE)

### Implement SpotIV method when we assume probit model with possibly invalid IV assumption and majority rule ###
RobustIV::ProbitControl(Y=Y, D=D, Z=Z, X=X, bs.Niter = 40, d1 = d1, d2 = d2, w0 = c(z0,x0), method='majority', intercept=TRUE)

```

### Reference
Li, S., Guo, Z. (2020), [Causal Inference for Nonlinear Outcome Models with Possibly Invalid Instrumental Variables](https://arxiv.org/abs/2010.09922), Preprint arXiv:2010.09922.
