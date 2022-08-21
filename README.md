# RobustIV
This package provides functions for Inference for the treatment effect with possibly invalid instrumental variables, including the Two-Stage Hard Thresholding(TSHT) method, the Endogeneity Testing method, and Searching-Sampling method.


# Installation
The package can be installed from Github using the following code:
```R
devtools::install_github("https://github.com/zijguo/RobustIV")
```
Before using the package, we can use the following code:
```R
library(RobustIV)
```

# Low-dimensional Examples

We use pseudodata provided by Youjin Lee, which is generated mimicing the structure of Framingham Heart Study data. We assume Y is a linear model of D,Z, and X, and D is a linear model of Z and X. 

## TSHT

```R
> data("lineardata")
> Y <- lineardata[,"Y"]
> D <- lineardata[,"D"]
> Z <- as.matrix(lineardata[,c("Z.1","Z.2","Z.3","Z.4","Z.5","Z.6","Z.7","Z.8")])
> X <- as.matrix(lineardata[,c("age","sex")])
> TSHT.model <- TSHT(Y=Y,D=D,Z=Z,X=X)
> summary(TSHT.model)

Relevant IVs: Z.3 Z.4 Z.5 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
 betaHat    Std.Error  CI(2.5%) CI(97.5%)  Valid IVs  
 0.05166598 0.02015546 0.012162 0.09116995 Z.3 Z.4 Z.5

```
### Reference
Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), [Confidence intervals for causal effects with invalid instruments by using two-stage hard thresholding with voting](https://doi.org/10.1111/rssb.12275), J. R. Stat. Soc. B, 80: 793-815. 

## Searching-Sampling

```R
> Y <- lineardata[,"Y"]
> D <- lineardata[,"D"]
> Z <- as.matrix(lineardata[,c("Z.1","Z.2","Z.3","Z.4","Z.5","Z.6","Z.7","Z.8")])
> X <- as.matrix(lineardata[,c("age","sex")])
> Searching.model <- SearchingSampling(Y,D,Z,X, Sampling = FALSE)
> summary(Searching.model)

Initial set of Valid Instruments: Z.3 Z.4 Z.5 

Plurality rule holds.
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Confidence Interval for Beta: [-0.0356797,0.1422332]
> Sampling.model <- SearchingSampling(Y,D,Z,X)
> summary(Sampling.model)

Initial set of Valid Instruments: Z.3 Z.4 Z.5 

Plurality rule holds.
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Confidence Interval for Beta: [-0.02297164,0.1295251]
```
### Reference
Guo, Z. (2021), [Causal  Inference  with  Invalid  Instruments: Post-selection Problems and A Solution Using Searching and Sampling](https://arxiv.org/abs/2104.06911), Preprint arXiv:2104.06911.

# High-dimensional examples (Simulated data)
In this section, we consider the following linear models.

Y = D\beta + Z\alpha + X\phi +u 

D = Z\gamma + X\psi + v

## TSHT

```R
> set.seed(1)
> n = 500; L = 600; s = 3; k = 10; px = 10;
> alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,k),rep(0,L-k))
> phi<-(1/px)*seq(1,px)+0.5; psi<-(1/px)*seq(1,px)+1
> epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
> Z = matrix(rnorm(n*L),n,L)
> X = matrix(rnorm(n*px),n,px)
> epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
> D =  0.5 + Z %*% gamma + X %*% psi + epsilon[,1]
> Y = -0.5 + Z %*% alpha + D * beta + X %*% phi + epsilon[,2]
> TSHT.model <- TSHT(Y,D,Z,X,method = "Fast.DeLasso")
> summary(TSHT.model)

Relevant IVs: 1 2 3 4 5 6 7 8 9 10 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
 betaHat   Std.Error CI(2.5%)  CI(97.5%) Valid IVs     
 0.9880546 0.0157771 0.9571321 1.018977  4 5 6 7 8 9 10

```

## Endogeneity test in high dimension
It uses same reduced form estimator as TSHT in each setting.


```R
> set.seed(1)
> n = 500; L = 600; s = 3; k = 10; px = 10;
> alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,k),rep(0,L-k))
> phi<-(1/px)*seq(1,px)+0.5; psi<-(1/px)*seq(1,px)+1
> epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
> Z = matrix(rnorm(n*L),n,L)
> X = matrix(rnorm(n*px),n,px)
> epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
> D =  0.5 + Z %*% gamma + X %*% psi + epsilon[,1]
> Y = -0.5 + Z %*% alpha + D * beta + X %*% phi + epsilon[,2]
> endo.test.model <- endo.test(Y,D,Z,X, invalid = TRUE)
> summary(endo.test.model)

Valid Instruments: 4 5 6 7 8 9 10 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
Estimated covariance: 0.7690532 
Test statistics Q =  13.47976 
P-value =  0 
'H0 : Sigma12 = 0' is rejected at the significance level 0.05 .
```

### Reference
Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), [Testing endogeneity with high dimensional covariates](https://www.sciencedirect.com/science/article/pii/S0304407618301325), Journal of Econometrics, Elsevier, vol. 207(1), pages 175-187.
