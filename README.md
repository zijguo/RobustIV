# RobustIV
This package provides functions for robust inference in the presence of potentially invalid instrumental variables, including the Two-Stage Hard Thresholding(TSHT) method, the Endogeneity Testing method, Searching-Sampling method, the Control Function method, and the SpotIV method.


# Installation
The package can be installed from Github using the following code:
```R
devtools::install_github("https://github.com/zijguo/RobustIV")
```

# Low-dimensional Examples

We use Mroz (1987) data to estimate the effect of education on log earnings. Here, the outcome Y is log earnings (lwage), the exposure D is years of schooling (educ), the candidates of instrument Z are father’s education (fatheduc), mother’s education (motheduc), husband’s education (huseduc), actual labor market experience (exper), its square (expersq), and the exogenous covariate X is age (age). We assume Y is a linear model of D,Z, and X, and D is a linear model of Z and X. 

## TSHT

```R
> data(mroz)
> Y <- mroz[,"lwage"]
> D <- mroz[,"educ"]
> Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc","exper","expersq")])
> X <- mroz[,"age"]
> TSHT.model <- TSHT(Y=Y,D=D,Z=Z,X=X)
> summary(TSHT.model)

Valid Instruments: motheduc fatheduc huseduc 

Relevant Instruments: motheduc fatheduc huseduc 
 
Thus, Majority rule holds 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

BetaHat: 0.08029083 

Confidence Interval of BetaHat: [0.03744511,0.1231365]

```
### Reference
Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), [Confidence intervals for causal effects with invalid instruments by using two-stage hard thresholding with voting](https://doi.org/10.1111/rssb.12275), J. R. Stat. Soc. B, 80: 793-815. 

## Searching-Sampling

```R
> Y <- mroz[,"lwage"]
> D <- mroz[,"educ"]
> Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc","exper","expersq")])
> X <- mroz[,"age"]
> Searching.model <- SearchingSampling(Y,D,Z,X, Sampling = FALSE)
> summary(Searching.model)

Initial set of Valid Instruments: motheduc fatheduc huseduc 

Relevant Instruments: motheduc fatheduc huseduc 
 
Thus, Majority rule holds 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Confidence Interval of BetaHat: [-0.2831715,0.2442571]

> SS.model <- SearchingSampling(Y,D,Z,X)
> summary(SS.model)

Initial set of Valid Instruments: motheduc fatheduc huseduc 

Relevant Instruments: motheduc fatheduc huseduc 
 
Thus, Majority rule holds 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Confidence Interval of BetaHat: [-0.09857149,0.1651428]
```
### Reference
Guo, Z. (2021), [Causal  Inference  with  Invalid  Instruments: Post-selection Problems and A Solution Using Searching and Sampling](https://arxiv.org/abs/2104.06911), Preprint arXiv:2104.06911.


## Control function method
We use the control function method in additive model with continuous outcome and valid IVs.
```R
> Y <- mroz[,"lwage"]
> D <- mroz[,"educ"]
> Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc")])
> X <- as.matrix(mroz[,c("exper","expersq","age")])
> cf.model <- cf(Y~D+I(D^2)+X|Z+I(Z^2)+X,d1 = c(median(D)+1,(median(D)+1)^2),d2 = c(median(D),median(D)^2))
> summary(cf.model)
 _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Coefficients of Control Function Estimators:

              Estimate    Std.Err t value Pr(>|t|)    
(Intercept)  1.2573907  0.7871438   1.597 0.055457 .  
D           -0.1434395  0.1102058   1.302 0.096884 .  
I(D^2)       0.0086426  0.0041004   2.108 0.017817 *  
Xexper       0.0438690  0.0131574   3.334 0.000465 ***
Xexpersq    -0.0008713  0.0003984   2.187 0.014631 *  
Xage        -0.0011636  0.0048634   0.239 0.405511    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Causal effect from D = 12 to 13 :  0.07262563 

> pretest.model <- pretest(Y~D+I(D^2)+X|Z+I(Z^2)+X)
> summary(pretest.model)
 
 _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Level 0.05 Pretest estimator is Control function estimator. 

Hausman Statistic :  1.313563 

P value =  0.2517505 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

Coefficients of Control Function Estimators:

              Estimate    Std.Err t value Pr(>|t|)    
(Intercept)  1.2573907  0.7871438   1.597 0.055457 .  
D           -0.1434395  0.1102058   1.302 0.096884 .  
I(D^2)       0.0086426  0.0041004   2.108 0.017817 *  
Xexper       0.0438690  0.0131574   3.334 0.000465 ***
Xexpersq    -0.0008713  0.0003984   2.187 0.014631 *  
Xage        -0.0011636  0.0048634   0.239 0.405511    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### References
Guo, Z. and D. S. Small (2016), [Control function instrumental variable estimation of nonlinear
causal effect models](https://www.jmlr.org/papers/volume17/14-379/14-379.pdf), The Journal of Machine Learning Research 17(1), 3448–3482.


## SpotIV and ProbitControl

```R
> Y <- mroz[,"lwage"]
> D <- mroz[,"educ"]
> Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc","exper","expersq")])
> X <- mroz[,"age"]
> Y0 <- as.numeric((Y>median(Y)))
> d2 = median(D); d1 = d2+1;
> w0 = apply(cbind(Z,X)[which(D == d2),], 2, mean)
> SpotIV.model <- SpotIV(Y0,D,Z[,-5],X,d1 = d1,d2 = d2,w0 = w0[-5])
> summary(SpotIV.model)

Valid Instruments: motheduc fatheduc huseduc 

Relevant Instruments: motheduc fatheduc huseduc 
 
Thus, Majority rule holds 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

BetaHat: 0.8791451 

CATEHat: 0.1298192 

Confidence Interval of CATEHat: [-0.6774876,0.9371261]
```

```R
> Probit.model <- ProbitControl(Y0,D,Z,X,d1 = d1,d2 = d2,w0 = w0)
> summary(Probit.model)

Valid Instruments: motheduc fatheduc huseduc 

Relevant Instruments: motheduc fatheduc huseduc 
 
Thus, Majority rule holds 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

BetaHat: 0.2118909 

CATEHat: 0.08435024 

Confidence Interval of CATEHat: [0.01999327,0.1487072]
```

### Reference
Li, S., Guo, Z. (2020), [Causal Inference for Nonlinear Outcome Models with Possibly Invalid Instrumental Variables](https://arxiv.org/abs/2010.09922), Preprint arXiv:2010.09922.

# High-dimensional examples with simulation
In this section, we consider the following linear models.

Y = D\beta + Z\alpha + X\phi +u 

D = Z\gamma + X\psi + v

## TSHT

```R
> n = 500; L = 600; s = 3; k = 10; px = 10;
> alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,k),rep(0,L-k))
> phi<-(1/px)*seq(1,px)+0.5; psi<-(1/px)*seq(1,px)+1
> epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
> Z = matrix(rnorm(n*L),n,L)
> X = matrix(rnorm(n*px),n,px)
> epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
> D =  0.5 + Z %*% gamma + X %*% psi + epsilon[,1]
> Y = -0.5 + Z %*% alpha + D * beta + X %*% phi + epsilon[,2]
> summary(TSHT(Y,D,Z,X, method = "DeLasso"))

Valid Instruments: 4 5 6 7 9 10 

Relevant Instruments: 1 2 3 4 5 6 7 8 9 10 
 
Thus, Majority rule holds 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

BetaHat: 1.018932 

Confidence Interval of BetaHat: [0.9841815,1.053682]

```

## Endogeneity test in high dimension
It uses same reduced form estimator as TSHT in each setting.


```R
> n = 500; L = 600; s = 3; k = 10; px = 10;
> alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,k),rep(0,L-k))
> phi<-(1/px)*seq(1,px)+0.5; psi<-(1/px)*seq(1,px)+1
> epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
> Z = matrix(rnorm(n*L),n,L)
> X = matrix(rnorm(n*px),n,px)
> epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
> D =  0.5 + Z %*% gamma + X %*% psi + epsilon[,1]
> Y = -0.5 + Z %*% alpha + D * beta + X %*% phi + epsilon[,2]
> endo.test.model <- endo.test(Y,D,Z,X)
> summary(endo.test.model)
 
Valid Instruments: 4 5 6 7 8 9 10 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
Test statistics Q =  13.71776 
'H0 : Sigma12 = 0' is rejected at level of significance 0.05 
P-value =  0 
Estimated covariance: 0.8922838 

```

### Reference
Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), [Testing endogeneity with high dimensional covariates](https://www.sciencedirect.com/science/article/pii/S0304407618301325), Journal of Econometrics, Elsevier, vol. 207(1), pages 175-187.
