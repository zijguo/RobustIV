##################################
### Date: 4/25/2022
### Author: Zhenyu WANG

library(SIHR)
library(MASS)
library(intervals)
source("helpers.R")

SearchingSampling <- function(Y, D, Z, X, intercept=TRUE, 
                              lowd=TRUE, 
                              robust=TRUE,
                              CI.init = NULL,
                              a=0.6,
                              Sampling=TRUE, 
                              rho=NULL, M=1000, prop=0.1, filtering=TRUE){
  ###############################################
  ## Arguments:
  ## lowd      If TRUE, use OLS method; else use SIHR LF method.
  ## robust    If TRUE, fit model for heteroscedastic model; else homoscedastic 
  ## CI.init   initial range [L, U], if not specified we provided a default method to construct
  ## a         the grid size is set as n^{-a}
  ## Sampling  If TRUE, use sampling approach; else use searching approach.
  ## rho       For sampling method, the initial rho (corresponding to lambda in paper)
  ## M         For sampling method, sampling times
  if(is.null(X)) W = Z else W = cbind(Z, X)
  n = length(Y); pz = ncol(Z); p = ncol(W)
  
  ## Preparation
  if(intercept) W = cbind(W, 1)
  
  if(lowd){
    covW = t(W)%*%W/n
    U = solve(covW) # precision matrix
    WUMat = (W%*%U)[,1:pz]
    ## OLS estimators
    qrW = qr(W)
    ITT_Y = qr.coef(qrW, Y)[1:pz]
    ITT_D = qr.coef(qrW, D)[1:pz]
    resid_Y = as.vector(qr.resid(qrW, Y))
    resid_D = as.vector(qr.resid(qrW, D))
    ## Computing (co)variance matrices
    if(robust){
      SigmaSqY = sum(resid_Y^2)/(n-1)
      SigmaSqD = sum(resid_D^2)/(n-1)
      SigmaYD = sum(resid_Y * resid_D)/(n-1)
      V.Gamma = SigmaSqY * U[1:pz, 1:pz]
      V.gamma = SigmaSqD * U[1:pz, 1:pz]
      C = SigmaYD * U[1:pz, 1:pz]
    }else{
      V.Gamma = (t(WUMat)%*%diag(resid_Y^2)%*%WUMat)/n
      V.gamma = (t(WUMat)%*%diag(resid_D^2)%*%WUMat)/n
      C = (t(WUMat)%*%diag(resid_Y * resid_D)%*%WUMat)/n
    }
  }else{
    ## LF estimators
    robust=FALSE # we only consider homoscedastic setting for LF estimator
    init_Y = Lasso.init(W, Y) 
    init_D = Lasso.init(W, D)
    ## residual
    resid_Y = as.vector(Y - W%*%init_Y)
    resid_D = as.vector(D - W%*%init_D)
    ## Debias
    ITT_Y = rep(NA, pz); ITT_D = rep(NA, pz)
    U = matrix(NA, nrow=ncol(W), ncol=pz)
    if(intercept) W_no_intercept = W[,-ncol(W)] else W_no_intercept = W_no_intercept
    for(i in 1:pz){
      loading = rep(0, p+as.integer(intercept))
      loading[i] = 1
      ITT_Y[i] = GLM_LF(W_no_intercept, Y, loading, intercept.loading=FALSE, intercept=intercept)$prop.est
      model_D = GLM_LF(W_no_intercept, D, loading=loading, intercept.loading=FALSE, intercept=intercept)
      ITT_D[i] = model_D$prop.est
      if(intercept){
        U[-nrow(U),i] = (model_D$proj)[-1]
        U[nrow(U),i] = (model_D$proj)[1]
      }else{
        U[,i] = model_D$proj
      }
    }
    WUMat = W%*%U
    
    SigmaSqY = sum(resid_Y^2)/(n-1)
    SigmaSqD = sum(resid_D^2)/(n-1)
    SigmaYD = sum(resid_Y * resid_D)/(n-1)
    
    Temp = t(WUMat)%*%WUMat / n
    V.Gamma = SigmaSqY * Temp
    V.gamma = SigmaSqD * Temp
    C = SigmaYD * Temp
  }
  
  TSHT.out <- TSHT.Init(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C)
  V0.hat = sort(TSHT.out$VHat)
  ## Construct range [L, U]
  if(is.vector(CI.init)){
    CI.init.union = t(as.matrix(sort(CI.init)))
  }else{
    ## current method to select initial [L, U]
    var.beta = 1/n * (diag(V.Gamma)/ITT_D^2 + diag(V.gamma)*ITT_Y^2/ITT_D^4 - 2*diag(C)*ITT_Y/ITT_D^3)
    var.beta = var.beta[V0.hat]
    CI.init = matrix(NA, nrow=length(V0.hat), ncol=2)
    CI.init[,1] = (ITT_Y/ITT_D)[V0.hat] - sqrt(log(n)*var.beta)
    CI.init[,2] = (ITT_Y/ITT_D)[V0.hat] + sqrt(log(n)*var.beta)
    uni = Intervals(CI.init)
    CI.init.union = as.matrix(interval_union(uni))
  }
  
  # Construct beta.grid
  beta.grid = grid.CI(CI.init.union, grid.size=n^{-a})
  
  if(Sampling){
    ## Sampling Method
    CI.sampling = Searching.CI.sampling(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet=V0.hat,
                                        beta.grid = beta.grid, rho=rho, M=M, prop=prop, filtering=filtering)
    CI=CI.sampling$CI
    rule=CI.sampling$rule
  }else{
    ## Searching Method
    CI.searching = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = V0.hat,
                                beta.grid = beta.grid)
    CI=CI.searching$CI
    rule=CI.searching$rule
  }
  returnList <- list(CI=CI, check=rule, VHat=V0.hat, SHat=TSHT.out$SHat)
  
  return(returnList)
}


###### example: for lowd setting
case = "homo" # "homo" or "hetero"
set.seed(0)
n = 500
VIO.str = 0.2
IV.str = 0.5
pi.value = IV.str*VIO.str
beta = 1
L = 10; px=10
s1 = 2; s2 = 4; s=s1+s2
alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/3)
gamma=rep(IV.str,L)
p=L+px # p stands for the number of total exogeneous variables 
phi<-rep(0,px)
psi<-rep(0,px)
phi[1:px]<-(1/px)*seq(1,px)+0.5
psi[1:px]<-(1/px)*seq(1,px)+1
rho=0.5
A1gen <- function(rho, p){
  A1 = matrix(0, nrow=p, ncol=p)
  for(i in 1:p) for(j in 1:p) A1[i, j] = rho^(abs(i-j))
  return(A1)
}
Cov<-(A1gen(rho,p))
W = mvrnorm(n, rep(0, p), Cov)
Z = W[, 1:L]
X = W[, (L+1):p]
if(case=="hetero"){
  epsilon1 = rnorm(n)
  tao1 = rep(NA, n); for(i.n in 1:n) tao1[i.n] = rnorm(n=1, mean=0, sd=0.25+0.5*(Z[i.n, 1])^2)
  tao2 = rnorm(n)
  epsilon2 = 0.3*epsilon1 + sqrt((1-0.3^2)/(0.86^4+1.38072^2))*(1.38072*tao1+0.86^2*tao2)
}else if(case=="homo"){
  epsilonSigma = matrix(c(1, 0.8, 0.8, 1), 2, 2)
  epsilon = mvrnorm(n, rep(0, 2), epsilonSigma)
  epsilon1 = epsilon[,1]
  epsilon2 = epsilon[,2]
}
D = 0.5 + Z %*% gamma+ X%*% psi + epsilon1
Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon2

out1 <- SearchingSampling(Y, D, Z, X, robust=TRUE, Sampling = FALSE)
out2 <- SearchingSampling(Y, D, Z, X, robust=TRUE, Sampling = TRUE)
out1$CI; out2$CI
