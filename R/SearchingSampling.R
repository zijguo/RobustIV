##################################
### Date: 4/28/2022
### Author: Zhenyu WANG

library(SIHR)
library(MASS)
library(intervals)

############# Main Function #############
#' SearchingSampling
#' @description The proposed searching/sampling method
#'
#' @param Y outcome vector
#' @param D treatment vector
#' @param Z instruments
#' @param X covariates (default=\code{NULL})
#' @param intercept fit model with intercept or not (default=\code{TRUE})
#' @param lowd if \code{TRUE}, fit low-dimensional model, else fit high-dimensional one (default=\code{TRUE})
#' @param robust if \code{TRUE}, fit the model in heteroscedastic way (default=\code{TRUE})
#' @param CI.init initial interval for beta. If \code{NULL}, it will be generated automatically. (default=\code{NULL})
#' @param a grid size for constructing beta grids (default=0.6)
#' @param Sampling if \code{TRUE}, use the proposed sampling method; else use the proposed searching method. (default=\code{TRUE})
#' @param rho initial value constructed for sampling method (default=\code{NULL})
#' @param M sampling times. (default=1000)
#' @param prop proportion of intervals kept when sampling. (default=0.1)
#' @param filtering filtering sampling or not (default=\code{TRUE})
#'
#' @return
#' \item{CI}{confidence interval for beta}
#' \item{check}{the plurality rule being checked TRUE or FALSE emprically}
#' \item{VHat}{valid instruments}
#' \item{SHat}{relevant instruments}
#' @export
#' @import intervals MASS SIHR
#' @examples
SearchingSampling <- function(Y, D, Z, X=NULL, intercept=TRUE, 
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


########### Example: for lowd setting ########### 
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

########### Helpers function #############
grid.CI <- function(CI.matrix, grid.size){
  d = dim(CI.matrix)[1]
  grid.seq = NULL
  for(l in 1:d) grid.seq = c(grid.seq, seq(CI.matrix[l, 1], CI.matrix[l, 2], by=grid.size))
  return(grid.seq)
}

Searching.CI <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet, beta.grid){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V.Gamma)[1]
  
  ## new rho method
  Tn = qnorm(1-0.05/(2*pz))
  
  ## valid grid
  valid.grid = rep(NA, n.beta)
  for(j in 1:n.beta){
    b = beta.grid[j]
    temp = sqrt(diag(V.Gamma + b^2*V.gamma - 2*b*C)/n)
    valid.grid[j] = sum(abs(ITT_Y[InitiSet] - b*ITT_D[InitiSet]) < Tn*temp[InitiSet])
  }
  
  ## select beta
  if(length(beta.grid[which(valid.grid > threshold.size)])==0){
    rule=FALSE
    warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
    sel.index = which(valid.grid==max(valid.grid))
  }else{
    rule=TRUE
    sel.index = which(valid.grid>threshold.size)
  }
  CI = t(as.matrix(c(min(beta.grid[sel.index]), max(beta.grid[sel.index]))))
  
  return(list(CI=CI, rule=rule))
}

Searching.CI.sampling <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet,
                                  beta.grid, rho=NULL, M=1000, prop=0.1, filtering=TRUE){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V.Gamma)[1]
  
  ## new rho method
  Tn = qnorm(1-0.05/(2*pz))
  
  ## Covariance Matrix
  Cov1 = cbind(V.Gamma/n, C/n)
  Cov2 = cbind(t(C/n), V.gamma/n)
  Cov.total = rbind(Cov1, Cov2)
  
  valid.grid.sample = matrix(NA, nrow=M, ncol=n.beta)
  Gen.mat = MASS::mvrnorm(M, rep(0, 2*pz), Cov.total)
  if(is.null(rho)) rho = (log(n)/M)^(1/(2*length(InitiSet)))/6 # initial rho if not specified
  
  if(filtering){
    temp1 = abs(t(t(Gen.mat[,InitiSet])/sqrt(diag(V.Gamma)[InitiSet]/n)))
    temp2 = abs(t(t(Gen.mat[,pz+InitiSet])/sqrt(diag(V.gamma)[InitiSet]/n)))
    temp = cbind(temp1, temp2)
    temp = apply(temp, MARGIN=1, FUN=max)
    temp = temp <= qnorm(1-0.05/(4*length(InitiSet)))
    Gen.mat = Gen.mat[temp,]
    M = sum(temp)
  }
  
  while(rho < 0.5){
    for(m in 1:M){
      ITT_Y.sample = ITT_Y - Gen.mat[m, 1:pz]
      ITT_D.sample = ITT_D - Gen.mat[m, (pz+1):(2*pz)]
      for(j in 1:n.beta){
        b = beta.grid[j]
        temp = sqrt(diag(V.Gamma + b^2*V.gamma - 2*b*C)/n)
        valid.grid.sample[m, j] = sum(abs(ITT_Y.sample[InitiSet] - b*ITT_D.sample[InitiSet])<
                                        rho*Tn*temp[InitiSet])
      }
    }
    CI = matrix(NA, nrow=M, ncol=2)
    for(m in 1:M){
      if(length(which(valid.grid.sample[m, ] > threshold.size))>0){
        CI[m, 1] = min(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
        CI[m, 2] = max(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
      }
    }
    CI = CI[!rowSums(is.na(CI)), , drop=FALSE] # CI = na.omit(CI)
    
    ## check CI dims, stop iterations if CI dim is big enough
    if(dim(as.matrix(CI))[1] >= prop*M) break
    rho = 1.25 * rho # increase rho with iterations
  }
  
  rule = TRUE
  if(dim(as.matrix(CI))[1] < prop*M){
    warning("Sampling Criterion not met, trasfer to Searching Method.")
    CI.searching = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet, beta.grid)
    rule = CI.searching$rule
    CI = CI.searching$CI
  }else{
    uni = Intervals(CI)
    CI = as.matrix(interval_union(uni))
    CI = t(as.matrix(c(min(CI[,1]), max(CI[,2]))))
  }
  
  return(list(CI=CI, rule=rule))
}

TSHT.Init <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C){
  pz = nrow(V.Gamma)
  ## First Stage
  Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)/2))
  SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak. 
            TSHT with these IVs will give misleading CIs, SEs, and p-values. 
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE
  
  ## Second Stage
  nCand = length(SHat)
  VHats.bool = matrix(FALSE, nCand, nCand)
  colnames(VHats.bool) = rownames(VHats.bool) = SHat
  
  for(j in SHat){
    beta.j = ITT_Y[j]/ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    Temp = V.Gamma + beta.j^2*V.gamma - 2*beta.j*C
    SE.j = rep(NA, pz)
    for(k in 1:pz){
      SE.j[k] = 1/n * (Temp[k,k] + (ITT_D[k]/ITT_D[j])^2*Temp[j,j] - 
                         2*(ITT_D[k]/ITT_D[j])*Temp[k,j])
    }
    PHat.bool.j = abs(pi.j) <= sqrt(SE.j)*sqrt(log(n))
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat), as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  diag(VHats.boot.sym) = 1
  
  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
  
  V.set<-NULL
  for(index in union(VM.m,VM.p)){
    V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
  }
  VHat<-NULL
  for(index in V.set){
    VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
  }
  VHat=sort(as.numeric(VHat))
  
  out <- list(SHat=SHat, VHat=VHat, voting.mat=VHats.boot.sym)
  return(out)
}

Lasso.init <- function(X, y, lambda = "CV.min", intercept = FALSE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  }
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}
