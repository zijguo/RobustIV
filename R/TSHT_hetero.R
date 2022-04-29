################ Source code ###################
TSHT_hetero <- function(Y, D, Z, X, intercept=TRUE, method=1){
  Y = as.numeric(Y)
  D = as.numeric(D)
  if(!missing(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))
    
    W = cbind(Z,X)
  } else {
    W = Z
  }
  # Derive Inputs for TSHT, consider only OLS setting now
  n = length(Y); pz=ncol(Z)
  inputs = TSHT.OLS_hetero(Y, D, W, pz, intercept)
  
  # Estimate Valid IVs
  SetHats = TSHT.VHat_hetero(n, inputs$ITT_Y, inputs$ITT_D, 
                             inputs$V.Gamma, inputs$V.gamma, inputs$C)
  VHat = SetHats$VHat; SHat = SetHats$SHat
  
  # Obtain point est, se and ci
  ### \Zhenyu{Currently, there are two methods to compute the result in hetero,
  ###         the first one is from CIIV package, the 2nd one is to utilize ivmodel}
  ### Method-1: CIIV
  if(method==1){
    out = CIIV_robust_method(Y, regressor=cbind(D, X, Z[,-VHat]), Exogenous = cbind(Z, X))
    CI = out$CI
    point = out$point
  }
  ### Method-2: ivmodel
  if(method==2){
    out = ivmodel(Y, D, Z[,VHat], cbind(Z[,-VHat], X), intercept=intercept, heteroSE=TRUE)
    CI = confint(out)["TSLS",]
    point = coef(out)["TSLS", "Estimate"]
  }
  return(list(point=point, CI=CI, VHat=VHat, SHat=SHat))
}

TSHT.OLS_hetero <- function(Y, D, W, pz, intercept=TRUE){
  n = nrow(W)
  if(intercept) W = cbind(W, 1)
  p = ncol(W)
  covW = t(W)%*%W/n
  U = solve(covW) # precision matrix
  WUMat = (W%*%U)[,1:pz]
  ## OLS estimators
  qrW = qr(W)
  ITT_Y = qr.coef(qrW, Y)[1:pz]
  ITT_D = qr.coef(qrW, D)[1:pz]
  resid_Y = as.vector(qr.resid(qrW, Y))
  resid_D = as.vector(qr.resid(qrW, D))
  V.Gamma = (t(WUMat)%*%diag(resid_Y^2)%*%WUMat)/n
  V.gamma = (t(WUMat)%*%diag(resid_D^2)%*%WUMat)/n
  C = (t(WUMat)%*%diag(resid_Y * resid_D)%*%WUMat)/n
  
  out <- list(ITT_Y = ITT_Y, 
              ITT_D = ITT_D, 
              V.Gamma = V.Gamma,
              V.gamma = V.gamma,
              C = C)
  return(out)
}

TSHT.VHat_hetero <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C){
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
  VHat = as.numeric(union(VM.m, VM.p))
  
  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }
  return(list(VHat=VHat, SHat=SHat))
}

CIIV_robust_method <- function(Y, regressor, Exogenous){
  length_regressor = ncol(regressor)
  Coefficients_CIM_GMM <- (CIM.HansenJTest(Y, regressor, Exogenous)[[2]])[1:length_regressor]
  sd_CIM_GMM <- (CIM.HansenJTest(Y, regressor, Exogenous)[[3]])[1:length_regressor];
  alpha = 0.05
  ci_CIM_GMM <- c(
    Coefficients_CIM_GMM[1] - qnorm(1-alpha/2) * sd_CIM_GMM[1],
    Coefficients_CIM_GMM[1] + qnorm(1-alpha/2) * sd_CIM_GMM[1]
  )
  return(list(point=Coefficients_CIM_GMM[1], CI = ci_CIM_GMM))
}

## Functions from the source code of CIIV
# Function for Two Step GMM and Hansen-J Test
CIM.HansenJTest <- function(Y,X,Z){
  res_FirstStep <- residuals(AER::ivreg(Y ~ X - 1 | Z));
  
  Weight_SecondStep <- crossprod(res_FirstStep * Z);
  
  Coef_SecondStep <- solve(
    t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
  ) %*% t(X) %*% Z %*% solve(Weight_SecondStep) %*% t(Z) %*% Y;
  
  res_SecondStep <- as.vector(Y - X %*% Coef_SecondStep);
  
  sd_SecondStep <- sqrt(diag(solve(
    t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
  ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)%*%crossprod(res_SecondStep * Z)%*%t(
    solve(
      t(X) %*% Z %*% solve(Weight_SecondStep) %*%t(Z) %*% X
    ) %*% t(X) %*% Z %*% solve(Weight_SecondStep)
  )));
  
  HansenJ_Stat <- t(res_SecondStep) %*% Z %*% solve(Weight_SecondStep) %*%
    t(Z) %*% res_SecondStep;
  
  list(HansenJ_Stat, Coef_SecondStep,sd_SecondStep)
}


################ Example ###################
library(MASS)
library(intervals)
library(ivmodel)
source("TSHT_hetero.R")

n=1000
VIO.str = 0.4
IV.str = 0.5
pi.value = IV.str*VIO.str
beta = 1

case = "hetero"
px = 10
L = 6; 
s1 = s2 = s3 = s4 = 1; s=s1+s2+s3+s4
alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(0.8,s2),-seq(0.4,s3),seq(0.6,s4))

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

set.seed(0)
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

out1 = TSHT_hetero(Y, D, Z, X, intercept=TRUE, method=1)
out1$point
out1$CI

out2 = TSHT_hetero(Y, D, Z, X, intercept=TRUE, method=2)
out2$point
out2$CI