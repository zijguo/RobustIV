
library(MASS)
library(intervals)
library(ivmodel)
################ Source code ###################
TSHT_hetero <- function(Y, D, Z, X, intercept=TRUE, method="OLS", tuning = 2.01){
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
  check = T
  if(length(VHat)< length(SHat)/2){
    cat('Majority rule fails.','\n')
    check=F
  }
  # Obtain point est, se and ci
  ITT_Y = inputs$ITT_Y;
  ITT_D = inputs$ITT_D;
  WUMat = inputs$WUMat;
  V.Gamma = inputs$V.Gamma
  V.gamma = inputs$V.gamma
  C = inputs$C
  
  if(method == "OLS") {
    A = t(WUMat) %*% WUMat / n
    AVHat = solve(A[VHat,VHat])
    betaHat = as.numeric((t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat])  %*% AVHat %*% ITT_D[VHat]))
    temp <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[VHat,VHat]
    AVHat = solve(temp)
    betaHat = as.numeric((t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat])  %*% AVHat %*% ITT_D[VHat]))
    temp <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[VHat,VHat]
    betaVar <- t(ITT_D[VHat])%*% AVHat %*%temp%*% AVHat %*%ITT_D[VHat] / (n*(t(ITT_D[VHat])%*% AVHat %*%ITT_D[VHat])^2)
    # U-2*betaHat*UV+betaHat^2*V
    ci = c(betaHat - qnorm(1-0.05/2) * sqrt(betaVar),betaHat + qnorm(1-0.05/2) * sqrt(betaVar))
    
  } else if (method == "DeLasso") {
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)
    A = diag(pz)
  }
  
  
  
  ### \Zhenyu{Currently, there are two methods to compute the result in hetero,
  ###         the first one is from CIIV package, the 2nd one is to utilize ivmodel}
  ### Method-1: CIIV
  if(method==1){
    out = CIIV_robust_method(Y, regressor=cbind(D, X, Z[,-VHat]), Exogenous = cbind(Z, X))
    ci = out$CI
    betaHat = out$point
  } else{
    out = ivmodel(Y, D, Z[,VHat], cbind(Z[,-VHat], X), intercept=intercept, heteroSE=TRUE, k = 1)
  }
  ### Method-2: ivmodel(TSLS)
  if(method==2){
    ci = confint(out)["TSLS",]
    betaHat = coef(out)["TSLS", "Estimate"]
  }
  ### Method-3: ivmodel(LIML)
  if (method==3) {
    ci = confint(out)["LIML",]
    betaHat = coef(out)["LIML", "Estimate"]
  }
  return(list(betaHat=betaHat,  ci=ci, VHat=VHat, SHat=SHat))
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
              C = C,
              WUMat = WUMat)
  return(out)
}

TSHT.VHat_hetero <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, tuning = 2.01){
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
    
    PHat.bool.j = abs(pi.j) <= sqrt(SE.j)*sqrt(tuning^2*log(pz))
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

