############# Main Function #############
#' SearchingSampling
#' @description Conduct Searching-Sampling method, which can construct uniformly valid confidence intervals for the causal effect, which are robust to the mistakes in separating valid and invalid instruments.
#'
#' @param Y A vector of outcomes.
#' @param D A continuous vector of endogenous variables.
#' @param Z A matrix of instruments.
#' @param X A matrix of exogenous covariates.
#' @param intercept Should the intercept be included? Default is \code{TRUE} and if so, you do not need to add a column of 1s in X.
#' @param method 
#' @param robust If \code{TRUE}, fit the model in heteroscedastic way (default=\code{TRUE})
#' @param Sampling If \code{TRUE}, use the proposed sampling method; else use the proposed searching method. (default=\code{TRUE})
#' @param alpha Significance level (default=0.05)
#' @param CI.init Initial interval for beta. If \code{NULL}, it will be generated automatically. (default=\code{NULL})
#' @param a Grid size for constructing beta grids (default=0.6)
#' @param rho Initial value constructed for sampling method (default=\code{NULL}). It works only for Sampling method.
#' @param M Sampling times (default = 1000). It works only for Sampling method.
#' @param prop Proportion of intervals kept when sampling (default=0.1). It works only for Sampling method.
#' @param filtering Filtering sampling or not (default=\code{TRUE}). It works only for Sampling method.
#'
#' @return
#' \item{ci}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#' \item{check}{True or False indicating whether the plurality rule  is satisfied or not emprically.}
#' \item{VHat}{a numeric vector denoting the set of valid and relevant IVs.}
#' \item{SHat}{a numeric vector denoting the set of relevant IVs.}
#' @export
#' @import intervals MASS SIHR
#' @examples
#'\dontrun{
#' Y <- mroz[,"lwage"]
#' D <- mroz[,"educ"]
#' Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc","exper","expersq")])
#' X <- mroz[,"age"]
#' Searching.model <- SearchingSampling(Y,D,Z,X, Sampling = FALSE)
#' summary(Searching.model)
#' SS.model <- SearchingSampling(Y,D,Z,X)
#' summary(SS.model)
#'
#' }

SearchingSampling <- function(Y, D, Z, X=NULL, intercept=TRUE,
                              method=c("OLS","DeLasso","Fast.DeLasso"),
                              robust=TRUE, Sampling=TRUE, alpha=0.05,
                              CI.init = NULL, a=0.6, 
                              rho=NULL, M=1000, prop=0.1, filtering=TRUE){
  method = match.arg(method)
  if(method %in% c("DeLasso", "Fast.DeLasso") && robust==TRUE){
    robust = FALSE
    cat(sprintf("For methods %s, robust is set FALSE, 
                as we only consider homoscedastic noise.\n", method))
  }

  if(is.null(X)) W = Z else W = cbind(Z, X)
  n = length(Y); pz = ncol(Z); p = ncol(W)

  ## centralize W
  W = scale(W, center=T, scale=F)
  
  if(method=="OLS"){
    
    if(robust){
      out = TSHT.OLS_hetero(Y, D, W, pz, intercept=intercept)
      ITT_Y = out$ITT_Y
      ITT_D = out$ITT_D
      V.Gamma = out$V.Gamma
      V.gamma = out$V.gamma
      C = out$C
    }else{
      out = TSHT.OLS(Y, D, W, pz, intercept=intercept)
      ITT_Y = out$ITT_Y
      ITT_D = out$ITT_D
      SigmaSqY = out$SigmaSqY
      SigmaSqD = out$SigmaSqD
      SigmaYD = out$SigmaYD
      WUMat = out$WUMat
      V.Gamma = SigmaSqY * t(WUMat)%*%WUMat / n
      V.gamma = SigmaSqD * t(WUMat)%*%WUMat / n
      C = SigmaYD * t(WUMat)%*%WUMat / n
    }
    
  }else if(method=="Fast.DeLasso"){
    
    out = TSHT.DeLasso(Y, D, W, pz, intercept=intercept)
    ITT_Y = out$ITT_Y
    ITT_D = out$ITT_D
    SigmaSqY = out$SigmaSqY
    SigmaSqD = out$SigmaSqD
    SigmaYD = out$SigmaYD
    WUMat = out$WUMat
    V.Gamma = SigmaSqY * t(WUMat)%*%WUMat / n
    V.gamma = SigmaSqD * t(WUMat)%*%WUMat / n
    C = SigmaYD * t(WUMat)%*%WUMat / n
    
  }else if(method=="DeLasso"){
    
    out = TSHT.SIHR(Y, D, W, pz, intercept=intercept)
    ITT_Y = out$ITT_Y
    ITT_D = out$ITT_D
    V.Gamma = out$V.Gamma
    V.gamma = out$V.gamma
    C = out$C
  }

  TSHT.out <- TSHT.Init(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C)
  V0.hat = sort(TSHT.out$VHat)
  ## Construct range [L, U]
  if(is.vector(CI.init)){
    CI.init.union = matrix(CI.init, ncol=2)
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
                                        beta.grid = beta.grid, alpha=alpha, rho=rho, M=M, prop=prop, filtering=filtering)
    CI=CI.sampling$CI
    rule=CI.sampling$rule
  }else{
    ## Searching Method
    CI.searching = Searching.CI(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet = V0.hat,
                                beta.grid = beta.grid, alpha=alpha)
    CI=CI.searching$CI
    rule=CI.searching$rule
  }
  VHat=V0.hat; SHat=TSHT.out$SHat;
  if (!is.null(colnames(Z))) {
    SHat = colnames(Z)[SHat]
    VHat = colnames(Z)[VHat]
  }
  returnList <- list(ci=CI, check=rule, VHat=VHat, SHat=SHat)
  class(returnList) <- "SS"
  return(returnList)
}
#' Summary of SS
#'
#' @param object SS object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.SS<- function(object,...){
  return(object)
}
#' Summary of SS
#'
#' @param x SS object
#' @param ...
#' @keywords internal
#' @return
#' @export
print.SS<- function(x,...){
  SS <- x
  cat("\nInitial set of Valid Instruments:", SS$VHat, "\n");
  cat("\nRelevant Instruments:", SS$SHat,"\n","\nThus, Majority rule",ifelse(SS$check,"holds.","does not hold."), "\n");
  cat(rep("_", 30), "\n")
  if (nrow(SS$ci)==1) {
    cat("\nConfidence Interval of BetaHat: [", SS$ci[1], ",", SS$ci[2], "]", "\n", sep = '');
  } else {
    cat("\nConfidence Interval of BetaHat:\n")
    for (i in 1:nrow(SS$ci)) {
      cat("[", SS$ci[i,1], ",", SS$ci[i,2], "]", "\n", sep = '')
    }
  }

}


########### Helpers function #############
grid.CI <- function(CI.matrix, grid.size){
  d = dim(CI.matrix)[1]
  grid.seq = NULL
  for(l in 1:d) grid.seq = c(grid.seq, seq(CI.matrix[l, 1], CI.matrix[l, 2], by=grid.size))
  return(grid.seq)
}

Searching.CI <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, InitiSet, beta.grid, alpha=0.05){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V.Gamma)[1]

  ## new rho method
  Tn = qnorm(1-alpha/(2*pz))

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
                                  beta.grid, alpha=0.05, rho=NULL, M=1000, prop=0.1, filtering=TRUE){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V.Gamma)[1]

  ## new rho method
  Tn = qnorm(1-alpha/(2*pz))

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

TSHT.SIHR <- function(Y, D, W, pz, method="OLS", intercept=TRUE){
  
  init_Y = Lasso(W, Y, lambda="CV.min", intercept=intercept)
  init_D = Lasso(W, D, lambda="CV.min", intercept=intercept)
  
  W_int = ifelse(intercept, cbind(1, W), W)
  resid_Y = as.vector(Y - W_int%*%init_Y)
  resid_D = as.vector(Y - W_int%*%init_D)
  
  loading.mat = matrix(0, nrow=ncol(W), ncol=pz)
  for(i in 1:pz) loading.mat[i, i] = 1
  out1 = LF(W, Y, loading.mat, model="linear", intercept=intercept, intercept.loading=FALSE, verbose=FALSE)
  out2 = LF(W, D, loading.mat, model="linear", intercept=intercept, intercept.loading=FALSE, verbose=FALSE)
  ITT_Y = out1$est.debias.vec
  ITT_D = out2$est.debias.vec
  U = out2$proj.mat
  WUMat = W_int%*%U
  
  SigmaSqY = sum(resid_Y^2)/n
  SigmaSqD = sum(resid_D^2)/n
  SigmaYD = sum(resid_Y * resid_D)/n
  Temp = t(WUMat)%*%WUMat / n
  V.Gamma = SigmaSqY * Temp
  V.gamma = SigmaSqD * Temp
  C = SigmaYD * Temp
  
  return(list(ITT_Y = ITT_Y,
              ITT_D = ITT_D,
              V.Gamma = V.Gamma,
              V.gamma = V.gamma,
              C = C))
}
