#' @title Searching-Sampling
#' @description Construct Searching and Sampling confidence intervals for the causal effect, which provides the robust inference of the treatment effect in the presence of invalid instrumental variables in both low-dimensional and high-dimensional settings. It is robust to the mistakes in separating valid and invalid instruments.
#'
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param method The method used to estimate the reduced form parameters. \code{"OLS"} stands for ordinary least squares, \code{"DeLasso"} stands for the debiased Lasso estimator, and \code{"Fast.DeLasso"} stands for the debiased Lasso estimator with fast algorithm. (default = \code{"OLS"})
#' @param robust If \code{TRUE}, the method is robust to heteroskedastic errors. If \code{FALSE}, the method assumes homoskedastic errors.  (default = \code{FALSE})
#' @param Sampling If \code{TRUE}, use the proposed sampling method; else use the proposed searching method. (default=\code{TRUE})
#' @param alpha The significance level (default=\code{0.05})
#' @param CI.init An initial range for beta. If \code{NULL}, it will be generated automatically. (default=\code{NULL})
#' @param a Grid size for constructing beta grids. (default=\code{0.6})
#' @param rho Initial value constructed for sampling method. (default=\code{NULL})
#' @param M Re-sampling size. (default = \code{1000})
#' @param prop Proportion of non-empty intervals for sampling. (default=\code{0.1})
#' @param filtering Filtering the re-sampled data or not. (default=\code{TRUE})
#' @param tuning.1st tuning parameter used in 1st stage to select relevant instruments. If \code{NULL}, it will be generated data-dependently, see Details. (default=\code{NULL})
#' @param tuning.2nd tuning parameter used in 2nd stage to select valid instruments. If \code{NULL}, it will be generated data-dependently, see Details. (default=\code{NULL})
#'
#' @details When \code{robust = TRUE}, only \code{’OLS’} can be input to \code{method}. For \code{rho}, \code{M}, \code{prop}, and \code{filtering}, they are required only for \code{Sampling = TRUE}.
#' As for tuning parameter in the 1st stage and 2nd stage, if do not specify, for method "OLS" we adopt \eqn{\sqrt{\log n}}, and for other methods
#' we adopt \eqn{\max{(\sqrt{2.01 \log p_z}, \sqrt{\log n})}}.
#'
#' @return
#' \code{SearchingSampling} returns an object of class "SS", which is a list containing the following components:
#' \item{ci}{1-alpha confidence interval for beta.}
#' \item{SHat}{The set of relevant IVs.}
#' \item{VHat}{The initial set of relevant and valid IVs.}
#' \item{check}{The indicator that the plurality rule is satisfied.}

#' @export
#' @import intervals MASS CVXR glmnet
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
#' @references {
#' Guo, Z. (2021), Causal Inference with Invalid Instruments: Post-selection Problems and A Solution Using Searching and Sampling, Preprint \emph{arXiv:2104.06911}. \cr
#' }
SearchingSampling <- function(Y, D, Z, X=NULL, intercept=TRUE,
                              method=c("OLS","DeLasso","Fast.DeLasso"),
                              robust=FALSE, Sampling=TRUE, alpha=0.05,
                              CI.init = NULL, a=0.6,
                              rho=NULL, M=1000, prop=0.1, filtering=TRUE,
                              tuning.1st=NULL, tuning.2nd=NULL){
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

    out = TSHT.OLS(Y, D, W, pz, intercept=intercept)
    ITT_Y = out$ITT_Y
    ITT_D = out$ITT_D
    if(robust){
      V.Gamma = out$V.Gamma
      V.gamma = out$V.gamma
      C = out$C
    }else{
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
    SigmaSqY = out$SigmaSqY
    SigmaSqD = out$SigmaSqD
    SigmaYD = out$SigmaYD
    WUMat = out$WUMat
    V.Gamma = SigmaSqY * t(WUMat)%*%WUMat / n
    V.gamma = SigmaSqD * t(WUMat)%*%WUMat / n
    C = SigmaYD * t(WUMat)%*%WUMat / n

  }

  TSHT.out <- TSHT.VHat(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, voting="Conservative", method=method,
                        tuning.1st=tuning.1st, tuning.2nd=tuning.2nd)
  V0.hat = sort(TSHT.out$VHat)
  if(length(V0.hat)==0) stop("No valid IVs selected.")
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
  returnList <- list(ci=CI, check=rule, SHat=SHat, VHat=VHat)
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
  SS <- object
  cat("\nInitial set of Valid Instruments:", SS$VHat, "\n");
  if (SS$check) {
    cat("\nPlurality rule holds.\n")
  } else {
    cat("\nPlurality rule does not hold.\n")
  }
  cat(rep("_", 30), "\n")
  if (nrow(SS$ci)==1) {
    cat("\nConfidence Interval for Beta: [", SS$ci[1], ",", SS$ci[2], "]", "\n", sep = '');
  } else {
    cat("\nConfidence Intervals for Beta:\n")
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
    warning("Sampling Criterion not met, transfer to Searching Method.")
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

# TSHT.Init <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C){
#   pz = nrow(V.Gamma)
#   ## First Stage
#   Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)))
#   SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]
#   if(length(SHat)==0){
#     warning("First Thresholding Warning: IVs individually weak.
#             TSHT with these IVs will give misleading CIs, SEs, and p-values.
#             Use more robust methods.")
#     warning("Defaulting to treating all IVs as strong.")
#     SHat= 1:pz
#   }
#   SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE
#
#   ## Second Stage
#   nCand = length(SHat)
#   VHats.bool = matrix(FALSE, nCand, nCand)
#   colnames(VHats.bool) = rownames(VHats.bool) = SHat
#
#   for(j in SHat){
#     beta.j = ITT_Y[j]/ITT_D[j]
#     pi.j = ITT_Y - ITT_D * beta.j
#     Temp = V.Gamma + beta.j^2*V.gamma - 2*beta.j*C
#     SE.j = rep(NA, pz)
#     for(k in 1:pz){
#       SE.j[k] = 1/n * (Temp[k,k] + (ITT_D[k]/ITT_D[j])^2*Temp[j,j] -
#                          2*(ITT_D[k]/ITT_D[j])*Temp[k,j])
#     }
#     PHat.bool.j = abs(pi.j) <= sqrt(SE.j)*sqrt(log(n))
#     VHat.bool.j = PHat.bool.j * SHat.bool
#     VHats.bool[as.character(SHat), as.character(j)] = VHat.bool.j[SHat]
#   }
#   VHats.boot.sym<-VHats.bool
#   for(i in 1:dim(VHats.boot.sym)[1]){
#     for(j in 1:dim(VHats.boot.sym)[2]){
#       VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
#     }
#   }
#   diag(VHats.boot.sym) = 1
#
#   VM= apply(VHats.boot.sym,1,sum)
#   VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
#   VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
#
#   V.set<-NULL
#   for(index in union(VM.m,VM.p)){
#     V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
#   }
#   VHat<-NULL
#   for(index in V.set){
#     VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
#   }
#   VHat=sort(as.numeric(VHat))
#
#   out <- list(SHat=SHat, VHat=VHat, voting.mat=VHats.boot.sym)
#   return(out)
# }


