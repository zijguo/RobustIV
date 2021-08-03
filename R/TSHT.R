#' @title Two-Stage Hard Thresholding
#' @description Two-Stage Hard Thresholding main function, which provides the robust inference of the treatment effect in the presence of invalid instrumental variables in both low-dimensional and high-dimensional settings.
#'
#' @param Y continuous and non-missing, n by 1 numeric outcome vector.
#' @param D continuous or discrete, non-missing, n by 1 numeric treatment vector.
#' @param Z continuous or discrete, non-missing, n by p_z numeric instrument matrix, containing p_z instruments..
#' @param X optional,continuous or discrete, n by p_x numeric covariate matrix, containing p_z covariates.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater 2, with default 2.01.
#' @param method a character scalar declaring the method used to estimate the inputs in TSHT, "OLS" works for ordinary least square and "DeLasso" works for high dimension. Default by "OLS".
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix, with default FALSE.
#'
#' @return
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs.}
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#'     \item{\code{betaHat}}{a numeric scalar denoting the estimate of treatment effect.}
#'     \item{\code{varHat}}{a numeric scalar denoting the estimated variance of betaHat.}
#'     \item{\code{ci}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#' @export
#'
#' @examples
#' \donttest{
#' ### Working Low Dimensional Example ###
#' library(RobustIV)
#' library(MASS)
#' n = 500; L = 10; s = 3
#' alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
#' epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#' Z = matrix(rnorm(n*L),n,L)
#' epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
#' D = 0.5 + Z %*% gamma + epsilon[,1]
#' Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]
#' TSHT(Y,D,Z)
#' TSHT(Y,D,Z,max_clique=TRUE)
#'
#'
#' ### Working High Dimensional Example ###
#' library(Matrix)
#' library(glmnet)
#' library(flare)
#' library(MASS)
#' n = 500; L = 600; s = 3; nRelevant = 10
#' alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,nRelevant),rep(0,L-nRelevant))
#' epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#' Z = matrix(rnorm(n*L),n,L)
#' epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
#' D =  0.5 + Z %*% gamma + epsilon[,1]
#' Y = -0.5 + Z %*% alpha + D * beta + epsilon[,2]
#' TSHT(Y,D,Z,method="DeLasso")
#' TSHT(Y,D,Z,method="DeLasso",max_clique=TRUE)
#' }
#'
#'
TSHT <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01,method="OLS",max_clique=FALSE) {
  # Check and Clean Input Type #
  # Check Y
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  Y = as.numeric(Y)

  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),(is.matrix(D) || is.data.frame(D)) && ncol(D) == 1)
  stopifnot(all(!is.na(D)))
  D = as.numeric(D)

  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),is.matrix(Z))
  stopifnot(all(!is.na(Z)))

  # Check dimesions
  stopifnot(length(Y) == length(D), length(Y) == nrow(Z))

  # Check X, if present
  if(!missing(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))

    W = cbind(Z,X)
  } else {
    W = Z
  }

  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  stopifnot(is.character(method))

  # Derive Inputs for TSHT
  n = length(Y); pz=ncol(Z)
  if(method == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept)
    A = t(inputs$WUMat) %*% inputs$WUMat / n
  } else if (method == "DeLasso") {
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)
    A = diag(pz)
  }

  # Estimate Valid IVs
  SetHats = TSHT.VHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                      SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,SigmaYD=inputs$SigmaYD,tuning=tuning,
                      max_clique = max_clique)
  VHat = SetHats$VHat; SHat = SetHats$SHat

  # Obtain point est, se, and ci
  AVHat = solve(A[VHat,VHat])
  betaHat = (t(inputs$ITT_Y[VHat]) %*% AVHat %*% inputs$ITT_D[VHat]) / (t(inputs$ITT_D[VHat]) %*% AVHat %*% inputs$ITT_D[VHat])
  SigmaSq = inputs$SigmaSqY + betaHat^2 * inputs$SigmaSqD - 2*betaHat * inputs$SigmaYD
  betaVarHat = SigmaSq * (t(inputs$ITT_D[VHat]) %*% AVHat %*% (t(inputs$WUMat) %*% inputs$WUMat/ n)[VHat,VHat] %*% AVHat %*% inputs$ITT_D[VHat]) / (t(inputs$ITT_D[VHat]) %*% AVHat %*% inputs$ITT_D[VHat])^2
  ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat / n),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat/n))

  return(list(VHat = VHat,SHat=SHat,betaHat=betaHat,betaVarHat = betaVarHat,ci=ci))
}


#' @title Two-Stage Hard Thresholding Ordinary Least Squares
#'
#' @description This function provides the estimates of the necessary inputs of TSHT method in a low-dimensional setting using Ordinary Least Squares.
#' @param Y continuous and non-missing, n by 1 numeric outcome vector.
#' @param D continuous or discrete, non-missing, n by 1 numeric treatment vector.
#' @param W continuous or discrete, non-missing, n by (p_z + p_x) numeric instrument-covariate matrix.
#' @param pz a numeric scalar denoting the number of instrument variables.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#'
#' @return
#'     \item{\code{ITT_Y}}{a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.}
#'     \item{\code{ITT_D}}{a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the outcome model.}
#'     \item{\code{WUMat}}{a numeric matrix denoting WU where U is the precision matrix of W and W is the instrument-covariate matrix.}
#'     \item{\code{SigmaSqY}}{a numeric scalar denoting the consistent estimator of the noise level in the outcome model.}
#'     \item{\code{SigmaSqY}}{a numeric scalar denoting the consistent estimator of the noise level in the treatment model.}
#'     \item{\code{SigmaYD}}{a numeric scalar denoting the consistent estimator of the covariance between the error term in the treatment model and the error term in the outcome model.}
#' @export
#'
TSHT.OLS <- function(Y,D,W,pz,intercept=TRUE) {
  # Include intercept
  if(intercept) {
    W = cbind(W,1)
  }
  p = ncol(W); n = nrow(W)

  # Compute covariance of W and W %*% U
  covW = t(W) %*% W /n #this should automatically turn covW into a matrix
  WUMat = W %*% (solve(covW))[,1:pz]

  # First Part (OLS Estimation)
  qrW = qr(W)
  ITT_Y = qr.coef(qrW,Y)[1:pz]
  ITT_D = qr.coef(qrW,D)[1:pz]
  SigmaSqY = sum(qr.resid(qrW,Y)^2)/(n -p)
  SigmaSqD = sum(qr.resid(qrW,D)^2)/(n -p)
  SigmaYD = sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / (n - p)

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}


#' @title Two-Stage Hard Thresholding Debiased LASSO
#'
#' @description This function provides the estimates of the necessary inputs of TSHT method in a high-dimensional setting using Debiased LASSO. This function shares the same inputs and outputs with the function \code{\link{TSHT.OLS}}.
#'
#' @param Y continuous and non-missing, n by 1 numeric outcome vector.
#' @param D continuous or discrete, non-missing, n by 1 numeric treatment vector.
#' @param W continuous or discrete, non-missing, n by (p_z + p_x) numeric instrument-covariate matrix.
#' @param pz a numeric scalar denoting the number of instrument variables.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#'
#' @return
#'     \item{\code{ITT_Y}}{a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.}
#'     \item{\code{ITT_D}}{a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.}
#'     \item{\code{WUMat}}{a numeric matrix denoting WU where U is the precision matrix of W and W is the instrument-covariate matrix.}
#'     \item{\code{SigmaSqY}}{a numeric scalar denoting the consistent estimator of the noise level in the outcome model.}
#'     \item{\code{SigmaSqY}}{a numeric scalar denoting the consistent estimator of the noise level in the treatment model.}
#'     \item{\code{SigmaYD}}{a numeric scalar denoting the consistent estimator of the covariance between the error term in the treatment model and the error term in the outcome model.}
#' @export
#'
#' @import Matrix
#' @importFrom flare slim
#' @importFrom glmnet glmnet
#' @importFrom stats qnorm coef
TSHT.DeLasso <- function(Y,D,W,pz,intercept=TRUE) {
  n = nrow(W)
  # Fit Reduced-Form Model for Y and D
  model_Y <- SSLasso(X=W,y=Y,intercept=intercept,verbose=FALSE)
  model_D = SSLasso(X=W,y=D,intercept=intercept,verbose=FALSE)
  ITT_Y = model_Y$unb.coef[1:(pz)]
  ITT_D = model_D$unb.coef[1:(pz)]
  resid_Y = model_Y$resid.lasso; resid_D = model_D$resid.lasso
  SigmaSqY=sum(resid_Y^2)/n
  SigmaSqD=sum(resid_D^2)/n
  SigmaYD =sum(resid_Y * resid_D)/n
  WUMat = model_D$WUMat[,1:(pz)]

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}


#' @title Two-Stage Hard Thresholding Instrumental Variable Selection
#'
#' @description Implementation of the Two-Stage Hard Thresholding with Voting procedures.This function only takes the necessary inputs of the functions \code{\link{TSHT.OLS}} or \code{\link{TSHT.DeLasso}}. Any other methods to compute the necessary inputs can be adopted by the users according to their preferences.
#'
#' @param ITT_Y a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.
#' @param ITT_D a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.
#' @param WUMat a numeric matrix denoting WU where U is the precision matrix of W and W is the instrument-covariate matrix (Z, X).
#' @param SigmaSqY a numeric scalar denoting the consistent estimator of the noise level in the outcome model.
#' @param SigmaSqD a numeric scalar denoting the consistent estimator of the noise level in the treatment model.
#' @param SigmaYD a numeric scalar denoting the consistent estimator of the covariance between the error term in the treatment model and the error term in the outcome model.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater 2, with default 2.01.
#' @param bootstrap a logical value, default by FALSE(What is this for in TSHT.Initial?).
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix.
#'
#' @return
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs.}
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#' @export
#'
#' @import igraph
#'
TSHT.VHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,bootstrap = FALSE,tuning = 2.01,max_clique) {
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)

  # Check WUMat
  stopifnot(!missing(WUMat),is.matrix(WUMat), nrow(WUMat) > 1, ncol(WUMat) == length(ITT_Y))
  stopifnot(all(!is.na(WUMat)))

  # Check Sigmas
  stopifnot(!missing(SigmaSqY), is.numeric(SigmaSqY), length(SigmaSqY) == 1, !is.na(SigmaSqY), SigmaSqY > 0)
  stopifnot(!missing(SigmaSqD), is.numeric(SigmaSqD), length(SigmaSqD) == 1, !is.na(SigmaSqD),SigmaSqD > 0)
  stopifnot(!missing(SigmaYD),is.numeric(SigmaYD), length(SigmaYD) == 1,!is.na(SigmaYD))

  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)

  # Constants
  n = nrow(WUMat);
  pz = length(ITT_Y)
  # First Stage
  if(bootstrap==TRUE){
    Tn<-min(cut.off.IVStr(SigmaSqD,WUMat,pz,cut.prob = 0.95),sqrt(log(n))) ### this can be modified by the user
    SE.norm<-(diag(solve(covW)/n)^{1/2})[1:pz]
    SHat<-(1:pz)[abs(ITT_D)>Tn*sqrt(SigmaSqD)*SE.norm]
  }else{
    SHat = (1:pz)[(abs(ITT_D) >= (sqrt(SigmaSqD * colSums(WUMat^2) /n) * sqrt(tuning*log(pz)/n)))]
  }
  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE

  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat
  for(j in SHat) {
    beta.j = ITT_Y[j] / ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums( (WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) *
      sqrt(tuning^2 * log(pz)/n)
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }

  # maximal clique
  if (max_clique) {
    voting.graph <- as.undirected(graph_from_adjacency_matrix(VHats.boot.sym))
    max_block <- largest_cliques(voting.graph)
    VHat <- as_ids(sort(max_block[[1]]))
    VHat <- as.numeric(VHat) # randomly pick the first one if multiple maximal cliques exist
  } else {
    # Voting
    #VM.1 = apply(VHats.bool,1,sum)
    #VM.2 = apply(VHats.bool,2,sum)
    #VM<-VM.1
    #for(l in 1:length(VM.1)){
    # VM[l]=min(VM.1[l],VM.2[l])
    #}
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
    VHat=as.numeric(VHat)
  }

  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }

  if (max_clique) {
    returnList <- list(SHat=SHat,VHat=VHat)
  } else {
    returnList <- list(SHat=SHat,VHat=VHat,V.set=V.set)
  }
  return(returnList)
}


cut.off.IVStr<-function(SigmaSqD,WUMat,pz,N=1000,cut.prob=0.99){
  n=nrow(WUMat)
  unit.matrix<-t(WUMat)%*%WUMat/n^2
  Cov.D<-SigmaSqD*unit.matrix
  max.vec<-rep(NA,N)
  for(j in 1:N){
    gamma.s<-mvrnorm(1, rep(0,pz), Cov.D)
    SE.norm<-diag(unit.matrix)^{1/2}
    max.vec[j]<-max(abs(gamma.s/SE.norm))
  }
  critical.val<-quantile(max.vec,probs=0.99)
  return(critical.val)
}
