#' @title Two-Stage Hard Thresholding
#' @description Two-Stage Hard Thresholding main function, which provides the robust inference of the treatment effect in the presence of invalid instrumental variables in both low-dimensional and high-dimensional settings.
#'
#' @param Y continuous and non-missing, n by 1 numeric outcome vector.
#' @param D continuous or discrete, non-missing, n by 1 numeric treatment vector.
#' @param Z continuous or discrete, non-missing, n by p_z numeric instrument matrix, containing p_z instruments..
#' @param X optional,continuous or discrete, n by p_x numeric covariate matrix, containing p_z covariates.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param boot.SHat a boolean scalar indicating to implement bootstrap to get threshold for Shat, with default FALSE.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater 2, with default 2.01.
#' @param method a character scalar declaring the method used to estimate the inputs in TSHT, "OLS" works for ordinary least square and "DeLasso" works for high dimension. Default by "OLS".
#' @param voting a character scalar declaring the voting option used to estimate Vhat, 'MP' works for majority and plurality voting, 'MaxClique' works for finding maximal clique in the IV voting matrix, and 'Conservative' works for conservative voting procedure, with default MaxClique.
#'
#' @return
#'     \item{\code{betaHat}}{a numeric scalar denoting the estimate of treatment effect.}
#'     \item{\code{beta.sdHat}}{a numeric scalar denoting the estimated standard deviation of betaHat.}
#'     \item{\code{ci}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs.}
#'     \item{\code{voting.mat}}{a numeric matrix denoting the votes among the candidates of valid and relevant IVs with components 0 and 1.}
#'     \item{\code{beta.clique}}{a numeric matrix where each row represents the estiamted betahat corresponding to each maximum clique. Only returns when \code{voting} is \code{'MaxClique'}.}
#'     \item{\code{beta.sd.clique}}{a numeric matrix where each row represents the estimated variance of betahat corresponding to each maximum clique. Only returns when \code{voting} is \code{'MaxClique'}}
#'     \item{\code{CI.clique}}{a numeric matrix where each row represents the CI corresponding to each maximum clique. Only returns when \code{voting} is \code{'MaxClique'}}
#'     \item{\code{max.clique}}{a numeric matrix denoting maximum cliques of voted as valid and relevant IVs. Only returns when \code{voting} is \code{'MaxClique'}}
#'
#' @export
#'
#' @examples
#' \dontrun{
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
TSHT <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05, boot.SHat = FALSE ,tuning=2.01,
                 method="OLS", voting = 'MaxClique') {
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),is.vector(Y)||(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  if (is.vector(Y)) {
    Y <- cbind(Y)
  }
  Y = as.numeric(Y)

  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),is.vector(D)||(is.matrix(D) || is.data.frame(D)) && ncol(D) == 1)
  stopifnot(all(!is.na(D)))
  if (is.vector(D)) {
    D <- cbind(D)
  }
  D = as.numeric(D)

  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),(is.vector(Z) || is.matrix(Z)))
  stopifnot(all(!is.na(Z)))
  if (is.vector(Z)) {
    Z <- cbind(Z)
  }
  # Check dimesions
  stopifnot(nrow(Y) == nrow(D), nrow(Y) == nrow(Z))

  # Check X, if present
  if(!missing(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),(is.vector(X))||(is.matrix(X) && nrow(X) == nrow(Z)))
    stopifnot(all(!is.na(X)))
    if (is.vector(X)) {
      X <- cbind(X)
    }
    W = cbind(Z,X)
  } else {
    W = Z
  }

  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(is.logical(boot.SHat))
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  stopifnot(method=='OLS' | method=='DeLasso')
  stopifnot(voting=='MP' | voting=='MaxClique' | voting == 'Conservative')


  # Derive Inputs for TSHT
  n = length(Y); pz=ncol(Z)
  if(method == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept)
    A = t(inputs$WUMat) %*% inputs$WUMat / n
  } else if (method == "DeLasso") {
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)
    A = diag(pz)
  }

  # Reduced form estimator
  ITT_Y = inputs$ITT_Y;
  ITT_D = inputs$ITT_D;
  WUMat = inputs$WUMat;
  SigmaSqD = inputs$SigmaSqD;
  SigmaSqY = inputs$SigmaSqY;
  SigmaYD=inputs$SigmaYD;
  covW=inputs$covW

  # Estimate Valid IVs
  SetHats = TSHT.VHat(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,
                      SigmaSqD = SigmaSqD,SigmaSqY = SigmaSqY,SigmaYD=SigmaYD,
                      covW=covW,boot.SHat = boot.SHat, tuning=tuning, voting = voting)
  VHat = SetHats$VHat; SHat = SetHats$SHat


  if (voting == 'MaxClique') {
    max.clique <- SetHats$max.clique
    max.clique.mat <- matrix(0,nrow = length(max.clique),ncol = length(max.clique[[1]]))
    CI.temp <- matrix(0,nrow = length(max.clique), ncol = 2)
    beta.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    betavar.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    for (i in 1:length(max.clique)) {
      temp <- SHat[sort(as.numeric(max.clique[[i]]))]
      max.clique.mat[i,] <- temp
      AVHat = solve(A[temp,temp])
      betaHat = (t(ITT_Y[temp]) %*% AVHat %*% ITT_D[temp]) / (t(ITT_D[temp]) %*% AVHat %*% ITT_D[temp])
      SigmaSq = SigmaSqY + betaHat^2 * SigmaSqD - 2*betaHat * SigmaYD
      betaVarHat = SigmaSq * (t(ITT_D[temp]) %*% AVHat %*% (t(WUMat) %*% WUMat/ n)[temp,temp] %*% AVHat %*% ITT_D[temp]) / (t(ITT_D[temp]) %*% AVHat %*% ITT_D[temp])^2
      ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat / n),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat/n))
      CI.temp[i,] <- ci
      beta.temp[i,] <- betaHat
      betavar.temp[i,] <- betaVarHat
    }
    uni<- intervals::Intervals(CI.temp)
    ###### construct the confidence interval by taking a union
    CI.union<-as.matrix(intervals::interval_union(uni))
    # CI.union <- t(as.matrix(c(min(CI.union),max(CI.union)))) # added
    # ci <- as.vector(CI.union)
  }

  # Obtain point est, se, and ci
  AVHat = solve(A[VHat,VHat])
  betaHat = (t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat]) %*% AVHat %*% ITT_D[VHat])
  SigmaSq = SigmaSqY + betaHat^2 * SigmaSqD - 2*betaHat * SigmaYD
  betaVarHat = SigmaSq * (t(ITT_D[VHat]) %*% AVHat %*% (t(WUMat) %*% WUMat/ n)[VHat,VHat] %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat]) %*% AVHat %*% ITT_D[VHat])^2
  if (voting != 'MaxClique') {
    ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat / n),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat/n))
    return(list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat/n),ci=ci,SHat=SHat,VHat = VHat,voting.mat=SetHats$voting.mat))
  } else {
    return(list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat/n),ci=CI.union,SHat=SHat,VHat = VHat,voting.mat=SetHats$voting.mat,
                beta.clique = beta.temp,beta.sd.clique = sqrt(betavar.temp), CI.clique = CI.temp,
                max.clique = max.clique.mat))
  }

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
#'     \item{\code{covW}}{a numeric, non-missing matrix that computes the sample covariance of W}
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
  qrW = Matrix::qr(W)
  ITT_Y = Matrix::qr.coef(qrW,Y)[1:pz]
  ITT_D = Matrix::qr.coef(qrW,D)[1:pz]
  SigmaSqY = sum(Matrix::qr.resid(qrW,Y)^2)/(n -p)
  SigmaSqD = sum(Matrix::qr.resid(qrW,D)^2)/(n -p)
  SigmaYD = sum(Matrix::qr.resid(qrW,Y) * Matrix::qr.resid(qrW,D)) / (n - p)

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD,covW=covW))
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
#'     \item{\code{covW}}{a numeric, non-missing matrix that computes the sample covariance of W}
#' @export
#'
#' @importFrom flare slim
#' @importFrom glmnet glmnet
#' @importFrom stats qnorm coef
TSHT.DeLasso <- function(Y,D,W,pz,intercept=TRUE) {
  n = nrow(W)
  covW = t(W) %*% W /n #this should automatically turn covW into a matrix
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

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD,covW=covW))
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
#' @param covW a numeric, non-missing matrix that computes the sample covariance of W
#' @param boot.SHat a boolean scalar indicating to implement bootstrap to get threshold for Shat, with default FALSE.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater 2, with default 2.01.
#' @param voting a character scalar declaring the voting option used to estimate Vhat, 'MP' works for majority and plurality voting, 'MaxClique' works for finding maximal clique in the IV voting matrix, and 'Conservative' works for conservative voting procedure, with default MaxClique.
#'
#' @return
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs.}
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#'     \item{\code{max.clique}}{a numeric list denoting the maximum cliques of valid and relevant IVs. Only return when \code{voting} is \code{MaxClique}.}
#'     \item{\code{voting.mat}}{a numeric matrix denoting the votes among the candidates of valid and relevant IVs with components 0 and 1.}
#' @export
#'
#'
TSHT.VHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,covW,
                      boot.SHat = FALSE,tuning = 2.01,voting = 'MaxClique') {
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
  stopifnot(voting=='MP' | voting=='MaxClique' | voting == 'Conservative')

  # Constants
  n = nrow(WUMat);
  pz = length(ITT_Y)
  # First Stage
  if(boot.SHat==TRUE){
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
  diag(VHats.boot.sym) <- 1
  # Voting method
  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners

  if (voting == 'MaxClique') {
    voting.graph <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(VHats.boot.sym))
    max.clique <- igraph::largest.cliques(voting.graph)
    VHat <- unique(igraph::as_ids(Reduce(c,max.clique))) # take the union if multiple max cliques exist
    VHat <- sort(as.numeric(VHat))
  } else if (voting == 'MP') {
    VHat <- sort(as.numeric(union(VM.m,VM.p))) # Union of majority and plurality winners
  } else if (voting == 'Conservative'){
    V.set<-NULL
    for(index in VM.p){
      V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
    }
    VHat<-NULL
    for(index in V.set){
      VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
    }
    VHat=sort(as.numeric(VHat))
  }

  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }

  if (voting == 'MaxClique') {
    returnList <- list(SHat=SHat,VHat=VHat,max.clique=max.clique,voting.mat=VHats.boot.sym)
  } else {
    returnList <- list(SHat=SHat,VHat=VHat,voting.mat=VHats.boot.sym)
  }
  return(returnList)
}


cut.off.IVStr<-function(SigmaSqD,WUMat,pz,N=1000,cut.prob=0.99){
  n=nrow(WUMat)
  unit.matrix<-t(WUMat)%*%WUMat/n^2
  Cov.D<-SigmaSqD*unit.matrix
  max.vec<-rep(NA,N)
  for(j in 1:N){
    gamma.s<-MASS::mvrnorm(1, rep(0,pz), Cov.D)
    SE.norm<-diag(unit.matrix)^{1/2}
    max.vec[j]<-max(abs(gamma.s/SE.norm))
  }
  critical.val<-quantile(max.vec,probs=0.99)
  return(critical.val)
}
