

#' @title Endogeneity-test
#' @description Conduct the endogeneity test in high-dimensional linear models in the presence of potentially invalid instrumental variables.
#'
#' @param Y A continuous vector of outcomes.
#' @param D A continuousvector of endogenous variables.
#' @param Z A matrix of instruments.
#' @param X A matrix of exogenous covariates.
#' @param intercept Should the intercept be included? Default is \code{TRUE} and if so, you do not need to add a column of 1s in X.
#' @param alpha The significance level for the confidence interval. (default = 0.05)
#' @param method The method which will be used to estimate the reduced form parameters in TSHT. "OLS" stands for ordinary least square and "DeLasso" stands for the debiased Lasso estimator. (default = "DeLasso")
#' @param invalid If \code{TRUE}, the method is robust to the presence of possibly invalid IVs; If \code{FALSE}, the method assumes all IVs to be valid. (default = TRUE)
#' @param voting The voting option used to estimate valid IVs. 'MP' stnads for majority and plurality voting, 'MaxClique' stands for maximum clique in the IV voting matrix. (default = 'MaxClique')
#'
#'
#' @return
#'     \code{endo.test} returns an object of class "endotest".
#'     An object class "endotest" is a list containing the following components:
#'    \item{\code{Q}}{Endogeneity test statistic.}
#'    \item{\code{Sigma12}}{Estimated covaraince.}
#'    \item{\code{VHat}}{The estimated set of relevant and vaild IVs.}
#'    \item{\code{alpha}}{The significance level.}
#' @export
#'
#' @examples
#' \dontrun{
#' n = 500; L = 600; s = 3; k = 10; px = 10;
#' alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,k),rep(0,L-k))
#' phi<-(1/px)*seq(1,px)+0.5; psi<-(1/px)*seq(1,px)+1
#' epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#' Z = matrix(rnorm(n*L),n,L)
#' X = matrix(rnorm(n*px),n,px)
#' epsilon = MASS::mvrnorm(n,rep(0,2),epsilonSigma)
#' D =  0.5 + Z %*% gamma + X %*% psi + epsilon[,1]
#' Y = -0.5 + Z %*% alpha + D * beta + X %*% phi + epsilon[,2]
#' TSHT(Y,D,Z,X, method = "DeLasso")
#' endo.test.model <- endo.test(Y,D,Z,X)
#' summary(endo.test.model)
#' }
#'
#'
#'
endo.test <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05, method="DeLasso",
                      invalid=TRUE, voting = 'MaxClique'){
  # Check and Clean Input Type #
  # Check Y

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
  stopifnot(method=='OLS' | method=='DeLasso')

  # Derive Inputs for Endogeneity test
  n = length(Y); pz=ncol(Z)
  if(method == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept) # recommend when n>>p

  } else if (method == "DeLasso"){
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)

  }

  # Estimate Relevant IVs

  if (invalid) {
    SetHats <- TSHT.VHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                         SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,
                         SigmaYD=inputs$SigmaYD,covW=inputs$covW, voting = voting)
    Set = SetHats$VHat
  } else {
    SetHats <- endo.SHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                         SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,
                         SigmaYD=inputs$SigmaYD,covW=inputs$covW)
    Set = SetHats
  }
  WUMat = inputs$WUMat

  # Obtain point est and our test statistic Q
  betaHat = t(inputs$ITT_Y[Set]) %*% inputs$ITT_D[Set] / (t(inputs$ITT_D[Set])  %*% inputs$ITT_D[Set])
  Sigma11 = inputs$SigmaSqY + betaHat^2 * inputs$SigmaSqD - 2*betaHat * inputs$SigmaYD # Sigma11.hat
  Sigma12 = inputs$SigmaYD-betaHat*inputs$SigmaSqD # Sigma12.hat

  Var1 = Sigma11 *(t(inputs$ITT_D[Set]) %*%  ((t(inputs$WUMat) %*% inputs$WUMat/ n)[Set,Set]) %*% inputs$ITT_D[Set]) / (t(inputs$ITT_D[Set]) %*% inputs$ITT_D[Set])^2
  Var2 = inputs$SigmaSqY*inputs$SigmaSqD-inputs$SigmaYD^2+2*(betaHat*inputs$SigmaSqD-inputs$SigmaYD)^2  # VarThetahat

  VarSig12 = inputs$SigmaSqD^2*Var1+Var2
  Q = sqrt(n)*Sigma12/sqrt(VarSig12) # our test statistic

  if (!is.null(colnames(Z))) {
    VHat = colnames(Z)[Set]
  }
  endo.test.model <- list(Q=Q,Sigma12=Sigma12,VHat=Set,alpha = alpha)
  class(endo.test.model) <- "endotest"
  return(endo.test.model)

}
#' Summary of endotest
#'
#' @param object endotest object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.endotest<- function(object,...){
  return(object)
}
#' Summary of endotest
#'
#' @param object endotest object
#' @param ...
#' @keywords internal
#' @return
#' @export
print.endotest<- function(x,...){
  endotest <- x
  cat("\nValid Instruments:", endotest$VHat, "\n");
  cat(rep("_", 30), "\n")
  cat("Test statistics Q = ",endotest$Q,"\n")
  if (abs(endotest$Q)>qnorm(1-endotest$alpha/2)) {
    cat("'H0 : Sigma12 = 0' is rejected at the significance level",endotest$alpha,".\n")
  } else {
    cat("'H0 : Sigma12 = 0' is not rejected at the significance level",endotest$alpha,".\n")
  }
  cat("P-value = ",1-pnorm(abs(endotest$Q)),"\n")
  cat("Estimated covariance:",endotest$Sigma12,"\n");
}

endo.SHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,covW) {
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


  # Constants
  n = nrow(WUMat);
  pz = length(ITT_Y)
  # First Stage
  Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)/2))

  SHat = (1:pz)[(abs(ITT_D) >= (sqrt(SigmaSqD * colSums(WUMat^2) /n) * sqrt(Tn^2/n)))]

  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  return(SHat)
}
