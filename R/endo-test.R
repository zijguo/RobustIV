

#' @title Endogeneity test in high dimensions
#' @description Conduct the endogeneity test with with high dimensional and possibly invalid instrumental variables.
#'
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param invalid If \code{TRUE}, the method is robust to the presence of possibly invalid IVs; If \code{FALSE}, the method assumes all IVs to be valid. (default = \code{FALSE})
#' @param method The method used to estimate the reduced form parameters. \code{"OLS"} stands for ordinary least squares, \code{"DeLasso"} stands for the debiased Lasso estimator, and \code{"Fast.DeLasso"} stands for the debiased Lasso estimator with fast algorithm. (default = \code{"Fast.DeLasso"})
#' @param voting The voting option used to estimate valid IVs. \code{'MP'} stnads for majority and plurality voting, \code{'MaxClique'} stands for maximum clique in the IV voting matrix. (default = \code{'MaxClique'})
#' @param alpha The significance level for the confidence interval. (default = \code{0.05})
#'
#'
#' @return
#'     \code{endo.test} returns an object of class "endotest", which is a list containing the following components:
#'    \item{\code{Q}}{Test statistic.}
#'    \item{\code{Sigma12}}{Estimated covaraince of the regression errors.}
#'    \item{\code{VHat}}{The set of vaild IVs.}
#'    \item{\code{p.value}}{The p-value of the endogeneity test.}
#'    \item{\code{check}}{Indicator for whether \eqn{H_0:\Sigma_{12}=0} is rejected or not.}
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
#' endo.test.model <- endo.test(Y,D,Z,X)
#' summary(endo.test.model)
#' }
#'
#' @references {
#' Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), Testing endogeneity with high dimensional covariates, \emph{Journal of Econometrics}, Elsevier, vol. 207(1), pages 175-187. \cr
#' }
#'
#'
#'
endo.test <- function(Y,D,Z,X,intercept=TRUE,invalid=FALSE, method=c("Fast.DeLasso","DeLasso","OLS"),
                       voting = c('MaxClique','MP'), alpha=0.05){
  # Check and Clean Input Type #
  # Check Y
  method = match.arg(method)
  voting = match.arg(voting)
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
  stopifnot(method=='OLS' | method=='DeLasso' | method =="Fast.DeLasso")
  stopifnot(voting=='MaxClique' | voting=='MP')

  # Derive Inputs for Endogeneity test
  n = length(Y); pz=ncol(Z)
  if(method == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept) # recommend when n>>p

  } else if (method == "Fast.DeLasso"){
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)

  } else if(method == "DeLasso"){
    inputs = TSHT.SIHR(Y, D, W, pz, intercept=intercept)
  }

  # Estimate Relevant IVs
  voting = match.arg(voting)
  if (invalid) {
    SetHats <- TSHT.VHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                         SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,
                         SigmaYD=inputs$SigmaYD, voting = voting)
    Set = SetHats$VHat
  } else {
    SetHats <- endo.SHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                         SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,
                         SigmaYD=inputs$SigmaYD)
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
  p.value <- 1-pnorm(abs(Q))
  check <- (abs(Q)>qnorm(1-alpha/2))
  endo.test.model <- list(Q=Q,Sigma12=Sigma12,VHat=Set,p.value = p.value,check = check,alpha = alpha)
  class(endo.test.model) <- "endotest"
  return(endo.test.model)

}
