

#' @title Two-Stage Hard Thresholding
#' @description Perform Two-Stage Hard Thresholding method, which provides the robust inference of the treatment effect in the presence of invalid instrumental variables.
#'
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param method The method used to estimate the reduced form parameters. \code{"OLS"} stands for ordinary least squares, \code{"DeLasso"} stands for the debiased Lasso estimator, and \code{"Fast.DeLasso"} stands for the debiased Lasso estimator with fast algorithm. (default = \code{"OLS"})
#' @param voting The voting option used to estimate valid IVs. \code{'MP'} stands for majority and plurality voting, \code{'MaxClique'} stands for finding maximal clique in the IV voting matrix, and \code{'Conservative'} stands for conservative voting procedure. Conservative voting is used to get an initial estimator of valid IVs in the Searching-Sampling method. (default= \code{'MaxClique'}).
#' @param robust If \code{TRUE}, the method is robust to heteroskedastic errors. If \code{FALSE}, the method assumes homoskedastic errors. (default = \code{FALSE})
#' @param alpha The significance level for the confidence interval. (default = \code{0.05})
#' @param tuning.1st tuning parameter used in 1st stage to select relevant instruments. If \code{NULL}, it will be generated data-dependently, see Details. (default=\code{NULL})
#' @param tuning.2nd tuning parameter used in 2nd stage to select valid instruments. If \code{NULL}, it will be generated data-dependently, see Details. (default=\code{NULL})
#'
#' @details When \code{robust = TRUE}, only \code{’OLS’} can be input to \code{method}.
#' When \code{voting = MaxClique} and there are multiple maximum cliques, \code{betaHat},\code{beta.sdHat},\code{ci}, and \code{VHat} will be list objects
#' where each element of list corresponds to each maximum clique.
#' As for tuning parameter in the 1st stage and 2nd stage, if do not specify, for method "OLS" we adopt \eqn{\sqrt{\log n}}, and for other methods
#' we adopt \eqn{\max{(\sqrt{2.01 \log p_z}, \sqrt{\log n})}}.
#
#' @return \code{TSHT} returns an object of class "TSHT", which is a list containing the following components:
#'     \item{\code{betaHat}}{The estimate of treatment effect.}
#'     \item{\code{beta.sdHat}}{The estimated standard error of \code{betaHat}.}
#'     \item{\code{ci}}{The 1-alpha confidence interval for \code{beta}.}
#'     \item{\code{SHat}}{The set of relevant IVs.}
#'     \item{\code{VHat}}{The set of relevant and valid IVs.}
#'     \item{\code{voting.mat}}{The voting matrix on whether the elements of each \code{SHat} are valid or not.}
#'     \item{\code{check}}{The indicator that the majority rule is satisfied.}
#' @export
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
#' Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), Confidence intervals for causal effects with invalid instruments by using two-stage hard thresholding with voting, \emph{J. R. Stat. Soc. B}, 80: 793-815. \cr
#' }
TSHT <- function(Y,D,Z,X,intercept=TRUE, method=c("OLS","DeLasso","Fast.DeLasso"),
                 voting = c('MaxClique','MP','Conservative'), robust = FALSE, alpha=0.05,
                 tuning.1st=NULL, tuning.2nd=NULL) {
  stopifnot(is.logical(robust))
  method = match.arg(method)
  if(method %in% c("DeLasso", "Fast.DeLasso") && robust==TRUE){
    robust = FALSE
    cat(sprintf("For methods %s, robust is set FALSE,
                as we only consider homoscedastic noise.\n", method))
  }
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

  # centralize W
  W = scale(W, center = T, scale = F)

  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(method=='OLS' | method=='DeLasso'|method =="Fast.DeLasso")
  stopifnot(voting=='MP' | voting=='MaxClique' | voting == 'Conservative')

  # Derive Inputs for TSHT
  n = length(Y); pz=ncol(Z)

  if(method == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept)
    A = t(inputs$WUMat) %*% inputs$WUMat / n
  } else if (method == "DeLasso") {
    inputs =  TSHT.SIHR(Y,D,W,pz,intercept)
    A = diag(pz)
  } else if (method =="Fast.DeLasso"){
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)
    A = diag(pz)
  }

  if (robust) {
    ITT_Y = inputs$ITT_Y;
    ITT_D = inputs$ITT_D;
    WUMat = inputs$WUMat;
    V.Gamma = inputs$V.Gamma
    V.gamma = inputs$V.gamma
    C = inputs$C
  } else {
    # Reduced form estimator
    ITT_Y = inputs$ITT_Y;
    ITT_D = inputs$ITT_D;
    WUMat = inputs$WUMat;
    SigmaSqD = inputs$SigmaSqD;
    SigmaSqY = inputs$SigmaSqY;
    SigmaYD=inputs$SigmaYD;
    V.Gamma = SigmaSqY * t(WUMat)%*%WUMat / n
    V.gamma = SigmaSqD * t(WUMat)%*%WUMat / n
    C = SigmaYD * t(WUMat)%*%WUMat / n
  }
  # Estimate Valid IVs
  SetHats = TSHT.VHat(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, voting, method=method, tuning.1st=tuning.1st, tuning.2nd=tuning.2nd)

  VHat = SetHats$VHat; SHat = SetHats$SHat
  check = T

  if(typeof(VHat)=="list"){
    if(length(VHat[[1]])< length(SHat)/2){
      cat('Majority rule fails.','\n')
      check=F
    }
    betaHat = beta.sdHat = ci = vector("list", length(VHat))
    # name of list
    names(VHat) = names(betaHat) = names(beta.sdHat) = names(ci) = paste0("MaxClique",1:length(VHat))
    for(i.VHat in 1:length(VHat)){
      VHat.i = VHat[[i.VHat]]
      AVHat = solve(A[VHat.i,VHat.i])
      betaHat.i = as.numeric((t(ITT_Y[VHat.i]) %*% AVHat %*% ITT_D[VHat.i]) / (t(ITT_D[VHat.i])  %*% AVHat %*% ITT_D[VHat.i]))
      temp <- (V.Gamma -2*betaHat.i*C+betaHat.i^2*V.gamma)[VHat.i,VHat.i]

      # 1-step iteration
      AVHat = solve(temp)
      betaHat.i = as.numeric((t(ITT_Y[VHat.i]) %*% AVHat %*% ITT_D[VHat.i]) / (t(ITT_D[VHat.i])  %*% AVHat %*% ITT_D[VHat.i]))
      temp <- (V.Gamma -2*betaHat.i*C+betaHat.i^2*V.gamma)[VHat.i,VHat.i]
      betaVarHat.i <- t(ITT_D[VHat.i])%*% AVHat %*%temp%*% AVHat %*%ITT_D[VHat.i] / (n*(t(ITT_D[VHat.i])%*% AVHat %*%ITT_D[VHat.i])^2)
      # U-2*betaHat*UV+betaHat^2*V
      ci.i = c(betaHat.i - qnorm(1-alpha/2) * sqrt(betaVarHat.i),betaHat.i + qnorm(1-alpha/2) * sqrt(betaVarHat.i))

      # store
      betaHat[[i.VHat]] = betaHat.i
      beta.sdHat[[i.VHat]] = as.numeric(sqrt(betaVarHat.i))
      ci[[i.VHat]] = ci.i
    }
    if(!is.null(colnames(Z))){
      SHat = colnames(Z)[SHat]
      VHat = lapply(VHat, FUN=function(x) colnames(Z)[x])
    }
  }else{

    if(length(VHat)< length(SHat)/2){
      cat('Majority rule fails.','\n')
      check=F
    }
    AVHat = solve(A[VHat,VHat])
    betaHat = as.numeric((t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat])  %*% AVHat %*% ITT_D[VHat]))
    temp <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[VHat,VHat]

    # 1-step iteration
    AVHat = solve(temp)
    betaHat = as.numeric((t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat])  %*% AVHat %*% ITT_D[VHat]))
    temp <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[VHat,VHat]
    betaVarHat <- t(ITT_D[VHat])%*% AVHat %*%temp%*% AVHat %*%ITT_D[VHat] / (n*(t(ITT_D[VHat])%*% AVHat %*%ITT_D[VHat])^2)
    # U-2*betaHat*UV+betaHat^2*V
    ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat))
    beta.sdHat = as.numeric(sqrt(betaVarHat))
    if (!is.null(colnames(Z))) {
      SHat = colnames(Z)[SHat]
      VHat = colnames(Z)[VHat]
    }
  }

  TSHTObject <- list(betaHat=betaHat,beta.sdHat=beta.sdHat,ci=ci,SHat=SHat,
                     VHat = VHat,voting.mat=SetHats$voting.mat,check = check,alpha=alpha)

  class(TSHTObject) <- 'TSHT'

  return(TSHTObject)
}


