

#' @title Endogeneity-test
#' @description Endogneity test function, which provides the robust inference of the treatment effect in the presence of invalid instrumental variables in both low-dimensional and high-dimensional settings.
#'
#' @param Y continuous and non-missing, n by 1 numeric outcome vector.
#' @param D continuous or discrete, non-missing, n by 1 numeric treatment vector.
#' @param Z continuous or discrete, non-missing, n by p_z numeric instrument matrix, containing p_z instruments.
#' @param X optional,continuous or discrete, n by p_x numeric covariate matrix, containing p_z covariates.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater 2, with default 2.01.
#' @param method a character scalar declaring the method used to estimate the inputs in TSHT, "OLS" works for ordinary least square and "DeLasso" works for high dimension.
#' @param invalid a boolean scalar asking to assume that there are some invalid instrument variables with TRUE/FALSE (default = TRUE)
#'
#' @return
#'    \item{\code{VHat}}{numeric vector : the estimated set of relevant and vaild IVs}
#'    \item{\code{Q}}{numeric value : our endogeneity test statistic}
#'    \item{\code{Sigma12}}{numeric value : estimated covaraince}
#' @export
#'
#' @examples
#'
#'
#' ### Working Low Dimensional Example ###
#'
#' library(MASS)
#' library(RobustIV)
#'
#' ## Generate the data
#' n <- 1000; pz <- 9; px <- 5;
#' p <- pz+px
#' s = 10; nRelevant = 7
#' beta <- 1
#' phi <- seq(0.6,1.0,length.out=px)
#' psi <- seq(1.1,1.5,length.out=px)
#' gamma = c(rep(1,nRelevant),rep(0,pz-nRelevant))
#' epsilonSigma = matrix(c(1.5,0.5,0.5,1.5),2,2)
#' W <- matrix(rnorm(n*p),n,p)
#' Z <- W[,1:pz]; X <- W[,(pz+1):p]
#' epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
#' D <- 1 + Z%*%gamma + X%*%psi + epsilon[,1]
#' Y <- -1 + D*beta + X%*%phi + epsilon[,2]
#'
#' ## Implement the endogeneity test
#'
#' endo.test(Y,D,Z,X)
#'
#'
#'
#'
endo.test <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01,method="OLS",invalid=TRUE){
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
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  stopifnot(is.character(method))

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
                         SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,SigmaYD=inputs$SigmaYD,tuning=tuning)
    Set = SetHats$VHat
  } else {
    SetHats <- endo.SHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                         SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,SigmaYD=inputs$SigmaYD,tuning=tuning)
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
  cat("Test result : our test statistics Q=",Q,"\n")
  if (abs(Q)>qnorm(1-alpha/2)) {
    cat("'H0 : Sigma12 = 0' is rejected","\n")
  } else {
    cat("'H0 : Sigma12 = 0' is not rejected","\n")
  }


  return(list(VHat=Set,Q=Q,Sigma12=Sigma12))
}

#' @title Relevant Instrumental Variable Selection for endogeneity test
#'
#' @description Implementation of Hard Thresholding for relevant instrumental variable selection. This function only takes the necessary inputs of the functions \code{\link{TSHT.OLS}} or \code{\link{TSHT.DeLasso}}. Any other methods to compute the necessary inputs can be adopted by the users according to their preferences.
#'
#' @param ITT_Y a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.
#' @param ITT_D a p_z by 1 numeric vector denoting the estimated coefficients of instruments in the treatment model.
#' @param WUMat a numeric matrix denoting WU where U is the precision matrix of W and W is the instrument-covariate matrix (Z, X).
#' @param SigmaSqY a numeric scalar denoting the consistent estimator of the noise level in the outcome model.
#' @param SigmaSqD a numeric scalar denoting the consistent estimator of the noise level in the treatment model.
#' @param SigmaYD a numeric scalar denoting the consistent estimator of the covariance between the error term in the treatment model and the error term in the outcome model.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater 2, with default 2.01.
#' @param bootstrap a logical value, default by FALSE.
#'
#' @return
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#' @export
#'
endo.SHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,bootstrap = FALSE,tuning = 2.01) {
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
  return(SHat)
}
