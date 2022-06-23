#' @title Two-Stage Hard Thresholding
#' @description Perform Two-Stage Hard Thresholding method, which provides the robust inference of the treatment effect in the presence of invalid instrumental variables in both low-dimensional and high-dimensional settings.
#'
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param alpha The significance level for the confidence interval. (default = 0.05)
#' @param method The method used to estimate the reduced form parameters. \code{"OLS"} stands for ordinary least squares, \code{"DeLasso"} stands for the debiased Lasso estimator, and \code{"Fast.DeLasso"} stands for the debiased Lasso estimator with fast algorithm. (default = \code{"OLS"})
#' @param voting The voting option used to estimate valid IVs. 'MP' stands for majority and plurality voting, 'MaxClique' stands for finding maximal clique in the IV voting matrix, and 'Conservative' stands for conservative voting procedure. Conservative voting is used to get an initial estimator of valid IVs in the Searching-Sampling method. (default= 'MaxClique').
#' @param robust If \code{TRUE}, the method is robust to heteroskedasticity errors. If \code{FALSE}, the method assumes homoskedasticity errors. When \code{robust = TRUE}, only \code{’OLS’} can be input to \code{method}. (default = \code{FALSE})
#'
#' @return
#'
#'     \code{TSHT} returns an object of class "TSHT".
#'     An object class "TSHT" is a list containing the following components:
#'     \item{\code{betaHat}}{The estimate of treatment effect.}
#'     \item{\code{beta.sdHat}}{The estimated standard error of \code{betaHat}.}
#'     \item{\code{ci}}{The 1-alpha confidence interval for \code{beta}.}
#'     \item{\code{SHat}}{The set of relevant IVs.}
#'     \item{\code{VHat}}{The set of relevant and valid IVs.}
#'     \item{\code{voting.mat}}{The voting matrix on whether the elements of each \code{SHat} are valid or not.}
#'     \item{\code{check}}{Whether the majority rule test is passed or not.}
#'     \item{\code{beta.clique}}{The estimates of treatment effect from each maximum clique. Only returns when \code{voting} is \code{'MaxClique'}.}
#'     \item{\code{beta.sd.clique}}{The estimated standard deviation of \code{betaHat} of each maximum clique. Only returns when \code{voting} is \code{'MaxClique'}}
#'     \item{\code{CI.clique}}{The 1-alpha confidence interval for \code{beta} of each maximum clique. Only returns when \code{voting} is \code{'MaxClique'}}
#'     \item{\code{max.clique}}{The maximum cliques voted as valid IVs. Only returns when \code{voting} is \code{'MaxClique'}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Y <- mroz[,"lwage"]
#' D <- mroz[,"educ"]
#' Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc","exper","expersq")])
#' X <- mroz[,"age"]
#' TSHT.model <- TSHT(Y=Y,D=D,Z=Z,X=X)
#' summary(TSHT.model)
#' }
#' @references {
#' Guo, Z., Kang, H., Tony Cai, T. and Small, D.S. (2018), Confidence intervals for causal effects with invalid instruments by using two-stage hard thresholding with voting, \emph{J. R. Stat. Soc. B}, 80: 793-815. \cr
#' }
TSHT <- function(Y,D,Z,X,intercept=TRUE, alpha=0.05,
                 method=c("OLS","DeLasso","Fast.DeLasso"), voting = 'MaxClique', robust = FALSE) {
  stopifnot(is.logical(robust))
  method = match.arg(method)
  if(method %in% c("DeLasso", "Fast.DeLasso") && robust==TRUE){
    robust = FALSE
    cat(sprintf("For methods %s, robust is set FALSE,
                as we only consider homoscedastic noise.\n", method))
  }
  if (robust == TRUE) {
    TSHTObject <- TSHT_hetero(Y = Y,D = D,Z = Z,X = X,intercept = intercept, alpha = alpha,
                                method = method, voting = voting)
  } else{
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

    # Reduced form estimator
    ITT_Y = inputs$ITT_Y;
    ITT_D = inputs$ITT_D;
    WUMat = inputs$WUMat;
    SigmaSqD = inputs$SigmaSqD;
    SigmaSqY = inputs$SigmaSqY;
    SigmaYD=inputs$SigmaYD;


    # Estimate Valid IVs
    SetHats = TSHT.VHat(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,
                        SigmaSqD = SigmaSqD,SigmaSqY = SigmaSqY,SigmaYD=SigmaYD,
                        voting = voting)
    VHat = SetHats$VHat; SHat = SetHats$SHat
    check = T
    if(length(VHat)< length(SHat)/2){
      cat('Majority rule fails.','\n')
      check=F
    }

    if (voting == 'MaxClique') {
      max.clique <- SetHats$max.clique
      max.clique.mat <- matrix(0,nrow = length(max.clique),ncol = length(max.clique[[1]]))
      CI.temp <- matrix(0,nrow = length(max.clique), ncol = 2)
      beta.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
      betavar.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
      for (i in 1:length(max.clique)) {
        temp <- SHat[sort(as.numeric(max.clique[[i]]))]
        max.clique.mat[i,] <- temp
        if (!is.null(colnames(Z))) {
          max.clique.mat[i,] <- colnames(Z)[max.clique.mat[i,]]
        }
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
    if (!is.null(colnames(Z))) {
      SHat = colnames(Z)[SHat]
      VHat = colnames(Z)[VHat]
    }
    if (voting != 'MaxClique') {
      ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat / n),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat/n))
      TSHTObject <- list( betaHat=betaHat,beta.sdHat = sqrt(betaVarHat/n),ci=ci,SHat=SHat,VHat = VHat,voting.mat=SetHats$voting.mat,check = check)
    } else {
      TSHTObject <- list( betaHat=betaHat,beta.sdHat = sqrt(betaVarHat/n),ci=CI.union,SHat=SHat,VHat = VHat,voting.mat=SetHats$voting.mat, check = check,
                  beta.clique = beta.temp,beta.sd.clique = sqrt(betavar.temp/n), CI.clique = CI.temp,
                  max.clique = max.clique.mat)
      }
  }
  class(TSHTObject) <- 'TSHT'

  TSHTObject
}



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

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}



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

TSHT.SIHR <- function(Y, D, W, pz, method="OLS", intercept=TRUE){
  n = nrow(W)
  covW = t(W)%*%W / n
  init_Y = Lasso(W, Y, lambda="CV.min", intercept=intercept)
  init_D = Lasso(W, D, lambda="CV.min", intercept=intercept)

  if(intercept) W_int = cbind(1, W) else W_int = W
  resid_Y = as.vector(Y - W_int%*%init_Y)
  resid_D = as.vector(D - W_int%*%init_D)

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
  # Temp = t(WUMat)%*%WUMat / n
  # V.Gamma = SigmaSqY * Temp
  # V.gamma = SigmaSqD * Temp
  # C = SigmaYD * Temp

  return(list(ITT_Y = ITT_Y,
              ITT_D = ITT_D,
              WUMat = WUMat,
              SigmaSqY = SigmaSqY,
              SigmaSqD = SigmaSqD,
              SigmaYD = SigmaYD))

}



TSHT.VHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,
                      voting = 'MaxClique') {
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
  stopifnot(voting=='MP' | voting=='MaxClique' | voting == 'Conservative')

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
      sqrt(log(n)/n)
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


#' Summary of TSHT
#'
#' @param object TSHT object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.TSHT <- function(object,...){
  return(object)
}

#' print of TSHT
#'
#' @param x TSHT object
#' @param ...
#' @keywords internal
#' @return
#' @export
print.TSHT <- function(x,...){
  TSHT <- x
  cat("\nRelevant Instruments:", TSHT$SHat, "\n");
  cat("\nValid Instruments:", TSHT$VHat,"\n","\nThus, Majority rule",ifelse(TSHT$check,"holds.","does not hold."), "\n");
  cat(rep("_", 30), "\n")
  cat("\nBetaHat:",TSHT$betaHat,"\n");
  cat("\nConfidence interval for Beta: [", TSHT$ci[1], ",", TSHT$ci[2], "]", "\n", sep = '');
  if (!is.null(TSHT$beta.clique)) {
    cat(rep("_", 30), "\n")
    cat("\nResults from each maximum clique\n")
    for (i in 1:nrow(TSHT$max.clique)) {
      cat("\nMaximum clique", i,":", TSHT$max.clique[i,], "\n");
      cat("\nBetaHat from the clique:",TSHT$beta.clique[i,],"\n");
      cat("\nConfidence interval for Beta from the clique: [", TSHT$CI.clique[i,1], ",", TSHT$CI.clique[i, 2], "]", "\n", sep = '');
      cat(rep("_", 30), "\n")
    }
  }

}
