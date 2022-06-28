#' @title SpotIV method for causal inference in semi-parametric outcome model
#' @description Perform causal inference in the semi-parametric outcome model with possibly invalid IVs.
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param invalid If TRUE, the method is robust to the presence of possibly invalid IVs; If FALSE, the method assumes all IVs to be valid. (default = \code{FALSE})
#' @param d1 A treatment value for computing CATE(d1,d2|w0).
#' @param d2 A treatment value for computing CATE(d1,d2|w0).
#' @param w0  A value of measured covariates and instruments for computing CATE(d1,d2|w0).
#' @param M.est If \code{TRUE}, estimate \code{M} based on BIC, otherwise \code{M} is not estimated and input or default value of \code{M} is used. (default = \code{TRUE})
#' @param M The dimension of indices in the outcome model, from 1 to 3. (default = \code{2})
#' @param bs.Niter The number of bootstrap resampling for computing the confidence interval. (default = \code{40})
#' @param bw  A (M+1) by 1 vector bandwidth specification. (default = \code{NULL})

#' @return
#'     \code{SpotIV} returns an object of class "SpotIV", which "SpotIV" is a list containing the following components:
#'     \item{\code{betaHat}}{The estimate of the model parameter in front of the treatment.}
#'     \item{\code{cateHat}}{The estimate of CATE(d1,d2|w0).}
#'     \item{\code{cate.sdHat}}{The estimated standard error of cateHat.}
#'     \item{\code{SHat}}{The set of relevant IVs.}
#'     \item{\code{VHat}}{The set of relevant and valid IVs.}
#'     \item{\code{Maj.pass}}{Indicator for whether the majority rule is satisfied or not.}
#' @import dr
#' @import orthoDr
#' @import foreach
#' @importFrom stats binomial glm median pnorm sd
#' @export
#'

#' @examples
#' \dontrun{
#' Y <- mroz[,"lwage"]
#' D <- mroz[,"educ"]
#' Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc","exper","expersq")])
#' X <- mroz[,"age"]
#' Y0 <- as.numeric((Y>median(Y)))
#' d2 = median(D); d1 = d2+1;
#' w0 = apply(cbind(Z,X)[which(D == d2),], 2, mean)
#' SpotIV.model <- SpotIV(Y0,D,Z[,-5],X,d1 = d1,d2 = d2,w0 = w0[-5])
#' summary(SpotIV.model)
#'}
#'
#' @references {
#' Li, S., Guo, Z. (2020), Causal Inference for Nonlinear Outcome Models with Possibly Invalid Instrumental Variables, Preprint \emph{arXiv:2010.09922}.\cr
#' }
#'
#'
#'

SpotIV<- function(Y, D, Z, X=NULL, intercept=TRUE, invalid=FALSE,  d1, d2 , w0,
                  M.est=TRUE, M=2, bs.Niter=40, bw=NULL){
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
  }

  # Check d1,d2, and w0
  stopifnot(!missing(d1),!missing(d2),!missing(w0),is.numeric(d1),is.numeric(d2),is.numeric(w0))

  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.logical(invalid))
  stopifnot(is.logical(M.est))

  pz<- ncol(Z)
  px<-0
  if(!is.null(X)){
    Z<-cbind(Z,X)
    X <- cbind(X)
    px<-ncol(X)
  }
  stopifnot(length(w0)==ncol(Z))
  n <- length(Y)
  Maj.pass=T

  if (invalid) {
    V <- NULL
  } else {
    V <- 1:pz
  }

  #first-stage regression
  if(intercept){
    gam.re<- lm(D ~ Z)
    gam.hat<-gam.re$coef[-1]
    var.gam <- n * as.matrix(vcov(gam.re)[2:(pz+1),2:(pz+1)])
  }else{
    gam.re<- lm(D ~ Z-1)
    gam.hat<-gam.re$coef
    var.gam <- n * vcov(gam.re)[1:pz,1:pz]
  }

  v.hat <- D-Z%*%gam.hat

  #voting and applying the majority rule
  if(is.null(V)){
    #get reduced-from
    SIR.re <-SIR.est(X.cov=cbind(Z,v.hat), Y, M= M, M.est=M.est)
    M <- ncol(SIR.re$theta.hat)
    Gam.hat<-as.matrix(SIR.re$theta.hat[1:ncol(Z),], ncol=M) #estimate Gamma
    ##voting
    Select.re<-Majority.test(n=n,ITT_Y=Gam.hat[1:pz,1], ITT_D=gam.hat[1:pz],
                         Cov.gGam=rbind(cbind(var.gam, diag(0,pz)),
                                            cbind(diag(0,pz), SIR.re$vGam[1:pz,1:pz])))
    SHat<-Select.re$SHat
    if(length(Select.re$VHat)< length(SHat)/2){
      cat('Majority rule fails.','\n')
      Maj.pass=F
    }
    VHat <- as.numeric(Select.re$VHat)
    beta.hat<-sapply(1:M, function(m) median(Gam.hat[SHat,m]/gam.hat[SHat]))
    beta.hat <- matrix(beta.hat,nrow=1,ncol=M)
    pi.hat<- Gam.hat - gam.hat %*% beta.hat
  }else{###oracle method
    SIR.re <-SIR.est(X.cov=cbind(Z%*%gam.hat,Z[,-V],v.hat), Y, M= M, M.est=M.est)
    M <- ncol(SIR.re$theta.hat)
    beta.hat<-matrix(SIR.re$theta.hat[1,],nrow=1,ncol=M)
    pi.hat<-matrix(0,nrow=ncol(Z), ncol=M)
    pi.hat[V,]<-0
    pi.hat[-V,]<-SIR.re$theta.hat[2:(ncol(Z)-length(V)+1),]
    SHat <- V
    VHat <- V
  }
  ##estimate cace##
  if(is.null(bw)){
    bw.z<- apply(cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat), 2,
                 function(x) 0.9*n^(-1/6)*min(sd(x), (quantile(x,0.75)-quantile(x,0.25))/1.34))
  }

  asf.dw<-Spot.ASF.est(d1=d1,d2=d2, z0=w0, beta.hat= beta.hat, pi.hat= pi.hat,
                  Y, D, Z, v.hat=v.hat, bw.z=bw.z)
  cace.hat = asf.dw$cace.hat #cace
  bw.z=asf.dw$bw.z
  boot_b<-list()
  for(i in 1: bs.Niter){
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    boot_b[[i]]<-Spot.boot.fun(data=bootstrap_data, M=M, d1=d1,d2=d2, w0=w0, SHat=SHat,
                               bw.z=bw.z, V=V, intercept=intercept, pz=pz)
  }
  cace.sd<-sqrt(mean((unlist(lapply(boot_b, function(x) x[1]))-cace.hat)^2))
  VHat <- as.numeric(VHat)
  if (!is.null(colnames(Z))) {
    SHat = colnames(Z)[SHat]
    VHat = colnames(Z)[VHat]
  }
  SpotIV.model <- list(betaHat = beta.hat, cateHat=cace.hat, cate.sdHat= cace.sd,
                       SHat=SHat, VHat = VHat, Maj.pass=Maj.pass)
  class(SpotIV.model) <- "SpotIV"
  return(SpotIV.model)
}
#' Summary of SpotIV
#'
#' @param object SpotIV object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.SpotIV <- function(object,...){
  SpotIV <- object
  cat("\nRelevant Instruments:", SpotIV$SHat, "\n");
  cat("\nValid Instruments:", SpotIV$VHat,"\n","\nThus, Majority rule",ifelse(SpotIV$Maj.pass,"holds.","does not hold."), "\n");
  cat(rep("_", 30), "\n")
  cat("\nBetaHat:",SpotIV$betaHat,"\n");
  if (!is.null(SpotIV$beta.sdHat)) {
    cat("\nStandard error of BetaHat:",SpotIV$beta.sdHat,"\n");
    ci.beta <- c(SpotIV$betaHat-qnorm(0.975)*SpotIV$beta.sdHat,SpotIV$betaHat+qnorm(0.975)*SpotIV$beta.sdHat)
    cat("\n95% Confidence Interval for Beta: [", ci.beta[1], ",", ci.beta[2], "]", "\n", sep = '');
  }
  cat("\nCATEHat:",SpotIV$cateHat,"\n");
  cat("\nStandard error of CATEHat:",SpotIV$cate.sdHat,"\n");
  ci <- c(SpotIV$cateHat-qnorm(0.975)*SpotIV$cate.sdHat,SpotIV$cateHat+qnorm(0.975)*SpotIV$cate.sdHat)
  cat("\n95% Confidence Interval for CATE: [", ci[1], ",", ci[2], "]", "\n", sep = '');
}


SIR.est<- function(X.cov,Y, M=2, M.est=TRUE){
  p<- ncol(X.cov)
  n<-length(Y)
  if(M.est){
    nslice=ifelse(length(table(Y))==2,2,8)
    SIR.re<-dr::dr(Y ~ X.cov -1, method='sir', numdir=M, nslices=nslice)
    evalues=SIR.re$evalues
    nobs.slice <- median(SIR.re$slice.info$slice.sizes)
    M <- which.max(
      sapply(1:M, function(m) sum(log(evalues[(m+1):p]+1)-evalues[(m+1):p])*n/2-log(n)*m*(2*p-m+1)/4/nobs.slice))
  }
  if(length(table(Y))==2){
    init.re<- glm(Y~X.cov-1,family=binomial(link='logit'))
  }else{
    init.re<-lm(Y~X.cov-1)
  }
  Gam.init<-init.re$coef
  vGam<-vcov(init.re)*n
  if(M==1){
    theta.hat <- orthoDr::orthoDr_reg(x=X.cov, y = Y, B.initial =as.matrix(Gam.init,ncol=1), ndr=1)$B
  }else{
    theta.hat<-dr::dr(Y ~ X.cov -1, method='sir', numdir=M)$evectors[,1:M] #using a faster computation
  }
  list(theta.hat=theta.hat, vGam=vGam)
}

Majority.test <- function(n, ITT_Y,ITT_D, Cov.gGam, tuning = 2.01, majority=TRUE) {
  Var.comp.est <- function(Cov.mat, gam.hat, j){
    diag(Cov.mat) + (gam.hat/gam.hat[j])^2 * Cov.mat[j,j] - 2*gam.hat/gam.hat[j] * Cov.mat[j,]
  }
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)
  # Check Sigmas
  stopifnot(!missing(Cov.gGam))

  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)

  # Constants
  pz = length(ITT_Y)

  # First Stage
  Var.gam.hat <- diag(Cov.gGam)[1:pz]
  SHat = (1:pz)[abs(ITT_D) >= (sqrt(Var.gam.hat) * sqrt(tuning*log(pz)/n))]
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
    #compute three components in eq(33)
    Sig1.j <- Var.comp.est(as.matrix(Cov.gGam[1:pz,1:pz]), ITT_D, j)
    Sig2.j <- Var.comp.est(as.matrix(Cov.gGam[(pz+1):(2*pz),(pz+1):(2*pz)]), ITT_D, j)
    Sig3.j <-  Var.comp.est(as.matrix(Cov.gGam[1:pz,(pz+1):(2*pz)]), ITT_D, j)
    sigmasq.j <- beta.j^2 *Sig1.j +  Sig2.j - 2* beta.j * Sig3.j
    PHat.bool.j <- abs(pi.j) <= sqrt(sigmasq.j) * tuning * sqrt(log(pz)/n)
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }

  # Voting
  diag(VHats.bool) <- rep(TRUE, nCand)
  VM = rowSums(VHats.bool)
  # cat(VM,'\n')
  VHat = rownames(VHats.bool)[VM > (0.5 * length(SHat))] # Majority winners
  return(list(VHat = VHat,SHat=SHat))
}



Spot.boot.fun<-function(data, M, d1, d2,w0, SHat,
                        bw.z=NULL, V=NULL, intercept, pz){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  if(intercept){
    gam.re<- lm(D ~ Z)
  }else{
    gam.re<- lm(D ~ Z-1)
  }
  gam.bs<-gam.re$coef[1:ncol(Z)]
  v.bs <- D- Z%*%gam.bs
  if(is.null(V)){
    SIR.bs.re <- SIR.est(cbind(Z,v.bs), Y, M=M, M.est=F)
    Gam.bs<-as.matrix(SIR.bs.re$theta.hat[1:ncol(Z),], ncol=M)
    beta.bs<-sapply(1:M, function(m) median(Gam.bs[SHat,m]/gam.bs[SHat]))
    pi.bs <- Gam.bs - gam.bs %*% matrix(beta.bs,nrow=1,ncol=M)
  }else{
    SIR.bs <-SIR.est(X.cov=cbind(Z%*%gam.bs,Z[,-V],v.bs), Y, M= M, M.est=F)
    M <- ncol(SIR.bs$theta.hat)
    beta.bs<-SIR.bs$theta.hat[1,]
    pi.bs<-matrix(0,nrow=ncol(Z), ncol=M)
    pi.bs[V,]<-0
    pi.bs[-V,]<-SIR.bs$theta.hat[2:(pz-length(V)+1),]
  }
  asf.dw<-Spot.ASF.est(d1=d1,d2=d2, z0=w0, beta.hat= beta.bs, pi.hat= pi.bs,
                  Y=Y, D=D, Z=Z, v.hat=D-Z%*%gam.bs, bw.z=bw.z)
  asf.dw$cace.hat

}

box.ker<-function(xx, X, h){
  apply(X, 1, function(xi) max(abs(xi-xx)/h)<1)/prod(h)
}

Spot.ASF.est <- function(d1, d2, z0, beta.hat, pi.hat,  Y, D, Z, v.hat, bw.z=NULL){
  beta.hat=as.vector(beta.hat)
  M=length(beta.hat)
  n=length(D)
  v.hat <- as.matrix(v.hat,ncol=1)
  D<-as.matrix(D,ncol=1)
  index <- cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat)

  #cace
  index1.z <- cbind(apply(matrix(d1%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat)
  index2.z <- cbind(apply(matrix(d2%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat)
  #cat(index1.z[1,1],index2.z[1,1],'\n')
  q1<-quantile(index[,1], 0.975)
  q2<-quantile(index[,1], 0.025)

  if(sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0 & sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0){
    index1.z<- index1.z[index1.z[,1]<= q1 & index1.z[,1]>= q2,]
    index2.z<- index2.z[index2.z[,1]<= q1 & index2.z[,1]>= q2,]
  }
  asf.dw1<- NA
  asf.dw2<-NA
  bw.z<-bw.z/1.5
  while(is.na(asf.dw1) | is.na(asf.dw2)){
    bw.z<-bw.z*1.5
    if(M==1){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] ])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2]])), na.rm=T)
    }else if(M==2){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
    }else if(M==3){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3] & abs(xx[4]-index[,4])<=bw.z[4]])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3] & abs(xx[4]-index[,4])<=bw.z[4]])), na.rm=T)
    }
  }


  cace.hat <- asf.dw1-asf.dw2


  list(cace.hat = cace.hat, ace.hat=0, bw=bw.z, bw.z=bw.z)
}

