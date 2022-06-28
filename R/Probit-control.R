#' @title Causal inference in probit outcome models with possibly invalid IVs
#' @description Perform causal inference in the probit outcome model with possibly invalid IVs under the majority rule.
#'
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param invalid If \code{TRUE}, the method is robust to the presence of possibly invalid IVs; If \code{FALSE}, the method assumes all IVs to be valid. (default = \code{FALSE})
#' @param d1 A treatment value for computing CATE(d1,d2|w0).
#' @param d2 A treatment value for computing CATE(d1,d2|w0).
#' @param w0  A vector for computing CATE(d1,d2|w0).
#' @param bs.Niter The number of bootstrap resampling for computing the confidence interval.
#'
#' @return
#'     \code{ProbitControl} returns an object of class "SpotIV", which is a list containing the following components:
#'     \item{\code{betaHat}}{The estimate of the model parameter in front of the treatment.}
#'     \item{\code{beta.sdHat}}{The estimated standard error of betaHat.}
#'     \item{\code{cateHat}}{The estimate of CATE(d1,d2|w0).}
#'     \item{\code{cate.sdHat}}{The estimated standard deviation of \code{cateHat}.}
#'     \item{\code{SHat}}{The estimated set of relevant IVs.}
#'     \item{\code{VHat}}{The estimated set of relevant and valid IVs.}
#'     \item{\code{Maj.pass}}{Indicator for whether the majority rule is satisfied or not.}
#'
#' @importFrom stats binomial glm median pnorm sd
#' @export
#'
#'
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
#' Probit.model <- ProbitControl(Y0,D,Z,X,d1 = d1,d2 = d2,w0 = w0)
#' summary(Probit.model)
#'}
#'
#' @references {
#' Li, S., Guo, Z. (2020), Causal Inference for Nonlinear Outcome Models with Possibly Invalid Instrumental Variables, Preprint \emph{arXiv:2010.09922}.\cr
#' }
#'
#'




ProbitControl<- function(Y, D, Z, X=NULL, intercept=TRUE, invalid=FALSE,
                         d1=NULL, d2=NULL , w0=NULL, bs.Niter=40){
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),is.vector(Y)||(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  if (is.vector(Y)) {
    Y <- cbind(Y)
  }
  Y = as.numeric(Y)
  stopifnot(length(table(Y))==2)
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

  # All the other argument
  stopifnot(is.logical(invalid))
  stopifnot(is.logical(intercept))

  pz<- ncol(Z)
  px<-0
  if(!is.null(X)){
    Z<-cbind(Z,X)
    X <- cbind(X)
    px<-ncol(X)
  }
  stopifnot(length(w0)==ncol(Z))

  if(intercept){ Z<-cbind(Z,1) }
  n <- length(Y)
  Maj.pass=T

  #first stage
  gam.re<-lm(D ~ Z -1)
  gam.hat<-gam.re$coef
  v.hat<-D-Z%*%gam.hat
  sig.v.hat<- mean(v.hat^2)
  gam.cov<-n*vcov(gam.re)[1:pz,1:pz]

  if(invalid){
    #reduced form
    Gam.re<-glm(Y~cbind(Z,v.hat)-1, family=binomial(link='probit'))
    Gam.hat<-Gam.re$coef[-length(Gam.re$coef)]
    lam.hat<-Gam.re$coef[length(Gam.re$coef)]
    Gam.cov<-as.matrix(n*vcov(Gam.re)[1:pz,1:pz]+lam.hat^2*gam.cov)
    Cov.gGam<-rbind(cbind(gam.cov,lam.hat^2*gam.cov),
                    cbind(lam.hat^2*Gam.cov, Gam.cov))
    #applying the majority rule
    Select.re<-Majority.test(n=n,ITT_Y=Gam.hat[1:pz], ITT_D=gam.hat[1:pz], Cov.gGam=Cov.gGam)
    SHat<-Select.re$SHat
    VHat <- Select.re$VHat
    if(length(Select.re$VHat)<= length(SHat)/2){
      cat('Majority rule fails.','\n')
      Maj.pass=F
    }
    beta.hat<-median(Gam.hat[SHat]/gam.hat[SHat])
    pi.hat<- Gam.hat - gam.hat * beta.hat
    kappa.hat<-lam.hat-beta.hat #coef of v.hat
  }else{ #valid IV method
    SHat=1:pz
    if((!intercept) && px==0){
      coef.re<-glm(Y~D+v.hat-1)$coef
      beta.hat<- coef.re[1]
      pi.hat<-c(rep(0,pz))
      kappa.hat<-coef.re[length(coef.re)]
    }else{
      if(intercept && (px>0)){
        X<-cbind(X,1)
      }else if(intercept && px==0){
        X<-matrix(rep(1,n),ncol=1)
      }
      coef.re<-glm(Y~D+X+v.hat-1)$coef
      beta.hat<- coef.re[1]
      pi.hat<-c(rep(0,pz),coef.re[2:(1+ncol(X))])
      kappa.hat<-coef.re[length(coef.re)]
    }
    VHat = SHat
  }

  cace.hat<-NA; cace.sd<-NA
  if(!is.null(d1) && !is.null(d2) && !is.null(w0)){ #compute cate and its sd
    if(intercept){ w0=c(w0,1)}
    cace.hat = mean(pnorm(as.numeric(d1*beta.hat+w0%*%pi.hat) + v.hat*kappa.hat))-
      mean(pnorm(as.numeric(d2*beta.hat+w0%*%pi.hat) + v.hat*kappa.hat))
    #bootstrap for sd
    bs.lst<-list()
    for(i in 1:bs.Niter){
      sample.true <- F
      while(sample.true==F)
        {tryCatch({
            bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),];
            bs.lst[[i]]<-Probit.boot.fun(data=bootstrap_data, pz=pz, d1=d1,d2=d2,
                            w0=w0, SHat=SHat, invalid=invalid, intercept=intercept);
            sample.true<-T;
            },error=function(e){
            },finally={})
        }
      }
    cace.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[1]))-cace.hat)^2))
    beta.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[2]))-beta.hat)^2))
  }else{ #only compute the sd of beta.hat
    bs.lst<-list()
    for(i in 1:bs.Niter){
      sample.true <- F
      while(sample.true==F)
      {tryCatch({
        bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),];
        bs.lst[[i]]<-Probit.boot.fun(data=bootstrap_data, pz=pz, d1=d1,d2=d2,
                                     w0=w0, SHat=SHat, invalid=invalid, intercept=intercept);
        sample.true<-T;
      },error=function(e){
      },finally={})
      }
      beta.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[2]))-beta.hat)^2))
    }
  }
  VHat <- as.numeric(VHat)
  if (!is.null(colnames(Z))) {
    SHat = colnames(Z)[SHat]
    VHat = colnames(Z)[VHat]
  }
  Probit.model <- list(betaHat=beta.hat, beta.sdHat=beta.sd, cateHat=cace.hat, cate.sdHat= cace.sd, SHat=SHat, VHat = VHat, Maj.pass=Maj.pass)
  class(Probit.model) <- "SpotIV"
  return(Probit.model)
}

Probit.boot.fun<-function(data, pz,d1=NULL, d2=NULL,w0=NULL, SHat, invalid=invalid, intercept=intercept){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  gam.bs<-lm(D~Z-1)$coef
  v.bs <- D- Z%*%gam.bs
  if(invalid){
    Gam.bs.re<-glm(Y~cbind(Z,v.bs)-1, family=binomial(link='probit'))
    Gam.bs<-Gam.bs.re$coef[-length(Gam.bs.re$coef)]
    lam.bs<-Gam.bs.re$coef[length(Gam.bs.re$coef)]
    beta.bs<-median(Gam.bs[SHat]/gam.bs[SHat])
    pi.bs <- Gam.bs - gam.bs *beta.bs
    kappa.bs<-lam.bs-beta.bs
  }else{
    if(pz==ncol(Z)){
      coef.bs<-glm(Y~D+v.bs-1)$coef
      beta.bs<- coef.bs[1]
      pi.bs<-c(rep(0,pz))
      kappa.bs<-coef.bs[length(coef.bs)]
    }else{
      X<- as.matrix(Z[,-(1:pz)])# the last column is the intercept
      coef.bs<-glm(Y~D+X+v.bs-1)$coef
      beta.bs<- coef.bs[1]
      pi.bs<-c(rep(0,pz),coef.bs[2:(1+ncol(X))])
      kappa.bs<-coef.bs[length(coef.bs)]
    }

  }

  cace.bs<-NA
  if(!is.null(d1) && !is.null(d2) && !is.null(w0)){
    cace.bs = mean(pnorm(as.numeric(d1*beta.bs+w0%*%pi.bs) + v.bs*kappa.bs))-
      mean(pnorm(as.numeric(d2*beta.bs+w0%*%pi.bs) + v.bs*kappa.bs))
  }
  c(cace.bs, beta.bs)
}

Majority.test <- function(n, ITT_Y,ITT_D, Cov.gGam, tuning = 2.01) {
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

