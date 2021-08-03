#' @title Causal inference in probit outcome models with possibly invalid IVs
#' @description Causal inference in probit outcome model with IVs, which provides the robust inference of the treatment effect in probit outcome models.
#'
#' @param Y binary and non-missing, n by 1 numeric outcome vector.
#' @param D continuous and non-missing, n by 1 numeric treatment vector.
#' @param Z continuous or discrete, non-missing, n by p_z numeric instrument matrix, containing p_z instruments.
#' @param X optional,continuous or discrete, n by p_x numeric covariate matrix, containing p_z covariates.
#' @param intercept a boolean scalar indicating to include the intercept or not. Default is TRUE.
#' @param method 'valid' or 'majority', indicating whether imposing the assumption of valid IVs (True) or possibly invalid IVs with majority rule (False). Default is True.
#' @param d1 a scalar for computing CATE(d1,d2|z0).
#' @param d2 a scalar for computing CATE(d1,d2|z0).
#' @param w0  a (pz+px) by 1 vector for computing CATE(d1,d2|z0).
#' @param bs.Niter a positive integer indicating the number of bootstrap resampling for computing the confidence interval.
#'
#' @return
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#'     \item{\code{betaHat}}{a numeric scalar denoting the estimate of beta/sig_u.}
#'     \item{\code{beta.sdHat}}{a numeric scalar denoting the estimated standard deviation of betaHat.}
#'     \item{\code{cateHat}}{a numeric scalar denoting the estimate of CATE(d1,d2|w0).}
#'     \item{\code{cate.sdHat}}{a numeric scalar denoting the estimated standard deviation of cateHat.}
#'     \item{\code{Maj.pass}}{True or False indicating whether the majority rule test is passed or not.}
#' @export
#'
#'
#'
#' @examples
#' \dontrun{
#' ### Working Low Dimensional Example ###
#' library(mvtnorm)
#' library(MASS)
#' n = 500; J = 7; s = 5; d1=-1; d2=1; z0=c(rep(0, J-1),0.1)
#' Z <- matrix(rnorm(n * J, 0, 1) , ncol = J, nrow = n)
#' gam <- c(rep(0.8, floor(J / 2)), rep(-0.8, J - floor(J / 2)))
#' cov.noise<-matrix(c(1,0.25, 0.25, 1),ncol=2)
#' noise.vec<-rmvnorm(n, rep(0,2), cov.noise)
#' v.vec<-noise.vec[,1]
#' X<-matrix(runif(n*2), ncol=2)
#' D = 0.5+Z %*% gam + X%*% rep(0.2,2) + v.vec
#'
#' pi0 <- c(rep(0, s), 0.4, 0.2)
#' beta0 <- 0.25
#' u.vec<- noise.vec[,2]
#' Y = (-0.5 + Z %*% pi0 + D * beta0 + u.vec>=0)
#' u1.r<-rnorm(2000,0,sd=1)
#' cace0 <- mean((as.numeric(-0.5+d1 * beta0 + z0 %*% pi0)+ u1.r )>=0) - mean((as.numeric(-0.5+d2 * beta0 + z0 %*% pi0) + u1.r)>=0)
#' ProbitControl(Y=Y, D=D, Z=Z, bs.Niter = 40, d1 = d1, d2 = d2, w0 = z0, method='valid', intercept=F)
#' ProbitControl(Y=Y, D=D, Z=Z, bs.Niter = 40, d1 = d1, d2 = d2, w0 = z0, method='majority', intercept=F)
#'}
#'
#'
#'
#'




ProbitControl<- function(Y, D, Z, X=NULL, d1=NULL, d2=NULL , w0=NULL, bs.Niter=40, intercept=T, method='majority'){
  stopifnot(method=='majority' | method=='valid')
  stopifnot(length(table(Y))==2)
  pz<- ncol(Z)
  px<-0
  if(!is.null(X)){
    Z<-cbind(Z,X)
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

  if(method=='majority'){
    #reduced form
    Gam.re<-glm(Y~cbind(Z,v.hat)-1, family=binomial(link='probit'))
    Gam.hat<-Gam.re$coef[-length(Gam.re$coef)]
    lam.hat<-Gam.re$coef[length(Gam.re$coef)]
    Gam.cov<-n*vcov(Gam.re)[1:pz,1:pz]+lam.hat^2*gam.cov
    Cov.gGam<-rbind(cbind(gam.cov,lam.hat^2*gam.cov),
                    cbind(lam.hat^2*Gam.cov, Gam.cov))
    #applying the majority rule
    Select.re<-Majority.test(n=n,ITT_Y=Gam.hat[1:pz], ITT_D=gam.hat[1:pz], Cov.gGam=Cov.gGam)
    SHat<-Select.re$SHat
    if(length(Select.re$VHat)<= length(SHat)/2){
      cat('Majority rule fails.','\n')
      Maj.pass=F
    }
    beta.hat<-median(Gam.hat[SHat]/gam.hat[SHat])
    pi.hat<- Gam.hat - gam.hat * beta.hat
    kappa.hat<-lam.hat-beta.hat #coef of v.hat
  }else{ #valid IV method
    Maj.pass=NA
    SHat=1:pz
    if((!intercept) & px==0){
      coef.re<-glm(Y~D+v.hat-1)$coef
      beta.hat<- coef.re[1]
      pi.hat<-c(rep(0,pz))
      kappa.hat<-coef.re[length(coef.re)]
    }else{
      if(intercept& px>0){
        X<-cbind(X,1)
      }else if(intercept & px==0){
        X<-matrix(rep(1,n),ncol=1)
      }
      coef.re<-glm(Y~D+X+v.hat-1)$coef
      beta.hat<- coef.re[1]
      pi.hat<-c(rep(0,pz),coef.re[2:(1+ncol(X))])
      kappa.hat<-coef.re[length(coef.re)]
    }
  }

  cace.hat<-NA; cace.sd<-NA
  if(!is.null(d1) & !is.null(d2) & !is.null(w0)){ #compute cate and its sd
    if(intercept){ w0=c(w0,1)}
    cace.hat = mean(pnorm(as.numeric(d1*beta.hat+w0%*%pi.hat) + v.hat*kappa.hat))-
      mean(pnorm(as.numeric(d2*beta.hat+w0%*%pi.hat) + v.hat*kappa.hat))
    #bootstrap for sd
    bs.lst<-list()
    for(i in 1:bs.Niter){
        bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
        bs.lst[[i]]<-Probit.boot.fun(data=bootstrap_data, pz=pz, d1=d1,d2=d2, w0=w0, SHat=SHat, method=method, intercept=intercept)
    }
    cace.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[1]))-cace.hat)^2))
    beta.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[2]))-beta.hat)^2))
  }else{ #only compute the sd of beta.hat
    bs.lst<-list()
    for(i in 1:bs.Niter){
      bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
      bs.lst[[i]]<-Probit.boot.fun(data=bootstrap_data, pz=pz, SHat=SHat, method=method, intercept=intercept)
      beta.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[2]))-beta.hat)^2))
    }
  }
  return(list(caceHat=cace.hat, cace.sdHat= cace.sd, betaHat=beta.hat, beta.sdHat=beta.sd, Maj.pass=Maj.pass))
}

Probit.boot.fun<-function(data, pz,d1=NULL, d2=NULL,w0=NULL, SHat, method=method, intercept=intercept){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  gam.bs<-lm(D~Z-1)$coef
  v.bs <- D- Z%*%gam.bs
  if(method=='majority'){
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
  if(!is.null(d1) & !is.null(d2) & !is.null(w0)){
    cace.bs = mean(pnorm(as.numeric(d1*beta.bs+w0%*%pi.bs) + v.bs*kappa.bs))-
      mean(pnorm(as.numeric(d2*beta.bs+w0%*%pi.bs) + v.bs*kappa.bs))
  }
  c(cace.bs, beta.bs)
}

Majority.test <- function(n, ITT_Y,ITT_D, Cov.gGam, tuning = 2.01, majority=T) {
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
    Sig1.j <- Var.comp.est(Cov.gGam[1:pz,1:pz], ITT_D, j)
    Sig2.j <- Var.comp.est(Cov.gGam[(pz+1):(2*pz),(pz+1):(2*pz)], ITT_D, j)
    Sig3.j <-  Var.comp.est(Cov.gGam[1:pz,(pz+1):(2*pz)], ITT_D, j)
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

