#' @title Searching-Sampling method
#' @description Construction of confidence intervals with Searching-Sampling method
#'
#' @param Y continuous and non-missing, n by 1 numeric outcome vector.
#' @param D continuous or discrete, non-missing, n by 1 numeric treatment vector.
#' @param Z continuous or discrete, non-missing, n by p_z numeric instrument matrix, containing p_z instruments.
#' @param X optional,continuous or discrete, n by p_x numeric covariate matrix, containing p_z covariates.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param Sampling a boolean scalar indicating to implement Sampling method, with default TRUE.
#' @param boot.value a boolean scalar indicating to implement bootstrap, with default TRUE.
#' @param M a positive integer indicating the number of bootstrap resampling for computing the confidence interval, with default 1000.
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix, with default FALSE.
#' @param alpha0 a numeric scalar value between 0 and 1 indicating the sampling threshold level for the generated samples, with default 0.01.
#'
#' @return
#' \item{\code{CI}}{a two dimensional numeric vector denoting the unioned 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#' \item{\code{CI.matrix}}{a numeric matrix denoting CI for appropriate candidates of beta. Only returns when \code{Sampling} is \code{FALSE}.}
#' \item{\code{rule}}{a boolean scalar denoting whether the identification condition is satisfied or not}
#' \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs.}
#' \item{\code{CI.clique}}{a numeric matrix where each row represents the CI corresponding to each maximum clique. Only returns when \code{max_clique} is \code{TRUE}.}
#' \item{\code{rule.clique}}{a boolean matrix where each row represents whether the identification condition is satisfied or not. Only returns when \code{max_clique} is \code{TRUE}.}
#' \item{\code{max.cliques}}{a numeric matrix where each row represents each maximum clique. Only returns when \code{max_clique} is \code{TRUE}.}
#' @export
#'
#' @examples
#' \dontrun{
#' ### example ###
#'
#' A1gen<-function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'    for(j in 1:p){
#'       A1[i,j]<-rho^(abs(i-j))
#'     }
#'   }
#'   A1
#' }
#'
#' ### Set the model parameter ###
#' n: sample size
#' IV.str: individual IV strength
#' VIO.str: violation strength
#' beta:   true treatment effect (set at 1)
#' px:number of covariates
#' L: number of candiate IVs
#' n = 500;
#' IV.str=0.5;
#' VIO.str=0.4;
#' pi.value<-IV.str*VIO.str
#' beta = 1;
#' px <- 10;
#' L = 10;
#' p=L+px
#' phi<-rep(0,px)
#' psi<-rep(0,px)
#' phi[1:px]<-(1/px)*seq(1,px)+0.5
#' psi[1:px]<-(1/px)*seq(1,px)+1
#' rho=0.5
#' Cov<-(A1gen(rho,p))
#'
#' ### Generate a setting where only the plurality rule holds
#' # 6 invalid IVs and 4 valid IVs
#' s1 = 3;
#' s2 = 3;
#' s=s1+s2
#' alpha = c(rep(0,L-s),rep(pi.value,s1),-seq(1,s2)/2);
#' gamma=rep(IV.str,L)
#'
#' epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#' W<-mvrnorm(n, rep(0, p), Cov)
#' Z=W[,1:L]
#' X=W[,(L+1):p]
#' epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
#' epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
#' D = 0.5 + Z %*% gamma+ X%*% psi + epsilon[,1]
#' Y = -0.5 + Z %*% alpha + D * beta + X%*%phi+ epsilon[,2]
#' Searching.Sampling(Y,D,Z,X)
#' Searching.Sampling(Y,D,Z,X,Sampling=FALSE)
#' }
Searching.Sampling <- function(Y, D, Z, X, intercept = TRUE, alpha = 0.05, alpha0 = 0.01,
                               Sampling=TRUE, boot.value=TRUE, M = 1000, max_clique=FALSE){
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
  stopifnot(is.logical(Sampling))

  n <- length(Y)
  pz<-ncol(Z)
  if (intercept) {
    W <- cbind(W,1)
  }

  inputs = TSHT.OLS(Y,D,W,pz,intercept=FALSE)
  # Compute covariance of W and W %*% U
  covW = t(W) %*% W /n  #this should automatically turn covW into a matrix
  WUMat = inputs$WUMat #W %*% (solve(covW))[,1:pz]
  # qrW = qr(W)
  ITT_Y = inputs$ITT_Y # qr.coef(qrW,Y)[1:pz]
  ITT_D = inputs$ITT_D #qr.coef(qrW,D)[1:pz]
  SigmaSqY = inputs$SigmaSqY #sum(qr.resid(qrW,Y)^2)/(n -p-1)
  SigmaSqD = inputs$SigmaSqD #sum(qr.resid(qrW,D)^2)/(n -p-1)
  SigmaYD = inputs$SigmaYD #sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / (n - p-1)



  ### screen out strongly invalid IVs and retain a set of valid and weakly invalid IVs ###
  TSHT.Init <- TSHT.Initial(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,covW,max_clique=max_clique,bootstrap = TRUE)

  if (max_clique) {

    ### if max_clique=TRUE, implement the method through all individual maximum cliques ###

    max.clique <- TSHT.Init$max.clique
    max.clique.mat <- matrix(0,nrow = length(max.clique),ncol = length(max.clique[[1]]))
    CI.clique <- matrix(0,nrow = length(max.clique), ncol = 2)
    rule.clique <- matrix(FALSE,nrow = length(max.clique),ncol = 1)
    for (i in 1:length(max.clique)) {
      V0.hat <- max.clique[[i]]
      temp<-(SigmaSqY)/(ITT_D[V0.hat]^2)+SigmaSqD*(ITT_Y[V0.hat]^2)/(ITT_D[V0.hat]^4)-
        2*SigmaYD*(ITT_Y[V0.hat])/(ITT_D[V0.hat]^3)
      var.beta<-(diag(solve(covW)/n)[1:pz])[V0.hat]*temp
      CI.initial<-matrix(NA,nrow=length(V0.hat),ncol=2)
      CI.initial[,1]<-(ITT_Y/ITT_D)[V0.hat]-sqrt(log(n)*var.beta)
      CI.initial[,2]<-(ITT_Y/ITT_D)[V0.hat]+sqrt(log(n)*var.beta)
      uni<- intervals::Intervals(CI.initial)
      CI.initial.union<-base::as.matrix(intervals::interval_union(uni))
      beta.grid.seq<-analysis.CI(CI.initial.union,grid.size=n^{-1})$grid.seq
      CI.sea<-Searching.CI(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD=SigmaSqD,SigmaSqY=SigmaSqY,
                           SigmaYD=SigmaYD,InitiSet = V0.hat,WUMat = WUMat,alpha = alpha,
                           beta.grid=beta.grid.seq,bootstrap=FALSE)
      CI.temp<-CI.sea$CI.search
      beta.grid<-analysis.CI(base::as.matrix(CI.sea$CI.search),beta,n^{-0.6})$grid.seq

      ### conduct the refined searching ###

      CI.sea.refined<-Searching.CI(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD=SigmaSqD,SigmaSqY=SigmaSqY,
                                   SigmaYD=SigmaYD,InitiSet = V0.hat,WUMat = WUMat,alpha = alpha,
                                   beta.grid,bootstrap=boot.value)
      CI.search<-CI.sea.refined$CI.search
      CI.clique[i,] <- t(base::as.matrix(c(min(CI.search),max(CI.search)))) # min and max of CIs
      max.clique.mat[i,] <- V0.hat
      rule.clique[i,] <- CI.sea.refined$rule
      if (Sampling) {
        CI.sampling<-Searching.CI.Sampling(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD=SigmaSqD,SigmaSqY=SigmaSqY,
                                           SigmaYD=SigmaYD,InitiSet = V0.hat,WUMat = WUMat,alpha = alpha,
                                           beta.grid,M=M,bootstrap=boot.value,alpha0=alpha0)
        CI.clique[i,]<-t(base::as.matrix(c(min(CI.sampling$CI.union[,1]),max(CI.sampling$CI.union[,2])))) # min and max of CIs
        rule.clique[i,] <- CI.sampling$rule
      }
    }
  }



  V0.hat<-TSHT.Init$VHat

  ### construct the initial range [L,U] ###

  temp<-(SigmaSqY)/(ITT_D[V0.hat]^2)+SigmaSqD*(ITT_Y[V0.hat]^2)/(ITT_D[V0.hat]^4)-
    2*SigmaYD*(ITT_Y[V0.hat])/(ITT_D[V0.hat]^3)
  var.beta<-(diag(solve(covW)/n)[1:pz])[V0.hat]*temp
  CI.initial<-matrix(NA,nrow=length(V0.hat),ncol=2)
  CI.initial[,1]<-(ITT_Y/ITT_D)[V0.hat]-sqrt(log(n)*var.beta)
  CI.initial[,2]<-(ITT_Y/ITT_D)[V0.hat]+sqrt(log(n)*var.beta)
  uni<- intervals::Intervals(CI.initial)
  CI.initial.union<-base::as.matrix(intervals::interval_union(uni))
  beta.grid.seq<-analysis.CI(CI.initial.union,grid.size=n^{-1})$grid.seq

  ### conduct the initial searching and output a refined range [L,U] ###

  CI.sea<-Searching.CI(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD=SigmaSqD,SigmaSqY=SigmaSqY,
                       SigmaYD=SigmaYD,InitiSet = V0.hat,WUMat = WUMat,alpha = alpha,
                       beta.grid=beta.grid.seq,bootstrap=FALSE)
  CI.temp<-CI.sea$CI.search
  beta.grid<-analysis.CI(base::as.matrix(CI.sea$CI.search),beta,n^{-0.6})$grid.seq

  ### conduct the refined searching ###

  CI.sea.refined<-Searching.CI(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD=SigmaSqD,SigmaSqY=SigmaSqY,
                               SigmaYD=SigmaYD,InitiSet = V0.hat,WUMat = WUMat,alpha = alpha,
                               beta.grid,bootstrap=boot.value)
  CI.search<-CI.sea.refined$CI.search
  CI.temp <- t(base::as.matrix(c(min(CI.search),max(CI.search))))
  ### conduct the refined sampling ###
  if (Sampling) {
    CI.sampling<-Searching.CI.Sampling(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD=SigmaSqD,SigmaSqY=SigmaSqY,
                                       SigmaYD=SigmaYD,InitiSet = V0.hat,WUMat = WUMat,alpha = alpha,
                                       beta.grid,M=M,bootstrap=boot.value,alpha0=alpha0)
    CI.temp<-t(base::as.matrix(c(min(CI.sampling$CI.union[,1]),max(CI.sampling$CI.union[,2]))))
    if (max_clique) {
      return(list(CI.union = CI.temp, rule = CI.sampling$rule,VHat = V0.hat,
                  CI.clique = CI.clique, rule.clique = rule.clique, max.cliques = max.clique.mat))
    } else{

      return(list(CI.union = CI.temp, rule = CI.sampling$rule,VHat = V0.hat))
    }

  } else {
    if (max_clique) {
      return(list(CI.union = CI.temp,CI.matrix = CI.search,
                  rule = CI.sea.refined$rule,VHat = V0.hat,
                  CI.clique = CI.clique, rule.clique = rule.clique, max.cliques = max.clique.mat))
    } else{
      return(list(CI.union = CI.temp,CI.matrix = CI.search,
                  rule = CI.sea.refined$rule,VHat = V0.hat))

    }

  }

}

###### bootstrap threshold
norm.diff<-function(vec,beta.grid,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,SE.norm,pz){
  gamma.s<-vec[1:pz]
  Gamma.s<-vec[-(1:pz)]
  norm.max<-0
  n.beta<-length(beta.grid)
  for(j in 1:n.beta){
    b<-beta.grid[j]
    se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
    temp.diff<-abs(Gamma.s[InitiSet]-b*gamma.s[InitiSet])/(se.b*SE.norm[InitiSet])
    norm.max<-max(norm.max,max(temp.diff))
  }
  return(norm.max)
}
cut.off<-function(SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,pz,alpha = 0.05,
                  beta.grid,N=1000){
  #N=1000
  n=nrow(WUMat)
  unit.matrix<-t(WUMat)%*%WUMat/n^2
  #(solve(covW)/n)[1:pz,1:pz]
  Cov1<-cbind(SigmaSqD*unit.matrix,SigmaYD*unit.matrix)
  Cov2<-cbind(SigmaYD*unit.matrix,SigmaSqY*unit.matrix)
  Cov.total<-rbind(Cov1,Cov2)
  Gen.mat<-MASS::mvrnorm(N, rep(0,2*pz), Cov.total)
  #se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
  SE.norm<-diag(unit.matrix)^{1/2}
  sample.sim<-rep(0,N)
  for(j in 1:N){
    sample.sim[j]<-norm.diff(Gen.mat[j,],beta.grid,SigmaSqD,SigmaSqY,SigmaYD,
                             InitiSet,SE.norm,pz)
  }
  critical.val<-quantile(sample.sim,probs=1-alpha)
  return(critical.val)
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
  critical.val<-quantile(max.vec,probs=cut.prob)
  return(critical.val)
}



##### Functions: handling the (possible) union of CIs
##### CI.matrix is a matrix, where each row represents a CI
analysis.CI<-function(CI.matrix,true.val=1,grid.size){
  #### number of rows in the CI.matrix
  d<-dim(CI.matrix)[1]
  CI.coverage<-0
  CI.len<-0
  grid.seq<-NULL
  for (l in 1: d){
    CI.len<-CI.len+CI.matrix[l,2]-CI.matrix[l,1]
    if((CI.matrix[l,2]>true.val)*(CI.matrix[l,1]<true.val)==1){
      CI.coverage<-1
    }
    grid.seq<-c(grid.seq,seq(CI.matrix[l,1],CI.matrix[l,2],by=grid.size))
  }
  return(list(CI.coverage=CI.coverage,CI.len=CI.len,grid.seq=grid.seq))
}

##### Functions: CI construction by searching method
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        InitiSet: a set of pre-selected IVs (for majority rule: it is the set of relevant IVs;
###        for plurality rule: it is a set of weakly violated IVs)
#' @title Searching method
#' @description CI construction by Searching method
#'
#' @param ITT_Y a numeric, non-missing vector estimating each instrument's effect on the outcome
#' @param ITT_D a numeric, non-missing vector estimating each instrument's effect on the treatment
#' @param SigmaSqD a numeric, non-missing, positive scalar value estimating the error variance of Y in the reduced-form model of Y
#' @param SigmaSqY  a numeric, non-missing, positive scalar value estimating the error variance of D in the reduced-form model of D
#' @param SigmaYD a numeric, non-missing, scalar value estimating the covariance of the error terms in the reduced-form models of Y and D
#' @param InitiSet a set of pre-selected IVs (for majority rule: it is the set of relevant IVs; for plurality rule: it is a set of weakly violated IVs)
#' @param WUMat a numeric matrix denoting WU where U is the precision matrix of W and W is the instrument-covariate matrix (Z, X).
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param beta.grid a numeric vector denoting candidates for betahat.
#' @param bootstrap a boolean scalar indicating to implement bootstrap, with default TRUE.
#'
#' @return
#' \item{\code{CI.search}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints constructed by Searching method.}
#' \item{\code{rule}}{a boolean scalar denoting whether the identification condition is satisfied or not}
#' \item{\code{valid.grid}}{a numeric vector denoting the candidates satisfying certain thresholding condition}
#' @export
#'
Searching.CI<-function(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,
                       alpha = 0.05,beta.grid=NULL,bootstrap=TRUE){
  pz<-ncol(WUMat)
  n = nrow(WUMat)
  if(length(InitiSet)==0){
    warning("No valid IV! OLS is used.")
    CI.search=t(base::as.matrix(stats::confint(lm(Y~D+W))[2,]))
    rule<-FALSE
    return(list(CI.search=CI.search,rule=rule))
  }else{
    threshold.size<-(length(InitiSet)/2)
    if(is.null(beta.grid)){
      beta.grid<-seq(-5,5,by=max(n,500)^{-1})
    }
    n.beta<-length(beta.grid)
    valid.grid<-rep(NA,n.beta)
    SE.norm<-sqrt(SigmaSqD * colSums(WUMat^2) /n^2)
    #(diag(solve(covW)/n)^{1/2})[1:pz]
    if(bootstrap==TRUE){
      Tn<-cut.off(SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,pz,alpha,beta.grid,
                  N=1000)
    }else{
      Tn<-sqrt(2.005*log(n.beta))
    }
    for(j in 1:n.beta){
      b<-beta.grid[j]
      se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
      valid.grid[j]<-sum(abs(ITT_Y[InitiSet]-b*ITT_D[InitiSet])<
                           Tn*se.b*SE.norm[InitiSet])
    }
    #CI.search<-c(min(beta.grid[which(valid.grid>threshold.size)]),max(beta.grid[which(valid.grid>threshold.size)]))
    if(length(beta.grid[which(valid.grid>threshold.size)])==0){
      ####### rule=FALSE indicates the corresponding rule fails
      rule=FALSE
      warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
      sel.index<-which(valid.grid==max(valid.grid))
    }else{
      rule=TRUE
      sel.index<-which(valid.grid>threshold.size)
    }
    CI<-matrix(NA,nrow=length(sel.index),ncol=2)
    CI[,1]<-beta.grid[sel.index]
    upper.index<-sel.index+1
    upper.index[length(upper.index)]<-min(upper.index[length(upper.index)],n.beta)
    CI[,2]<-beta.grid[upper.index]
    if(dim(base::as.matrix(CI))[1]==1)
    {
      CI.search<-base::as.matrix(CI)
    }else{
      uni<- intervals::Intervals(CI)
      ###### construct the confidence interval by taking a union
      CI.search<-base::as.matrix(intervals::interval_union(uni))
      # CI.search <- t(base::as.matrix(c(min(CI.search),max(CI.search)))) #added
    }
    return(list(CI.search=CI.search,rule=rule,valid.grid=valid.grid))
  }
}

###### Sampling and Searching
##### Functions: CI construction by sampling and searching method
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        InitiSet: a set of pre-selected IVs (for majority rule: it is the set of relevant IVs;
###        for plurality rule: it is a set of weakly violated IVs)
#' @title Searching-Sampling method
#' @description CI construction by Searching-Sampling method
#'
#' @param ITT_Y a numeric, non-missing vector estimating each instrument's effect on the outcome
#' @param ITT_D a numeric, non-missing vector estimating each instrument's effect on the treatment
#' @param SigmaSqD a numeric, non-missing, positive scalar value estimating the error variance of Y in the reduced-form model of Y
#' @param SigmaSqY a numeric, non-missing, positive scalar value estimating the error variance of D in the reduced-form model of D
#' @param SigmaYD a numeric, non-missing, scalar value estimating the covariance of the error terms in the reduced-form models of Y and D
#' @param InitiSet a set of pre-selected IVs (for majority rule: it is the set of relevant IVs; for plurality rule: it is a set of weakly violated IVs)
#' @param WUMat a numeric matrix denoting WU where U is the precision matrix of W and W is the instrument-covariate matrix (Z, X).
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param alpha0 a numeric scalar value between 0 and 1 indicating the threshold level for the samples for Sampling method, with default 0.01.
#' @param beta.grid a numeric vector denoting candidates for betahat.
#' @param rho a numeric scalar denoting thresholding level used in the sampling property.
#' @param M a positive integer indicating the number of bootstrap resampling for computing the confidence interval, with default 1000.
#' @param bootstrap a boolean scalar indicating to implement bootstrap, with default TRUE.
#'
#' @return
#' \item{\code{CI.union}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints constructed by Searching-Sampling method.}
#' \item{\code{rule}}{a boolean scalar denoting whether the identification condition is satisfied or not.}
#' \item{\code{CI}}{a numeric matrix denoting the confidence intervals for betaHat constructed by valid candidates for betaHat.}
#' @export
#'
Searching.CI.Sampling<-function(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,InitiSet,
                                WUMat,alpha = 0.05,alpha0 = 0.01,beta.grid=NULL,rho=NULL,
                                M=1000,bootstrap=TRUE){
  pz<-ncol(WUMat)
  n = nrow(WUMat)
  if(length(InitiSet)==0){
    warning("No valid IV! OLS is used.")
    CI.union=t(base::as.matrix(stats::confint(lm(Y~D+W))[2,]))
    rule<-FALSE
    CI=CI.union
    return(list(CI.union=CI.union,rule=rule,CI=CI))
  }else{
    threshold.size<-(length(InitiSet)/2)
    if(is.null(rho)){
      rho<-(log(n)/M)^{1/(2*length(InitiSet))}/6
    }
    if(is.null(beta.grid)){
      beta.grid<-seq(-5,5,by=n^{-0.8})
    }
    n.beta<-length(beta.grid)
    valid.grid.sample<-matrix(NA,nrow=M,ncol=n.beta)
    SE.norm<-sqrt(SigmaSqD * colSums(WUMat^2) /n^2)
    #(diag(solve(covW)/n)^{1/2})[1:pz]
    #Cov.Y<-SigmaSqD*solve(covW)[1:pz,1:pz]/n
    #Cov.D<-SigmaSqY*solve(covW)[1:pz,1:pz]/n
    if(bootstrap==TRUE){
      Tn<-cut.off(SigmaSqD,SigmaSqY,SigmaYD,InitiSet,WUMat,pz,alpha,beta.grid,N=1000)
    }else{
      Tn<-sqrt(2.005*log(n.beta))
    }
    unit.matrix<-t(WUMat)%*%WUMat/n^2
    #(solve(covW)/n)[1:pz,1:pz]
    Cov1<-cbind(SigmaSqD*unit.matrix,SigmaYD*unit.matrix)
    Cov2<-cbind(SigmaYD*unit.matrix,SigmaSqY*unit.matrix)
    Cov.total<-rbind(Cov1,Cov2)

    for(m in 1:M){
      Gen.mat<-MASS::mvrnorm(1, rep(0,2*pz), Cov.total)
      ITT_Y.sample<-ITT_Y-Gen.mat[(pz+1):(2*pz)]
      ITT_D.sample<-ITT_D-Gen.mat[1:pz]
      ###### generating indepedent copies
      #ITT_Y.sample<-ITT_Y-mvrnorm(1,rep(0,pz),Cov.Y)
      #ITT_D.sample<-ITT_D-mvrnorm(1,rep(0,pz),Cov.D)

      for(j in 1:n.beta){
        b<-beta.grid[j]
        se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
        valid.grid.sample[m,j]<-sum(abs(ITT_Y.sample[InitiSet]-b*ITT_D.sample[InitiSet])<rho*Tn*se.b*SE.norm[InitiSet])
      }

    }
    CI<-matrix(NA,nrow=M,ncol=2)
    for(m in 1:M){
      if(length(which(valid.grid.sample[m,]>threshold.size))>0){
        CI[m,1]<-min(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
        CI[m,2]<-max(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
      }
    }
    CI<-stats::na.omit(CI)

    delta <- 0.05
    while((dim(base::as.matrix(CI))[1]<min(0.05*M,50)) && (rho<0.5)){
      M0 <- NULL
      #print(base::as.matrix(CI)[1])
      #print(rho)
      rho<-1.25*rho
      for(m in 1:M){
        Gen.mat<-MASS::mvrnorm(1, rep(0,2*pz), Cov.total)
        ITT_Y.sample<-ITT_Y-Gen.mat[(pz+1):(2*pz)]
        ITT_D.sample<-ITT_D-Gen.mat[1:pz]

        #ITT_Y.sample<-ITT_Y-mvrnorm(1,rep(0,pz),Cov.Y)
        #ITT_D.sample<-ITT_D-mvrnorm(1,rep(0,pz),Cov.D)
        for(j in 1:n.beta){
          b<-beta.grid[j]
          se.b<-sqrt(SigmaSqY+b^2*SigmaSqD-2*b*SigmaYD)
          valid.grid.sample[m,j]<-sum(abs(ITT_Y.sample[InitiSet]-b*ITT_D.sample[InitiSet])<rho*Tn*se.b*SE.norm[InitiSet])
        }

        if (t(Gen.mat)%*%solve(Cov.total,Gen.mat)<=(1+delta)*stats::qchisq(p = 1-alpha0,df = 2*pz)) {
          M0 <- c(M0,m)
        }

      }

      CI<-matrix(NA,nrow=M,ncol=2)
      for(m in M0){
        if(length(which(valid.grid.sample[m,]>threshold.size))>0){

          CI[m,1]<-min(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
          CI[m,2]<-max(beta.grid[which(valid.grid.sample[m,]>threshold.size)])
        }
      }
      CI<-stats::na.omit(CI)
    }
    rule<-TRUE
    if(dim(base::as.matrix(CI))[1]==0){
      rule<-FALSE
      CI.union<-t(base::as.matrix(c(min(beta.grid),max(beta.grid))))
    }else if(dim(base::as.matrix(CI))[1]==1)
    {
      CI.union<-base::as.matrix(CI)
    }else{
      uni<- intervals::Intervals(CI)
      ###### construct the confidence interval by taking a union
      CI.union<-base::as.matrix(intervals::interval_union(uni))
      # CI.union <- t(base::as.matrix(c(min(CI.union),max(CI.union)))) # added

    }
    #CI.upper<-quantile(CI[,2],0.975)
    #CI.lower<-quantile(CI[,1],0.025)
    #CI.quantile<-c(CI.lower,CI.upper)
    return(list(CI.union=CI.union,rule=rule,CI=CI))
  }
}

### TSHT.Initial (This modifies the original TSHT function and produces an initial esitmator for sampling and searching)
### FUNCTION: Estimates the set of valid and relevant IVs
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        WUMat: an n by pz numeric, non-missing matrix that computes the precision of W
###               (ex1 t(WUMat) %*% WUMat / n = (W^T W/n)^(-1))
###               (ex2 U = (W^T W/n)^(-1))
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        tuning, a numeric scalar value tuning parameter for TSHT greater
###                than 2 (default = 2.01)
### OUTPUT: a list (a) VHat (numeric vector denoting the set of valid and relevant IVs)
###                (b) SHat (numeric vector denoting the set of relevant IVs)
#' @title Two Stage Hard Thresholding for initial estimator.
#' @description This modifies the original TSHT function and produces an initial esitmator for sampling and searching
#'
#' @param ITT_Y a numeric, non-missing vector estimating each instrument's effect on the outcome
#' @param ITT_D a numeric, non-missing vector estimating each instrument's effect on the treatment
#' @param WUMat an n by pz numeric, non-missing matrix that computes the precision of W
#' @param SigmaSqY a numeric, non-missing, positive scalar value estimating the error variance of Y in the reduced-form model of Y
#' @param SigmaSqD a numeric, non-missing, positive scalar value estimating the error variance of D in the reduced-form model of D
#' @param SigmaYD a numeric, non-missing, scalar value estimating the covariance of the error terms in the reduced-form models of Y and D
#' @param covW a numeric, non-missing matrix that computes the sample covariance of W
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05.
#' @param bootstrap a boolean scalar indicating to implement bootstrap, with default TRUE.
#' @param tuning a numeric scalar value tuning parameter for TSHT greater than 2 (default = 2.01)
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix, with default FALSE.
#'
#' @return
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs.}
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs.}
#'     \item{\code{max.clique}}{a numeric list denoting the maximum cliques of valid and relevant IVs. Only return when \code{max_clique} is \code{TRUE}.}
#' @export
#'
#'
TSHT.Initial <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,covW, alpha = 0.05,
                         bootstrap=FALSE,tuning = 2.01,max_clique=FALSE) {
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
    Tn<-min(cut.off.IVStr(SigmaSqD,WUMat,pz),sqrt(log(n))) ### this can be modified by the user
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

  # maximal clique
  if (max_clique) {
    voting.graph <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(VHats.boot.sym))
    max.clique <- igraph::largest.cliques(voting.graph)
    VHat <- unique(igraph::as_ids(Reduce(c,max.clique))) # take the union if multiple max cliques exist
    VHat <- sort(as.numeric(VHat))
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
    VHat=sort(as.numeric(VHat))
  }

  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }
  if (max_clique) {
    returnList <- list(SHat=SHat,VHat=VHat,max.clique=max.clique)
  } else {
    returnList <- list(SHat=SHat,VHat=VHat)
  }
  return(returnList)
}
