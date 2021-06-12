#' @title Causal inference in nonlinear outcome models with valid control functions
#' @description Causal inference in nonlinear outcome model with IVs, which assumes all the IVs are valid.
#'
#' @param Y binary and non-missing, n by 1 numeric outcome vector.
#' @param D continuous and non-missing, n by 1 numeric treatment vector.
#' @param Z continuous or discrete, non-missing, n by pz numeric instrument matrix, containing p_z instruments.
#' @param intercept a boolean scalar indicating to include the intercept or not, with default TRUE.
#' @param d1 a scalar for computing CATE(d1,d2|z0).
#' @param d2 a scalar for computing CATE(d1,d2|z0).
#' @param z0  a pz by 1 vector for computing CATE(d1,d2|z0).
#' @param bs.Niter a positive integer indicating the number of bootstrap resampling for computing the confidence interval.
#' @param bw  a 2 by 1 vector bandwidth specification; default is NULL and the bandwidth is chosen by rule of thumb.
#'
#' @return
#'     \item{\code{cateHat}}{a numeric scalar denoting the estimate of CATE(d1,d2|z0).}
#'     \item{\code{cate.sdHat}}{a numeric scalar denoting the estimated standard deviation of cateHat.}
#' @export
#'
#'
#'
#' @examples
#' \dontrun{
#' ### Working Low Dimensional Example ###
#' library(mvtnorm)
#' library(MASS)
#' n = 1000; J = 7; s = 5; d1=-1; d2=1; z0=c(rep(0, J-1),0.1)
#' Z <- matrix(rnorm(n * J, 0, 1) , ncol = J, nrow = n)
#' gam <- c(rep(0.8, floor(J / 2)), rep(-0.8, J - floor(J / 2)))
#' cov.noise<-matrix(c(1,0.25, 0.25, 1),ncol=2)
#' noise.vec<-rmvnorm(n, rep(0,2), cov.noise)
#' v.vec<-noise.vec[,1]
#' D = 0.5+Z %*% gam  + v.vec
#'
#' pi0 <- c(rep(0, s), 0.4, 0.2)
#' beta0 <- 0.25
#' u.vec<- noise.vec[,2]
#' Y = (-0.5 + Z %*% pi0 + D * beta0 + u.vec>=0)
#' u1.r<-rnorm(2000,0,sd=1)
#' cace0 <- mean((as.numeric(-0.5+d1 * beta0 + z0 %*% pi0)+ u1.r )>=0) - mean((as.numeric(-0.5+d2 * beta0 + z0 %*% pi0) + u1.r)>=0)
#' SemiControl.valid(Y=Y, D=D, Z=Z, bs.Niter = 40, d1 = d1, d2 = d2, z0 = z0)
#'}
#'
#'



SemiControl.valid<- function(Y, D, Z, bs.Niter=40, intercept=T, d1, d2 , z0, bw=NULL){
  pz<- ncol(Z)
  stopifnot(length(z0)==ncol(Z))
  if(intercept){ Z<-cbind(Z,1) }
  n <- length(Y)
  gam.hat<- lm(D ~ Z-1)$coef
  v.hat<-D-Z%*%gam.hat
  bw.z<- apply(cbind(D,v.hat), 2, function(x) 0.9*n^(-1/6)*min(sd(x), (quantile(x,0.75)-quantile(x,0.25))/1.34))
  asf.dw<- VC.ASF.est(d1=d1,d2=d2, z0=z0, Y, D, v.hat=v.hat,bw.z=bw.z)
  cace.hat = asf.dw$cace.hat #cate
  bw.z=asf.dw$bw.z
  cat(bw.z,'\n')
  ####bootstrap####
  boot_b<-list()
  for(i in 1: bs.Niter){
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    boot_b[[i]]<-VC.boot.fun(data=bootstrap_data, d1=d1,d2=d2, z0=z0, bw.z=bw.z)
  }
  cace.sd<-sqrt(mean((unlist(lapply(boot_b, function(x) x[1]))-cace.hat)^2))
  return(list(cace.hat=cace.hat, cace.sd= cace.sd))
}


VC.boot.fun<-function(data, d1, d2,z0, bw.z=NULL){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  pz<-ncol(Z)
  v.bs <- lm(D~Z-1)$res
  asf.dw<-VC.ASF.est(d1=d1,d2=d2, z0=z0, Y=Y, D=D, v.hat=v.bs, bw.z=bw.z)
  asf.dw$cace.hat

}

box.ker<-function(xx, X, h){
  apply(X, 1, function(xi) max(abs(xi-xx)/h)<1)/prod(h)
}

VC.ASF.est <- function(d1, d2, z0, Y, D, v.hat, bw.z=NULL){
  n=length(D)
  v.hat <- as.matrix(v.hat,ncol=1)
  D<-as.matrix(D,ncol=1)
  index <- cbind(D,v.hat)
  #cace
  index1.z <- cbind(rep(d1,n),v.hat)
  index2.z <- cbind(rep(d2,n),v.hat)
  asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] ])), na.rm=T)
  asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2]])), na.rm=T)

  cace.hat <- asf.dw1-asf.dw2
  list(cace.hat = cace.hat, bw.z=bw.z)
}


