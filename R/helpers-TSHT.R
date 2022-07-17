TSHT.OLS <- function(Y,D,W,pz,intercept=TRUE) {

  n = nrow(W)
  if(intercept) W = cbind(W, 1)
  p = ncol(W)
  covW = t(W)%*%W/n
  U = solve(covW) # precision matrix
  WUMat = (W%*%U)[,1:pz]
  ## OLS estimators
  qrW = qr(W)
  ITT_Y = Matrix::qr.coef(qrW, Y)[1:pz]
  ITT_D = Matrix::qr.coef(qrW, D)[1:pz]
  resid_Y = as.vector(Matrix::qr.resid(qrW, Y))
  resid_D = as.vector(Matrix::qr.resid(qrW, D))
  SigmaSqY = sum(Matrix::qr.resid(qrW,Y)^2)/(n -p)
  SigmaSqD = sum(Matrix::qr.resid(qrW,D)^2)/(n -p)
  SigmaYD = sum(Matrix::qr.resid(qrW,Y) * Matrix::qr.resid(qrW,D)) / (n - p)
  ## V and C below are results for robust=TRUE
  # V.Gamma = (t(WUMat)%*%diag(resid_Y^2)%*%WUMat)/n
  # V.gamma = (t(WUMat)%*%diag(resid_D^2)%*%WUMat)/n
  # C = (t(WUMat)%*%diag(resid_Y * resid_D)%*%WUMat)/n
  V.Gamma = crossprod(resid_Y*WUMat)
  V.gamma = crossprod(resid_D*WUMat)
  C = crossprod(resid_Y*WUMat, resid_D*WUMat)

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,V.gamma = V.gamma, V.Gamma = V.Gamma, C = C, SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
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

TSHT.SIHR <- function(Y, D, W, pz, intercept=TRUE){
  n = nrow(W)
  covW = t(W)%*%W / n
  init_Y = Lasso(W, Y, lambda="CV.min", intercept=intercept)
  init_D = Lasso(W, D, lambda="CV.min", intercept=intercept)

  if(intercept) W_int = cbind(1, W) else W_int = W
  resid_Y = as.vector(Y - W_int%*%init_Y)
  resid_D = as.vector(D - W_int%*%init_D)

  loading.mat = matrix(0, nrow=ncol(W), ncol=pz)
  for(i in 1:pz) loading.mat[i, i] = 1
  out1 = LF(W, Y, loading.mat, model="linear", intercept=intercept, intercept.loading=FALSE, verbose=TRUE)
  out2 = LF(W, D, loading.mat, model="linear", intercept=intercept, intercept.loading=FALSE, verbose=TRUE)
  ITT_Y = out1$est.debias.vec
  ITT_D = out2$est.debias.vec
  U = out2$proj.mat
  WUMat = W_int%*%U

  SigmaSqY = sum(resid_Y^2)/n
  SigmaSqD = sum(resid_D^2)/n
  SigmaYD = sum(resid_Y * resid_D)/n

  return(list(ITT_Y = ITT_Y,
              ITT_D = ITT_D,
              WUMat = WUMat,
              SigmaSqY = SigmaSqY,
              SigmaSqD = SigmaSqD,
              SigmaYD = SigmaYD))

}


TSHT.VHat <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, voting = 'MaxClique', method='OLS'){
  pz = nrow(V.Gamma)
  if(method=="OLS"){
    Tn = sqrt(log(n))
  }else{
    Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)))
  }
  ## First Stage
  # Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)))
  SHat = (1:pz)[abs(ITT_D) > (Tn * sqrt(diag(V.gamma)/n))]

  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak.
            TSHT with these IVs will give misleading CIs, SEs, and p-values.
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE, pz); SHat.bool[SHat] = TRUE

  ## Second Stage
  nCand = length(SHat)
  VHats.bool = matrix(FALSE, nCand, nCand)
  colnames(VHats.bool) = rownames(VHats.bool) = SHat

  for(j in SHat){
    beta.j = ITT_Y[j]/ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    Temp = V.Gamma + beta.j^2*V.gamma - 2*beta.j*C
    SE.j = rep(NA, pz)
    for(k in 1:pz){
      SE.j[k] = 1/n * (Temp[k,k] + (ITT_D[k]/ITT_D[j])^2*Temp[j,j] -
                         2*(ITT_D[k]/ITT_D[j])*Temp[k,j])
    }

    PHat.bool.j = abs(pi.j) <= sqrt(SE.j)*Tn #sqrt(log(n))
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat), as.character(j)] = VHat.bool.j[SHat]
  }
  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }
  diag(VHats.boot.sym) = 1

  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
  VHat = as.numeric(union(VM.m, VM.p))

  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }
  if (voting == 'MaxClique') {
    voting.graph <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(VHats.boot.sym))
    max.clique <- igraph::largest.cliques(voting.graph)
    # VHat <- unique(igraph::as_ids(Reduce(c,max.clique))) # take the union if multiple max cliques exist
    # VHat <- sort(as.numeric(VHat))
    VHat = lapply(max.clique, FUN=function(x) sort(as.numeric(x)))
    n.VHat = length(VHat[[1]])
    if(length(VHat)==1) VHat = VHat[[1]]
  } else if (voting == 'MP') {
    VHat <- sort(as.numeric(union(VM.m,VM.p))) # Union of majority and plurality winners
    n.VHat = length(VHat)
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
    n.VHat = length(VHat)
  }

  # Error check
  if(n.VHat == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }

  returnList <- list(SHat=SHat,VHat=VHat,voting.mat=VHats.boot.sym)
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
  TSHT1 <- object
  if (typeof(TSHT1$VHat)=="list") {
    result<-matrix(NA, ncol=5, nrow=length(TSHT1$VHat))
    result <- data.frame(result)
    colnames(result)<-c("betaHat","Std.Error",paste("CI(",round(TSHT1$alpha/2*100, digits=2), "%)", sep=""),
                        paste("CI(",round((1-TSHT1$alpha/2)*100, digits=2), "%)", sep=""),"Valid IVs")
    rownames(result)<-paste0("MaxClique",1:length(TSHT1$VHat))
    result[,1] <- unlist(TSHT1$betaHat)
    result[,2] <- unlist(TSHT1$beta.sdHat)
    result[,3:4] <- matrix(unlist(TSHT1$ci),nrow = length(TSHT1$VHat),ncol = 2,byrow = T)
    for (i in 1:length(TSHT1$VHat)) {
      result[i,5] <- paste(TSHT1$VHat[[i]], collapse = " ")
    }
    cat("\nRelevant IVs:", TSHT1$SHat, "\n");
    cat(rep("_", 30), "\n")
    print(result,right=F)
  } else {
    cat("\nRelevant IVs:", TSHT1$SHat, "\n");
    cat(rep("_", 30), "\n")
    result<-matrix(NA, ncol=5, nrow=1)
    result <- data.frame(result)
    colnames(result)<-c("betaHat","Std.Error",paste("CI(",round(TSHT1$alpha/2*100, digits=2), "%)", sep=""),
                        paste("CI(",round((1-TSHT1$alpha/2)*100, digits=2), "%)", sep=""),"Valid IVs")
    rownames(result)<-""
    result[,1] <- TSHT1$betaHat
    result[,2] <- TSHT1$beta.sdHat
    result[,3:4] <- TSHT1$ci
    result[,5] <- paste(TSHT1$VHat,collapse = " ")
    print(result,right=F)
  }
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);

  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }

  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;

  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){

    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }

    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }

  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }

  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

Lasso <- function( X, y, lambda = NULL, intercept = TRUE){
  #
  # Compute the Lasso estimator:
  # - If lambda is given, use glmnet and standard Lasso
  # - If lambda is not given, use square root Lasso
  #
  p <- ncol(X);
  n <- nrow(X);

  if  (is.null(lambda)){
    lambda <- sqrt(qnorm(1-(0.1/p))/n);
    outLas <- flare::slim(X,y,lambda=c(lambda),method="lq",q=2,verbose=FALSE);
    # Objective : sqrt(RSS/n) +lambda *penalty
    if (intercept==TRUE) {
      return (c(as.vector(outLas$intercept),as.vector(outLas$beta)))
    }  else {
      return (as.vector(outLas$beta));
    }
  } else {
    # outLas <- cv.glmnet(X, y, family ="gaussian", alpha =1, intercept = intercept, standardize=T);
    # # Objective :1/2 RSS/n +lambda *penalty
    # htheta = if(lambda=="CV.min"){
    #   as.vector(coef(outLas, s=outLas$lambda.min))
    # }else if(lambda=="CV"){
    #   as.vector(coef(outLas, s=outLas$lambda.1se))
    # }else{
    #   as.vector(coef(outLas, s=lambda))
    # }
    htheta <- if (lambda == "CV.min") {
      outLas <- glmnet::cv.glmnet(X, y, family = "gaussian", alpha = 1,
                                  intercept = intercept, standardize = T)
      as.vector(coef(outLas, s = outLas$lambda.min))
    } else if (lambda == "CV") {
      outLas <- glmnet::cv.glmnet(X, y, family = "gaussian", alpha = 1,
                                  intercept = intercept, standardize = T)
      as.vector(coef(outLas, s = outLas$lambda.1se))
    } else {
      outLas <- glmnet::glmnet(X, y, family = "gaussian", alpha = 1,
                               intercept = intercept, standardize = T)
      as.vector(coef(outLas, s = lambda))
    }
    if (intercept==TRUE){
      return (htheta);
    } else {
      return (htheta[2:(p+1)]);
    }
  }
}

SSLasso <- function (X, y, lambda = NULL, mu = NULL, intercept = TRUE,
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = TRUE) {
  #
  # Compute confidence intervals and p-values.
  #
  # Args:
  #   X     :  design matrix
  #   y     :  response
  #   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
  #   mu    :  Linfty constraint on U (if null, searches)
  #   intercept: Should the intercept term be included?
  #   resol :  step parameter for the function that computes U
  #   maxiter: iteration parameter for computing U
  #   threshold : tolerance criterion for computing U
  #   verbose : verbose?
  #
  # Returns:
  #   coef    : Lasso estimated coefficients
  #   unb.coef: Unbiased coefficient estimates
  #   WUMat: projection of the inverse covariance matrix.
  #   resid.lasso: residual based on Lasso
  #
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);

  # Solve Lasso problem using FLARE package
  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);

  # Format design matrix to include intercept and standardize.
  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  resid.lasso = (y - Xb %*% htheta)
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);

  # Estimation of U (or M in Javanard and Montanari's original paper)
  # Check to see if this is a low dimensional problem
  if ((n>=2*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }

  # If low-dimensional problem, use inverse of covariance as an estiamte of precision matrix
  # Otherwise, solve the convex optimizatio problem for U
  if ((n>=2*p)&&(tmp>=1e-4)){
    U <- solve(sigma.hat)
  }else{
    U <- InverseLinfty(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  }

  # Debias Lasso
  unbiased.Lasso <- as.numeric(htheta + (U%*%t(Xb)%*%(y - Xb %*% htheta))/n);

  # Scale them back to the original scaling.
  htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm;
  WUMat = Xb %*% t(U) %*% diag(col.norm)

  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    WUMat <- WUMat[,2:pp]
  }

  returnList <- list("coef" = htheta,
                     "unb.coef" = unbiased.Lasso,
                     "WUMat" = WUMat,
                     "resid.lasso" = resid.lasso)
  return(returnList)
}
