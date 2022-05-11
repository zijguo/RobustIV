
TSHT_hetero <- function(Y, D, Z, X, intercept=TRUE, alpha = 0.05, tuning=2.01,
                        method="OLS", voting = 'MaxClique'){
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
  stopifnot(method=='OLS')
  stopifnot(voting=='MP' | voting=='MaxClique' | voting == 'Conservative')


  # Derive Inputs for TSHT, consider only OLS setting now
  n = nrow(Z); pz=ncol(Z)
  inputs = TSHT.OLS_hetero(Y, D, W, pz, intercept)

  # Estimate Valid IVs
  SetHats = TSHT.VHat_hetero(n, inputs$ITT_Y, inputs$ITT_D,
                             inputs$V.Gamma, inputs$V.gamma, inputs$C, tuning, voting)
  VHat = SetHats$VHat; SHat = SetHats$SHat
  check = T
  if(length(VHat)< length(SHat)/2){
    cat('Majority rule fails.','\n')
    check=F
  }
  # Obtain point est, se and ci
  ITT_Y = inputs$ITT_Y;
  ITT_D = inputs$ITT_D;
  WUMat = inputs$WUMat;
  V.Gamma = inputs$V.Gamma
  V.gamma = inputs$V.gamma
  C = inputs$C
  A = t(WUMat) %*% WUMat / n

  if (voting == 'MaxClique') {
    max.clique <- SetHats$max.clique
    max.clique.mat <- matrix(0,nrow = length(max.clique),ncol = length(max.clique[[1]]))
    CI.temp <- matrix(0,nrow = length(max.clique), ncol = 2)
    beta.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    betavar.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    for (i in 1:length(max.clique)) {
      temp <- SHat[sort(as.numeric(max.clique[[i]]))]
      max.clique.mat[i,] <- temp
      AVHat = solve(A[temp,temp])
      betaHat = as.numeric((t(ITT_Y[temp]) %*% AVHat %*% ITT_D[temp]) / (t(ITT_D[temp])  %*% AVHat %*% ITT_D[temp]))
      temp2 <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[temp,temp]
      AVHat = solve(temp2)
      betaHat = as.numeric((t(ITT_Y[temp]) %*% AVHat %*% ITT_D[temp]) / (t(ITT_D[temp])  %*% AVHat %*% ITT_D[temp]))
      temp2 <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[temp,temp]
      betaVar <- t(ITT_D[temp])%*% AVHat %*%temp2%*% AVHat %*%ITT_D[temp] / (n*(t(ITT_D[temp])%*% AVHat %*%ITT_D[temp])^2)
      ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVar),betaHat + qnorm(1-alpha/2) * sqrt(betaVar))
      CI.temp[i,] <- ci
      beta.temp[i,] <- betaHat
      betavar.temp[i,] <- betaVar
    }
    uni<- intervals::Intervals(CI.temp)
    ###### construct the confidence interval by taking a union
    CI.union<-as.matrix(intervals::interval_union(uni))
    # CI.union <- t(as.matrix(c(min(CI.union),max(CI.union)))) # added
    # ci <- as.vector(CI.union)
  }


  AVHat = solve(A[VHat,VHat])
  betaHat = as.numeric((t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat])  %*% AVHat %*% ITT_D[VHat]))
  temp <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[VHat,VHat]
  AVHat = solve(temp)
  betaHat = as.numeric((t(ITT_Y[VHat]) %*% AVHat %*% ITT_D[VHat]) / (t(ITT_D[VHat])  %*% AVHat %*% ITT_D[VHat]))
  temp <- (V.Gamma -2*betaHat*C+betaHat^2*V.gamma)[VHat,VHat]
  betaVar <- t(ITT_D[VHat])%*% AVHat %*%temp%*% AVHat %*%ITT_D[VHat] / (n*(t(ITT_D[VHat])%*% AVHat %*%ITT_D[VHat])^2)
  # U-2*betaHat*UV+betaHat^2*V
  ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVar),betaHat + qnorm(1-alpha/2) * sqrt(betaVar))
  if (voting != 'MaxClique') {
    TSHT.model <- list( betaHat=betaHat,beta.sdHat = sqrt(betaVar),ci=ci,SHat=SHat,VHat = VHat,voting.mat=SetHats$voting.mat,check = check)
  } else {
    TSHT.model <- list( betaHat=betaHat,beta.sdHat = sqrt(betaVar),ci=CI.union,SHat=SHat,VHat = VHat,voting.mat=SetHats$voting.mat, check = check,
                        beta.clique = beta.temp,beta.sd.clique = sqrt(betavar.temp), CI.clique = CI.temp,
                        max.clique = max.clique.mat)
  }
  return(TSHT.model)
}

TSHT.OLS_hetero <- function(Y, D, W, pz, intercept=TRUE){
  n = nrow(W)
  if(intercept) W = cbind(W, 1)
  p = ncol(W)
  covW = t(W)%*%W/n
  U = solve(covW) # precision matrix
  WUMat = (W%*%U)[,1:pz]
  ## OLS estimators
  qrW = qr(W)
  ITT_Y = qr.coef(qrW, Y)[1:pz]
  ITT_D = qr.coef(qrW, D)[1:pz]
  resid_Y = as.vector(qr.resid(qrW, Y))
  resid_D = as.vector(qr.resid(qrW, D))
  V.Gamma = (t(WUMat)%*%diag(resid_Y^2)%*%WUMat)/n
  V.gamma = (t(WUMat)%*%diag(resid_D^2)%*%WUMat)/n
  C = (t(WUMat)%*%diag(resid_Y * resid_D)%*%WUMat)/n

  out <- list(ITT_Y = ITT_Y,
              ITT_D = ITT_D,
              V.Gamma = V.Gamma,
              V.gamma = V.gamma,
              C = C,
              WUMat = WUMat)
  return(out)
}

TSHT.VHat_hetero <- function(n, ITT_Y, ITT_D, V.Gamma, V.gamma, C, tuning = 2.01,voting = 'MaxClique'){
  pz = nrow(V.Gamma)
  ## First Stage
  Tn = max(sqrt(tuning*log(pz)), sqrt(log(n)/2))
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

    PHat.bool.j = abs(pi.j) <= sqrt(SE.j)*sqrt(tuning^2*log(pz))
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



