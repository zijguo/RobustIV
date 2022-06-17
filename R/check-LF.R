check.args.LF <- function(X=NULL, y=NULL, loading.mat=NULL, model=NULL, intercept=NULL, intercept.loading=NULL,
                          lambda=NULL, mu=NULL, init.step=NULL, resol=NULL, maxiter=NULL, alpha=NULL,
                          verbose=NULL){
  if(is.null(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if(is.null(y) || !is.numeric(y)) stop("y must be a numeric vector")
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) must match")
  if(is.null(loading.mat) || !is.numeric(loading.mat)) stop("loading must be a numeric matrix")
  if(ncol(X) != nrow(loading.mat)) stop("ncol(X) and nrow(loading) must match")
  if(is.null(model) || !model%in%c("linear","logistic","logistic_alternative","probit")){
    stop("specified model has not been developed yet")
  }
  if(is.null(intercept) || length(intercept)!=1 || !is.logical(intercept)){
    stop("intercept must be a Boolean")
  }
  if(is.null(intercept.loading) || length(intercept.loading)!=1 || !is.logical(intercept)){
    stop("intercept must be a Boolean")
  }
  if(!is.null(lambda)){
    if(!(lambda %in% c("CV.min","CV")) || length(lambda)!=1 || !is.numeric(lambda)) stop("lambda must be a number, except from 'CV.min','CV'.")
  }
  if(!is.null(mu)) if(!is.numeric(mu) || length(mu)!=1) stop("mu must be a number")
  if(!is.null(init.step) && (!is.integer(init.step) || length(init.step)!=1 || init.step <= 0)){
    stop("If specified, init.step must be a positive integer")
  }
  if(is.null(resol) || !is.numeric(resol) || length(resol)!=1 || resol < 0){
    stop("resol must be a nonnegative number")
  }
  if(is.null(maxiter) || !is.numeric(maxiter) || length(maxiter)!=1 || maxiter <= 0){
    stop("maxiter must be a positive integer")
  }
  if(is.null(alpha) || length(alpha)!= 1 || !is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("alpha must be a number between 0 and 1")
  }
  if(is.null(verbose) || !is.logical(verbose) || length(verbose)!=1 ) stop("verbose must be a Boolean")
}

check.args.QF <- function(X=NULL, y=NULL, G=NULL, A=NULL, model=NULL, 
                          intercept=NULL, tau.vec=NULL, lambda=NULL, mu=NULL, init.step=NULL, resol=NULL, 
                          maxiter=NULL, alpha=NULL, verbose=NULL){
  if(is.null(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if(is.null(y) || !is.numeric(y)) stop("y must be a numeric vector")
  if(nrow(X) != length(y)) stop("nrow(X) and length(y) must match")
  if(is.null(G) || !is.numeric(G)) stop("G must be a numeric vector")
  if(any(G<1) || any(G>ncol(X)) || length(G)!=length(unique(G))) stop("G must be distinct indexes between 1 and ncol(X)")
  if(!is.null(A)) if(!is.numeric(A) || nrow(A)!=ncol(A) || nrow(A)!=length(G)) stop("A must be a numeric square matrix")
  if(is.null(model) || !model%in%c("linear","logistic","logistic_alternative","probit")){
    stop("specified model has not been developed yet")
  }
  if(is.null(intercept) || length(intercept)!=1 || !is.logical(intercept)){
    stop("intercept must be a Boolean")
  }
  if(is.null(tau.vec) || !is.numeric(tau.vec) || any(tau.vec < 0)) stop("tau.vec must be a numeric vector")
  if(!is.null(lambda)){
    if(!(lambda %in% c("CV.min","CV")) || length(lambda)!=1 || !is.numeric(lambda)) stop("lambda must be a number, except from 'CV.min','CV'.")
  }
  if(!is.null(mu)) if(!is.numeric(mu) || length(mu)!=1) stop("mu must be a number")
  if(!is.null(init.step) && (!is.integer(init.step) || length(init.step)!=1 || init.step <= 0)){
    stop("If specified, init.step must be a positive integer")
  }
  if(is.null(resol) || !is.numeric(resol) || length(resol)!=1 || resol < 0){
    stop("resol must be a nonnegative number")
  }
  if(is.null(maxiter) || !is.numeric(maxiter) || length(maxiter)!=1 || maxiter <= 0){
    stop("maxiter must be a positive integer")
  }
  if(is.null(alpha) || length(alpha)!= 1 || !is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("alpha must be a number between 0 and 1")
  }
  if(is.null(verbose) || !is.logical(verbose) || length(verbose)!=1 ) stop("verbose must be a Boolean")
}
