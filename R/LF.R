LF <- function(X, y, loading.mat, model=c("linear","logistic","logistic_alternative","probit"),
               intercept=TRUE, intercept.loading=TRUE, lambda=NULL,
               mu=NULL, init.step=NULL, resol=1.5, maxiter=6, alpha=0.05,
               verbose=TRUE){
  model = match.arg(model)
  X = as.matrix(X)
  y = as.vector(y)
  loading.mat = as.matrix(loading.mat)
  nullmu = ifelse(is.null(mu), TRUE, FALSE)

  ### Check arguments ###
  if(!is.logical(verbose)) verbose=TRUE
  if(intercept==FALSE && intercept.loading==TRUE){
    intercept.loading = FALSE
    cat("Argument 'intercept.loading' is set to FALSE, because 'intercept' is FALSE")
  }
  check.args.LF(X=X, y=y, loading.mat=loading.mat, model=model, intercept=intercept,
             intercept.loading=intercept.loading, lambda=lambda, mu=NULL, init.step=NULL,
             resol=resol, maxiter=maxiter, alpha=alpha, verbose=verbose)
  null_initstep = ifelse(is.null(init.step), TRUE, FALSE)

  ### specify relevant functions ###
  funs.all = relevant.funs(intercept=intercept, model=model)
  train.fun = funs.all$train.fun
  pred.fun = funs.all$pred.fun
  deriv.fun = funs.all$deriv.fun
  weight.fun = funs.all$deriv.fun
  cond_var.fun = funs.all$cond_var.fun

  ### centralize X ###
  X_means = colMeans(X)
  X = scale(X, center=TRUE, scale=F)

  ### Initial lasso estimator of beta ###
  beta.init = train.fun(X, y, lambda=lambda)$lasso.est

  ### prepare values ###
  if(intercept) X = cbind(1, X)
  n = nrow(X); p = ncol(X)
  pred = as.vector(pred.fun(X%*%beta.init))
  deriv = as.vector(deriv.fun(X%*%beta.init))
  weight = as.vector(weight.fun(X%*%beta.init))
  cond_var = as.vector(cond_var.fun(pred, y))

  ### storing infos ###
  n.loading = ncol(loading.mat)
  est.plugin.vec = rep(NA, n.loading)
  est.debias.vec = rep(NA, n.loading)
  se.vec = rep(NA, n.loading)
  ci.mat = matrix(NA, nrow = n.loading, ncol = 2)
  colnames(ci.mat) = c("lower","upper"); rownames(ci.mat) = paste("loading",1:n.loading,sep="")
  proj.mat = matrix(NA, nrow=p, ncol=n.loading)

  for(i.loading in 1:n.loading){
    ### adjust loading ###
    loading = as.vector(loading.mat[,i.loading])
    if(intercept){
      if(intercept.loading){
        loading = loading - X_means
        loading = c(1, loading)
      }else{
        loading = c(0, loading)
      }
    }
    loading.norm = sqrt(sum(loading^2))

    ############### Correction Direction #################
    if(verbose) cat(sprintf("Computing LF for loading (%i/%i)... \n", i.loading, n.loading))
    if (n >= 6*p){
      temp = weight*deriv*X
      Sigma.hat = t(temp)%*%temp / n
      direction = solve(Sigma.hat) %*% loading / loading.norm
    } else {
      ### find init.step ###
      if(null_initstep){
        step.vec <- rep(NA,3)
        for(t in 1:3){
          index.sel <- sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
          Direction.Est.temp <-  Direction_searchtuning(X[index.sel,], loading, weight = weight[index.sel], deriv.vec = deriv[index.sel], resol, maxiter)
          step.vec[t] <- Direction.Est.temp$step
        }
        init.step<- getmode(step.vec)
      }
      ### for loop to find direction ###
      for(step in init.step:1){
        if(verbose) cat(sprintf("---> Finding Direction with step: %s \n", step))
        if(nullmu) mu = sqrt(2.01*log(p)/n)*resol^{-(step-1)}
        Direction.Est <-  Direction_fixedtuning(X, loading, mu = mu, weight = weight, deriv.vec = deriv)
        if(is.na(Direction.Est)|| length(Direction.Est$proj)==0){
          step = step - 1
        }else{
          if(verbose) cat(sprintf("---> Direction is identified at step: %s \n", step))
          direction <- Direction.Est$proj
          break
        }
      }
    }

    ############## Bias Correction ###############
    est.plugin = sum(beta.init * loading)
    correction = mean((weight * (y-pred) * X) %*% direction)
    est.debias = est.plugin + correction*loading.norm

    ############## Compute SE and Construct CI ###############
    V = sum(((sqrt(weight^2 * cond_var) * X) %*% direction)^2)/n * loading.norm^2
    se = sqrt(V/n)
    ci = c(est.debias - qnorm(1-alpha/2)*se, est.debias + qnorm(1-alpha/2)*se)

    ############## Store Infos ###############
    est.plugin.vec[i.loading] = est.plugin
    est.debias.vec[i.loading] = est.debias
    se.vec[i.loading] = se
    ci.mat[i.loading, ] = ci
    proj.mat[, i.loading] = direction
  }

  obj <- list(est.plugin.vec = est.plugin.vec,
              est.debias.vec = est.debias.vec,
              se.vec         = se.vec,
              ci.mat         = ci.mat,
              proj.mat       = proj.mat)
  class(obj) = "LF"
  obj
}

