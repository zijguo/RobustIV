getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
}

relevant.funs <- function(intercept=TRUE, model=c("linear","logistic","logistic_alternative","probit")){
  model = match.arg(model)

  ### init step function ###
  if(model=="linear"){

    train.fun <- function(X, y, lambda=NULL){
      if(is.null(lambda)) lambda = "CV.min"
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
      if(intercept==FALSE) htheta = htheta[2:(p+1)]

      return(list(lasso.est = htheta))
    }

    cond_var.fun <- function(pred, y){
      n = length(y)
      sigma.sq = mean((y - pred)^2)
      return(rep(sigma.sq, n))
    }

  }else{

    train.fun <- function(X, y, lambda=NULL){
      if(is.null(lambda)) lambda = "CV.min"
      htheta <- if (lambda == "CV.min") {
        outLas <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                            intercept = intercept, standardize = T)
        as.vector(coef(outLas, s = outLas$lambda.min))
      } else if (lambda == "CV") {
        outLas <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                            intercept = intercept, standardize = T)
        as.vector(coef(outLas, s = outLas$lambda.1se))
      } else {
        outLas <- glmnet(X, y, family = "binomial", alpha = 1,
                         intercept = intercept, standardize = T)
        as.vector(coef(outLas, s = lambda))
      }
      if(intercept==FALSE) htheta = htheta[2:(p+1)]

      return(list(lasso.est = htheta))
    }

    cond_var.fun <- function(pred, y=NULL){
      cond_var = pred * (1 - pred)
      return(cond_var)
    }
  }

  if(model=="linear"){
    pred.fun = function(x) x
    deriv.fun = function(x) rep(1, length(x))
    weight.fun = function(x) rep(1, length(x))
  }else if(model=="logistic"){
    pred.fun = function(x) exp(x)/(1+exp(x))
    deriv.fun = function(x) exp(x)/(1+exp(x))^2
    weight.fun = function(x) (1+exp(x))^2/exp(x)
  }else if(model=="logistic_alternative"){
    pred.fun = function(x) exp(x)/(1+exp(x))
    deriv.fun = function(x) exp(x)/(1+exp(x))^2
    weight.fun = function(x) rep(1, length(x))
  }else if(model=="probit"){
    pred.fun = function(x) pnorm(x)
    deriv.fun = function(x) dnorm(x)
    weight.fun = function(x){
      out = c()
      out[abs(x)<5]=(fp(x)/f(x)/(1-f(x)))[abs(x)<5]
      out[abs(x)>5]=((abs(x)+sqrt(x^2+8/pi))/2)[abs(x)>5]
      return(out)
    }
  }

  return(list(train.fun    = train.fun,
              pred.fun     = pred.fun,
              deriv.fun    = deriv.fun,
              weight.fun   = weight.fun,
              cond_var.fun = cond_var.fun))
}

stars.pval <- function(p.value){
  unclass(symnum(p.value, corr = FALSE, na = FALSE,
                 cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")))
}
