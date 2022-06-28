#' @title Control-Function
#' @description Implement the control function method for estimation and inference of nonlinear treatment effects.
#' @param formula A formula describing the model to be fitted.
#' @param d1 The transformation of baseline treatment value.
#' @param d2 The transformation of target treatment value.
#'
#' @return
#'    \code{cf} returns an object of class "cf", which is a list containing the following components:
#'    \item{\code{coefficients}}{The estimate of the coefficients in the outcome model.}
#'    \item{\code{vcov}}{The estimated covariance matrix of coefficients.}
#'    \item{\code{CausalEffect}}{The causal effect of increasing the treatment from \code{d1} to \code{d2}.}
#'
#' @details For example, the formula \code{Y ~ D + I(D^2)+X|Z+I(Z^2)+X} describes the model where
#' \eqn{Y = \alpha_0 + D\beta_1 + D^2\beta_2 + X\phi + u}
#' and
#' \eqn{D = \gamma_0 + Z\gamma_1 + Z^2\gamma_2 + X\psi + v}.
#' Here, the outcome is \code{Y}, the endogenous variables are \code{D} and \code{I(D^2)}, the baseline covariates are \code{X}, and the instrument variables are \code{Z}. The formula environment follows
#' the formula environment in the ivreg function in the AER package. The linear term of the endogenous variable, for example, \code{D}, must be included in the first of the right side of the formula.
#' If one of \code{d1} or \code{d2} is \code{NULL}, \code{CausalEffect} returns \code{NULL}.
#' @export
#'
#'
#' @importFrom Formula as.Formula
#'
#' @examples
#' \dontrun{
#' Y <- mroz[,"lwage"]
#' D <- mroz[,"educ"]
#' Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc")])
#' X <- as.matrix(mroz[,c("exper","expersq","age")])
#' cf.model <- cf(Y~D+I(D^2)+X|Z+I(Z^2)+X,d1 = c(median(D),median(D)^2),
#' d2 =c(median(D)+1,(median(D)+1)^2) )
#' summary(cf.model)
#' }
#' @references {
#' Guo, Z. and D. S. Small (2016), Control function instrumental variable estimation of nonlinear causal effect models, \emph{The Journal of Machine Learning Research} 17(1), 3448–3482. \cr
#' }
#'
#'


cf <- function(formula,d1 = NULL,d2 = NULL){
  if(!inherits(formula,"formula")) {
    stop("method is only for formula objects!")
  }
  # code gratefully lifted from ivreg() (package AER) and ivmodelFormula (package ivmodel).

  mf = match.call()
  m <- match(c("formula"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  formula <- as.Formula(formula)

  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in%
              1:2)
  has_dot <- function(formula) inherits(try(terms(formula),silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2)) {
      formula <- as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))
    }
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf, "numeric");
  n <- length(Y)
  Y = matrix(as.numeric(Y),length(Y),1)
  mt <- terms(formula)
  mtX <- terms(formula, rhs = 1)
  X <- model.matrix(mtX, mf)

  mtZ <- delete.response(terms(formula, rhs = 2))
  Z <- model.matrix(mtZ, mf)

  if("(Intercept)" %in% colnames(X)) {
    intercept=TRUE
    X = X[,!(colnames(X) %in% "(Intercept)"),drop=FALSE]
    Z = Z[,!(colnames(Z) %in% "(Intercept)"),drop=FALSE]
    if(dim(Z)[2] < 1) stop("There aren't any instruments!")
  } else{
    intercept=FALSE
  }

  # Parse X and Z into D, X, and Z

  whichD = !(colnames(X) %in% colnames(Z))
  D = X[,whichD,drop=FALSE]
  if (sum(!whichD) == 0) {
    if (intercept) {
      first.model <- lm(D~Z)
      e1 <- first.model$residuals[,1]
      second.model <- lm(Y~D+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef)[2:length(cf.coef)] = colnames(D)
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov)[2:length(cf.coef)] = colnames(cf.vcov)[2:length(cf.coef)] = colnames(D)
      if (!is.null(d1)&!is.null(d2)) {
        stopifnot(ncol(D)==length(d1))
        stopifnot(ncol(D)==length(d2))
        CausalEffect <- (d2-d1)%*%cf.coef[-1]
      } else {
        CausalEffect = NULL
      }

    } else {
      first.model <- lm(D~0+Z)
      e1 <- first.model$residuals[,1]
      second.model <- lm(Y~0+D+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef) = colnames(D)
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov) = colnames(cf.vcov) = colnames(D)
      if (!is.null(d1)&!is.null(d2)) {
      stopifnot(ncol(D)==length(d1))
      stopifnot(ncol(D)==length(d2))
      CausalEffect <- (d2-d1)%*%cf.coef
      } else {
        CausalEffect = NULL
      }
    }


  } else {
    X = X[,!whichD,drop=FALSE]
    whichZ = !(colnames(Z) %in% colnames(X))
    Z = Z[,whichZ,drop=FALSE]
    if (intercept) {
      first.model <- lm(D~Z+X)
      e1 <- first.model$residuals[,1]
      second.model <- lm(Y~D+X+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef)[2:length(cf.coef)] = c(colnames(D),colnames(X))
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov)[2:length(cf.coef)] = colnames(cf.vcov)[2:length(cf.coef)] = c(colnames(D),colnames(X))
      if (!is.null(d1)&!is.null(d2)) {
        stopifnot(ncol(D)==length(d1))
        stopifnot(ncol(D)==length(d2))
        CausalEffect <- (d2-d1)%*%cf.coef[-c(1,seq(length(cf.coef)-ncol(X)+1,length =ncol(X)))]
      } else{
        CausalEffect <- NULL
      }

    } else {
      first.model <- lm(D~0+Z+X)
      e1 <- first.model$residuals[,1]
      second.model <- lm(Y~0+D+X+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef) = c(colnames(D),colnames(X))
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov) = colnames(cf.vcov) = c(colnames(D),colnames(X))
      if (!is.null(d1)&!is.null(d2)) {
        stopifnot(ncol(D)==length(d1))
        stopifnot(ncol(D)==length(d2))
        CausalEffect <- (d2-d1)%*%cf.coef[(1:length(cf.coef))[-seq(length(cf.coef)-ncol(X)+1,length =ncol(X))]]
      } else{
        CausalEffect <- NULL
      }

    }

  }
  out <- list(coefficients=cf.coef,vcov =cf.vcov,CausalEffect=CausalEffect,n=n,d1 = d1,d2 = d2)
  class(out) = 'cf'
  return(out)
}

#' Summary of cf
#'
#' @param object cf object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.cf<- function(object,...){
  cf <- object
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of Control Function Estimators:\n\n")
  coeff <- cf$coefficients
  std <- sqrt(diag(cf$vcov))
  t.value <- abs(coeff/std)
  pr.t <- 1-pt(t.value,df = (cf$n)-1)
  cmat <- cbind(coeff,std,t.value,pr.t)
  colnames(cmat) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  printCoefmat(cmat, digits = max(3L, getOption("digits") - 3L))
  if (!is.null(cf$d2)&!is.null(cf$d1)) {
    cat("Causal effect from D =",cf$d1[1],"to",cf$d2[1],": ",cf$CausalEffect,"\n" )
  }
}



#' @title Prestest estimator
#' @description This function implements the pretest approach  for estimation and inference of nonlinear treatment effects.
#'
#' @param formula A formula describing the model to be fitted.
#' @param alpha The significant level. (default = \code{0.05})
#'
#' @details For example, the formula \code{Y ~ D + I(D^2)+X|Z+I(Z^2)+X} describes the model where
#' \eqn{Y = \alpha_0 + D\beta_1 + D^2\beta_2 + X\phi + u}
#' and
#' \eqn{D = \gamma_0 + Z\gamma_1 + Z^2\gamma_2+X\psi + v}.
#' Here, the outcome is \code{Y}, the endogenous variables are \code{D} and \code{I(D^2)}, the baseline covariates are \code{X}, and the instrument variables are \code{Z}. The formula environment follows
#' the formula environment in the ivreg function in the AER package. The linear term of the endogenous variable, for example, \code{D}, must be included in the first of the right side of the formula.

#' @return
#'    \code{pretest} returns an object of class "pretest", which is a list containing the following components:
#'    \item{\code{coefficients}}{The estimate of the coefficients in the outcome model.}
#'    \item{\code{vcov}}{The estimated covariance matrix of coefficients.}
#'    \item{\code{Hausman.stat}}{The Hausman test statistic used to test the validity of control function.}
#'    \item{\code{p.value}}{The p-value of Hausman test.}
#'    \item{\code{cf.check}}{Indicator for whether the control function is valid or not.}
#' @export
#'
#'
#' @examples
#' \dontrun{
#' Y <- mroz[,"lwage"]
#' D <- mroz[,"educ"]
#' Z <- as.matrix(mroz[,c("motheduc","fatheduc","huseduc")])
#' X <- as.matrix(mroz[,c("exper","expersq","age")])
#' pretest.model <- pretest(Y~D+I(D^2)+X|Z+I(Z^2)+X)
#' summary(pretest.model)
#' }
#' @references {
#' Guo, Z. and D. S. Small (2016), Control function instrumental variable estimation of nonlinear causal effect models, \emph{The Journal of Machine Learning Research} 17(1), 3448–3482. \cr
#' }
#'
pretest <- function(formula,alpha = 0.05){


  tsls.model <- AER::ivreg(formula)
  n <- tsls.model$n
  iv.coef <- coef(tsls.model)
  iv.vcov <- vcov(tsls.model)

  cf.model <- cf(formula)

  cf.coef <- cf.model$coefficients
  cf.vcov <- cf.model$vcov

  diff <- t(iv.coef-cf.coef)%*%solve(iv.vcov-cf.vcov)%*%(iv.coef-cf.coef)
  prob.larger.than.diff <- pchisq(diff, df = 1, lower.tail = FALSE)
  prob.larger.than.diff
  if (prob.larger.than.diff>alpha) {
    cf.check = TRUE
    pretest.coef <- cf.coef
    pretest.vcov <- cf.vcov
  } else{

    cf.check = FALSE
    pretest.coef <- iv.coef
    pretest.vcov <- iv.vcov
  }

  pretest.val <- list(
    coefficients = pretest.coef,
    vcov = pretest.vcov,
    Hausman.stat = diff,
    p.value = prob.larger.than.diff,
    alpha = alpha,
    cf.check = cf.check,
    n = n
  )
  class(pretest.val) = 'pretest'
  return(pretest.val)
}
#' Summary of pretest
#'
#' @param object pretest object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.pretest<- function(object,...){
  pretest <- object
  cat(rep("_", 30), "\n")
  cat("\nHausman Statistic : ",pretest$Hausman.stat,"\n")
  cat("\nP value = ",pretest$p.value,"\n")
  if (pretest$cf.check) {
    cat("\nH0 : Augmented instrumental variables from Control function are valid, is not rejected.","\n")
    cat("\nLevel",pretest$alpha, "Pretest estimator is Control function estimator.","\n")
  } else{
    cat("\nH0 : Augmented instrumental variables from Control function are valid, is rejected.","\n")
    cat("\nLevel",pretest$alpha, "Pretest estimator is Two Stage Least Square estimator.","\n")
  }
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of Pretest Estimators:\n\n")
  coeff <- pretest$coefficients
  std <- sqrt(diag(pretest$vcov))
  t.value <- abs(coeff/std)
  pr.t <- 1-pt(t.value,df = (pretest$n)-1)
  cmat <- cbind(coeff,std,t.value,pr.t)
  colnames(cmat) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  printCoefmat(cmat, digits = max(3L, getOption("digits") - 3L))
}
