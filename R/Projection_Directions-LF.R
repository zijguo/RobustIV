# Constructs the projection direction with a fixed tuning parameter
#
# @description Constructs the projection direction, used for bias-correction, with a fixed tuning parameter
#
# @param X Design matrix, of dimension \eqn{n} x \eqn{p}
# @param loading Loading, of length \eqn{p}
# @param mu The dual tuning parameter used in the construction of the projection direction
# @param weight The weight vector of length \eqn{n}
# @param deriv.vec The first derivative vector of length \eqn{n}
# @return
# \item{proj}{The projection direction, of length \eqn{p}}
# @export
# @examples
# n <- 100
# p <- 400
# set.seed(1203)
# X <- matrix(sample(-2:2,n*p,replace = TRUE),nrow = n,ncol = p)
# resol <- 1.5
# step <- 3
# mu=sqrt(2.01*log(p)/n)*resol^{-(step-1)}
# weight = rep(1, n)
# deriv.vec = rep(1, n)
# Est <- Direction_fixedtuning(X,loading=c(1,rep(0,(p-1))),mu=mu,weight=weight,deriv.vec=deriv.vec)

Direction_fixedtuning <- function(X, loading, mu = NULL, weight = NULL, deriv.vec = NULL){
  pp <- ncol(X)
  n <- nrow(X)
  if(is.null(mu)){
    mu <- sqrt(2.01*log(pp)/n)
  }
  loading.norm <- sqrt(sum(loading^2))

  if (loading.norm == 0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  v <- CVXR::Variable(pp+1)
  obj <- 1/4*sum(((X%*%H%*%v)^2)*weight*deriv.vec)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob <- CVXR::Problem(CVXR::Minimize(obj))
  result <- CVXR::solve(prob)
  if(result$status=="optimal" || result$status == "unbounded"){
    opt.sol<-result$getValue(v)
    cvxr_status<-result$status
    direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  }else{
    direction <- numeric(0)
  }
  returnList <- list("proj"=direction)
  return(returnList)
}

# Searches for the best step size and computes the projection direction with the searched best step size
#
# @description Searches for the best step size and computes the projection direction with the searched best step size
#
# @param X Design matrix, of dimension \eqn{n} x \eqn{p}
# @param loading Loading, of length \eqn{p}
# @param weight The weight vector of length \eqn{n}
# @param deriv.vec The first derivative vector of length \eqn{n}
# @param resol The factor by which \code{mu} is increased/decreased to obtain the smallest \code{mu}
# such that the dual optimization problem for constructing the projection direction converges (default = 1.5)
# @param maxiter Maximum number of steps along which \code{mu} is increased/decreased to obtain the smallest \code{mu}
# such that the dual optimization problem for constructing the projection direction converges (default = 6)
#
# @return
# \item{proj}{The projection direction, of length \eqn{p}}
# \item{step}{The best step size}
# @export
#
# @examples
# n <- 100
# p <- 400
# set.seed(1203)
# X <- matrix(sample(-2:2,n*p,replace = TRUE),nrow = n,ncol = p)
# weight = rep(1, n)
# deriv.vec = rep(1, n)
# Est <- Direction_searchtuning(X,loading=c(1,rep(0,(p-1))), weight=weight, deriv.vec=deriv.vec)
Direction_searchtuning <- function(X, loading, weight = NULL, deriv.vec = NULL, resol = 1.5, maxiter = 6){
  pp <- ncol(X)
  n <- nrow(X)
  tryno <- 1
  opt.sol <- rep(0,pp+1)
  lamstop <- 0
  cvxr_status <- "optimal"
  mu <- sqrt(2.01*log(pp)/n)
  while (lamstop == 0 && tryno < maxiter){
    lastv <- opt.sol;
    lastresp <- cvxr_status;
    loading.norm <- sqrt(sum(loading^2))
    if (loading.norm == 0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    v <- CVXR::Variable(pp+1)
    obj <- 1/4*sum(((X%*%H%*%v)^2)*weight*deriv.vec)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))

    prob <- CVXR::Problem(CVXR::Minimize(obj))
    result <- CVXR::solve(prob)
    cvxr_status <- result$status
    if(tryno == 1){
      if(cvxr_status == "optimal"){
        incr = 0
        mu=mu/resol
        opt.sol <- result$getValue(v)
        temp.vec <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        initial.sd <- sqrt(sum(((X%*% temp.vec)^2)*weight*deriv.vec)/(n)^2)*loading.norm
        temp.sd <- initial.sd
      } else {
        incr <- 1
        mu <- mu*resol
      }
    } else {
      if(incr == 1){
        if(cvxr_status == "optimal"){
          opt.sol <- result$getValue(v)
          lamstop <- 1
        } else {
          mu <- mu*resol
        }
      } else {
        if(cvxr_status == "optimal" && temp.sd < 3*initial.sd){
          mu <- mu/resol
          opt.sol <- result$getValue(v)
          temp.vec <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)

          temp.sd <- sqrt(sum(((X%*% temp.vec)^2)*weight*deriv.vec)/(n)^2)*loading.norm
        } else {
          mu <- mu*resol
          opt.sol <- lastv
          lamstop <- 1
          tryno <- tryno-1
        }
      }
    }
    tryno = tryno + 1
  }
  direction <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step <- tryno-1
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}
