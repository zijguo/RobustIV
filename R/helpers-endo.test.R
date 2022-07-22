#' Summary of endotest
#'
#' @param object endotest object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.endotest<- function(object,...){
  endotest <- object
  if (typeof(endotest$VHat)=="list") {
    result<-matrix(NA, ncol=5, nrow=length(endotest$VHat))
    result <- data.frame(result)
    colnames(result)<-c("Estimated Cov","Test statistics","P value","Test result","Valid IVs")
    rownames(result)<-paste0("MaxClique",1:length(endotest$VHat))
    result[,1] <- unlist(endotest$Sigma12)
    result[,2] <- unlist(endotest$Q)
    result[,3] <- unlist(endotest$p.value)
    result[,4] <- ifelse(unlist(endotest$check),"H0 rejected","H0 not rejected")
    for (i in 1:length(endotest$VHat)) {
      result[i,5] <- paste(endotest$VHat[[i]], collapse = " ")
    }
    cat("Test result with significance level",endotest$alpha,"\n")
    cat(rep("_", 30), "\n")
    print(result,right=F)
  } else {
    cat("\nValid Instruments:", endotest$VHat, "\n");
    cat(rep("_", 30), "\n")
    cat("\nEstimated covariance:",endotest$Sigma12,"\n");
    cat("Test statistics Q = ",endotest$Q,"\n")
    cat("P-value = ",endotest$p.value,"\n")
    if (endotest$check) {
      cat("'H0 : Sigma12 = 0' is rejected at the significance level",endotest$alpha,".\n")
    } else {
      cat("'H0 : Sigma12 = 0' is not rejected at the significance level",endotest$alpha,".\n")
    }
  }

}
endo.SHat <- function(n, ITT_D, V.gamma, method='OLS', tuning.1st=NULL, tuning.2nd=NULL){
  pz = nrow(V.gamma)
  if(method=="OLS"){
    Tn1 = Tn2 = sqrt(log(n))
  }else{
    Tn1 = Tn2 = max(sqrt(2.01*log(pz)), sqrt(log(n)))
  }
  if(!is.null(tuning.1st)) Tn1 = tuning.1st
  if(!is.null(tuning.2nd)) Tn2 = tuning.1st
  ## First Stage
  SHat = (1:pz)[abs(ITT_D) > (Tn1 * sqrt(diag(V.gamma)/n))]
  ## First Stage
  SHat = (1:pz)[abs(ITT_D) > (Tn1 * sqrt(diag(V.gamma)/n))]

  if(length(SHat)==0){
    warning("First Thresholding Warning: IVs individually weak.
            TSHT with these IVs will give misleading CIs, SEs, and p-values.
            Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  return(SHat)
}
