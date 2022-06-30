#' Summary of endotest
#'
#' @param object endotest object
#' @param ...
#' @keywords internal
#' @return
#' @export
summary.endotest<- function(object,...){
  endotest <- object
  cat("\nValid Instruments:", endotest$VHat, "\n");
  cat(rep("_", 30), "\n")
  cat("Estimated covariance:",endotest$Sigma12,"\n");
  cat("Test statistics Q = ",endotest$Q,"\n")
  cat("P-value = ",endotest$p.value,"\n")
  if (endotest$check) {
    cat("'H0 : Sigma12 = 0' is rejected at the significance level",endotest$alpha,".\n")
  } else {
    cat("'H0 : Sigma12 = 0' is not rejected at the significance level",endotest$alpha,".\n")
  }

}
endo.SHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD) {
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


  # Constants
  n = nrow(WUMat);
  pz = length(ITT_Y)
  # First Stage
  Tn = max(sqrt(2.01*log(pz)), sqrt(log(n)/2))

  SHat = (1:pz)[(abs(ITT_D) >= (sqrt(SigmaSqD * colSums(WUMat^2) /n) * sqrt(Tn^2/n)))]

  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  return(SHat)
}
