AUC_ME <-
function(cont.var, time, status, predict.time, plot=FALSE, ci=TRUE, alpha=0.05, ...){
    ## checks
    if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
    if (!is.numeric(cont.var)) stop("'cont.var' must be numeric")
    if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
    if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
    if (length(cont.var)!=length(time)) stop("'cont.var' and 'time' should have the same length")
    if (length(cont.var)!=length(status)) stop("'cont.var' and 'status' should have the same length")
    if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
    # vector with missings
    status.predict.time.NA <- ifelse(time>predict.time | (time==predict.time & status==0), 0, ifelse(status==1, 1, NA))
    if (sum(is.na(status.predict.time.NA))==0) warning("There are no NAs in the status at 'predict.time'. The estimation proceeds without missing exclusion.")
    if (dim(table(status.predict.time.NA))!=2){
      warning("Missing exclusion leads to one single group. Returning NA as AUC")
      out <- NA
    }else{
      rocobj <- roc(status.predict.time.NA, cont.var, quiet=TRUE)
      auc <- as.numeric(rocobj$auc)
      if (ci){
        auc.var <- var(rocobj)
        se <- sqrt(auc.var)
        CI <- expit(logit(auc)+c(-1, 1)*qnorm(1-alpha/2)*se/(auc*(1-auc)))
      }else{
        CI <- NA
      }
      out <- list(AUC=auc, CI=CI)
      if (plot){
        plot(rocobj, ...)
      }
    }
    out$data <- data.frame(cont.var, time, status, statusNA=status.predict.time.NA)
    class(out) <- "AUC_ME"
    return(out)
    
}

print.AUC_ME <- function(x, ...){
  cat("AUC: ", x$AUC)
  cat("\nConfidence interval: ", x$CI[1], " - ",  x$CI[2])
}
