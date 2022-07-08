AUC_ICT <-
function(cont.var, time, status, predict.time, m=10, ci=TRUE, alpha=0.05, range=3){
    ## checks
    if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
    if (!is.numeric(cont.var)) stop("'cont.var' must be numeric")
    if (length(cont.var)!=length(time)) stop("'cont.var' and 'time' should have the same length")
    if (length(cont.var)!=length(status)) stop("'cont.var' and 'status' should have the same length")
    if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
    if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
    if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
    # vector with missings
    status.predict.time.NA <- ifelse(time>predict.time | (time==predict.time & status==0), 0, ifelse(status==1, 1, NA))
    if (sum(is.na(status.predict.time.NA))==0) stop("There are no NAs in the status at 'predict.time'. It is not possible to apply ICT.")
    
    # imputation with InformativeCensoring
    df <- data.frame(id=1:length(cont.var), cont.var, time, status, status.predict.time.NA)
    df$arm <- 1
    df$arm <- factor(df$arm, levels=c("0", "1"))
    df$DCO.time <- Inf
    df$to.impute <- is.na(df$status.predict.time.NA)
    col.control <- col.headings(has.event="status", time="time", Id="id", arm="arm",
                                DCO.time="DCO.time", to.impute="to.impute")
    imputed.data.sets <- ScoreImpute(data=df, event.model=~cont.var,
                                     col.control=col.control, m=m,
                                     bootstrap.strata=df$arm,
                                     NN.control=NN.options(NN=10, w.censoring=0.2))
    # estimation for each imputated dataset
    #result <- matrix(rep(NA, m*2), ncol=m, nrow=2)
    result.AUC <- matrix(rep(NA, m*2), ncol=m, nrow=2)
    for (k in 1:m){
      imputed.k <- ExtractSingle(imputed.data.sets, index=k)
      status.new <- with(imputed.k$data, ifelse(impute.time>predict.time | (impute.time==predict.time & impute.event==0), 0, ifelse(impute.event==1, 1, NA)))
      # estimacio
      if (dim(table(status.new))==2 & all(table(status.new)>1)){
        #out.aux <- with(df, minCostThresholdROC(cont.var, status.new, rho, costs, method, var.equal, ci=TRUE, plot=FALSE, alpha=alpha, B=B))
        #out <- c(out.aux$T$thres, out.aux$CI$se)
        rocobj <- roc(status.new, df$cont.var, quiet=TRUE)
      }else{
        out <- c(NA, NA)
      }
      #result[1, k] <- out[1]
      #result[2, k] <- out[2]
      #result[3, k] <- out[3]
      result.AUC[1, k] <- as.numeric(rocobj$auc)
      result.AUC[2, k] <- sqrt(var(rocobj))
    }
    # m estimates --> rubin
    est <- mean(result.AUC[1, ], na.rm=TRUE)
    # amb estimador windsoritzat
    bp.est <- boxplot(result.AUC[1, ], range=range, plot=FALSE)
    if (length(bp.est$out)>0) {
      est <- winsor.mean(result.AUC[1, ], trim=0.25)
    }
    if (ci){
      # IC
      var.w <- mean(result.AUC[2, ]^2, na.rm=TRUE)
      var.b <- var(result.AUC[1, ], na.rm=TRUE)
      # amb estimador windsoritzat
      bp.se <- boxplot(result.AUC[2, ], range=range, plot=FALSE)
      if (length(bp.se$out)>0){
        var.w <- winsor.mean(result.AUC[2, ]^2, trim=0.25)
      }
      if (length(bp.est$out)>0) {
        var.b <- winsor.var(result.AUC[1, ], trim=0.25)
      }
      var <- var.w+var.b+var.b/m
      IC <- expit(logit(est)+c(-1, 1)*qnorm(1-alpha/2)*sqrt(var)/(est*(1-est)))
      out <- list(estimate=est, se=sqrt(var), CI=IC)
    }else{
      out <- list(estimate=est, se=NULL, CI=NULL)
    }
    out$data <- data.frame(cont.var, time, status, statusNA=status.predict.time.NA)
    
    
    # # plot -> right now no plot for this method
    # if (plot){
    #   
    # }
    
    class(out) <- "AUC_ICT"
    return(out)
}


print.AUC_ICT <- function(x, ...){
  cat("AUC: ", x$estimate)
  if(!is.null(x$se)){
    cat("\nStandard error: ", x$se)
    cat("\nConfidence interval: ", x$CI[1], " - ", x$CI[2])
  }
}

