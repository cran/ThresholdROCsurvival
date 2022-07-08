th_survivalROC <-
function(cont.var, time, status, predict.time, costs=NULL, R=NULL, lambda=0.05, plot=FALSE, ci=FALSE, alpha=0.05, B=1000, ...){
    ## checks
    if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
    if (!is.numeric(cont.var)) stop("'cont.var' must be numeric")
    if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
    if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
    if (length(cont.var)!=length(time)) stop("'cont.var' and 'time' should have the same length")
    if (length(cont.var)!=length(status)) stop("'cont.var' and 'status' should have the same length")
    if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
    # rho estimation
    # KM
    fit <- survfit(Surv(time, status) ~ 1)
    rho <- 1-summary(fit, times=predict.time)$surv # rho=1-S(predict.time)
    
    # costs/R
    if (is.null(costs) & is.null(R)){
      costs2 <- matrix(c(0, 0, 1, (1-rho)/rho), 2, 2, byrow=TRUE)
    }else if (is.null(costs) & !is.null(R)){
      if (!is.numeric(R) | length(R)!=1){
        stop("R must be a single number")
      }else{
        costs2 <- matrix(c(0, 0, 1, (1-rho)/(R*rho)), 2, 2, byrow=TRUE)
      }
    }else if (!is.null(costs) & is.null(R)){
      costs2 <- costs
      if (!is.matrix(costs)){
        stop("'costs' must be a matrix")
      }
      if (dim(costs)[1] != 2 | dim(costs)[2] != 2){
        stop("'costs' must be a 2x2 matrix")
      }
    }else if (!is.null(costs) & !is.null(R)){
      stop("either 'costs' or 'R' (or both) must be NULL")
    }
    
    out <- minCostSurvivalROC(cont.var, time, status, predict.time, rho, costs2, lambda, plot, ci, alpha=alpha, B=B, ...)
    
    if (ci){
      out2 <- list(T=list(thres=as.numeric(out$cut.est$optimum[1]), prev=rho, costs=costs2, R=(1-rho)/rho*(costs2[1, 2]-costs2[2, 1])/(costs2[1, 1]-costs2[2, 2])),
                  CI=list(lower.norm=out$CI[1, 1], upper.norm=out$CI[1, 2], se=out$se, lower.perc=out$CI[2, 1], upper.perc=out$CI[2, 2], alpha=alpha,
                          ci.method="boot", B=B))
    }else{
      out2 <- list(T=list(thres=as.numeric(out$cut.est$optimum[1]), prev=rho, costs=costs2, R=(1-rho)/rho*(costs2[1, 2]-costs2[2, 1])/(costs2[1, 1]-costs2[2, 2])),
                   CI=NULL)
    }
    
    class(out2) <- "th_survivalROC"
    return(out2)
}


print.th_survivalROC <- function(x, ...){
  cat("\nEstimate:")
  cat("\n  Threshold: ", x$T$thres)
  cat("\n")
  if (!is.null(x$CI)){
    cat("\nConfidence intervals (bootstrap):")
    cat("\n  CI based on normal distribution: ", x$CI$lower.norm, " - ", x$CI$upper.norm)
    cat("\n  CI based on percentiles: ", x$CI$lower.perc, " - ", x$CI$upper.perc)
    cat("\n  Bootstrap resamples: ", x$CI$B)
    cat("\n")
  }
  cat("\nParameters used:")
  cat("\n  Disease prevalence:", x$T$prev)
  cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$T$costs)
  cat("\n  R:", x$T$R)
  if (!is.null(x$CI)){
    cat("\n  Significance Level: ", x$CI$alpha)
  }
  cat("\n")
}
