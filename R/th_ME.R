th_ME <-
function(cont.var, time, status, predict.time, costs=NULL, R=NULL, method=c("normal", "empirical"), var.equal=FALSE, plot=FALSE, ci=TRUE, alpha=0.05, B=1000, ...){
    ## checks
    if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
    if (!is.numeric(cont.var)) stop("'cont.var' must be numeric")
    if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
    if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
    if (length(cont.var)!=length(time)) stop("'cont.var' and 'time' should have the same length")
    if (length(cont.var)!=length(status)) stop("'cont.var' and 'status' should have the same length")
    if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
    method <- match.arg(method)
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
    
    # vector with missings
    status.predict.time.NA <- ifelse(time>predict.time | (time==predict.time & status==0), 0, ifelse(status==1, 1, NA))
    if (sum(is.na(status.predict.time.NA))==0) warning("There are no NAs in the status at 'predict.time'. The estimation proceeds without missing exclusion.")
    out <- minCostThresholdROC(cont.var, status.predict.time.NA, rho, costs2, method, var.equal, plot, ci, alpha, B, ...)
    if (dim(table(status.predict.time.NA))!=2){
      warning("Missing exclusion leads to one single group. Returning NA as threshold")
    }
    out$data <- data.frame(cont.var, time, status, statusNA=status.predict.time.NA)
    return(out)
    
}
