th_ICT <-
function(cont.var, time, status, predict.time, costs=NULL, R=NULL, method=c("normal", "empirical"), var.equal=FALSE, m=10, ci=TRUE, alpha=0.05, B=1000, range=3){
    ## checks
    if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
    if (!is.numeric(cont.var)) stop("'cont.var' must be numeric")
    if (length(cont.var)!=length(time)) stop("'cont.var' and 'time' should have the same length")
    if (length(cont.var)!=length(status)) stop("'cont.var' and 'status' should have the same length")
    if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
    if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
    if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
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
    result <- matrix(rep(NA, m*2), ncol=m, nrow=2)
    #result.AUC <- matrix(rep(NA, m*2), ncol=m, nrow=2)
    result.sens <- matrix(rep(NA, m*2), ncol=m, nrow=2)
    result.spec <- matrix(rep(NA, m*2), ncol=m, nrow=2)
    for (k in 1:m){
      imputed.k <- ExtractSingle(imputed.data.sets, index=k)
      status.new <- with(imputed.k$data, ifelse(impute.time>predict.time | (impute.time==predict.time & impute.event==0), 0, ifelse(impute.event==1, 1, NA)))
      # estimacio
      if (dim(table(status.new))==2 & all(table(status.new)>1)){
        out.aux <- with(df, minCostThresholdROC(cont.var, status.new, rho, costs2, method, var.equal, ci=TRUE, plot=FALSE, alpha=alpha, B=B))
        out <- c(out.aux$T$thres, out.aux$CI$se)
        cont.var.cat <- ifelse(cont.var>=out.aux$T$thres, ">=", "<")
        tab <- table(cont.var.cat, status.new)[2:1, 2:1]
        sens.est <- tab[1, 1]/(tab[1, 1]+tab[2, 1])
        sens.se <- sqrt(sens.est*(1-sens.est)/(tab[1, 1]+tab[2, 1]))
        spec.est <- tab[2, 2]/(tab[1, 2]+tab[2, 2])
        spec.se <- sqrt(spec.est*(1-spec.est)/(tab[1, 2]+tab[2, 2]))
        out.sens <- c(sens.est, sens.se)
        out.spec <- c(spec.est, spec.se)
        #rocobj <- roc(status.new, df$cont.var, quiet=TRUE)
      }else{
        out <- out.sens <- out.spec <- c(NA, NA)
      }
      result[1, k] <- out[1]
      result[2, k] <- out[2]
      result.sens[1, k] <- out.sens[1]
      result.sens[2, k] <- out.sens[2]
      result.spec[1, k] <- out.spec[1]
      result.spec[2, k] <- out.spec[2]
      #result[3, k] <- out[3]
      #result.AUC[1, k] <- as.numeric(rocobj$auc)
      #result.AUC[2, k] <- sqrt(var(rocobj))
    }
    ## threshold
    # m estimates --> rubin
    est <- mean(result[1, ], na.rm=TRUE)
    # amb estimador windsoritzat
    bp.est <- boxplot(result[1, ], range=range, plot=FALSE)
    if (length(bp.est$out)>0) {
      est <- winsor.mean(result[1, ], trim=0.25)
    }
    ## sens, spec
    est.sens <- mean(result.sens[1, ], na.rm=TRUE)
    est.spec <- mean(result.spec[1, ], na.rm=TRUE)
    # amb estimador windsoritzat
    bp.est.sens <- boxplot(result.sens[1, ], range=range, plot=FALSE)
    if (length(bp.est.sens$out)>0) {
      est.sens <- winsor.mean(result.sens[1, ], trim=0.25)
    }
    bp.est.spec <- boxplot(result.spec[1, ], range=range, plot=FALSE)
    if (length(bp.est.spec$out)>0) {
      est.spec <- winsor.mean(result.spec[1, ], trim=0.25)
    }
    if (ci){
      ## IC threshold
      var.w <- mean(result[2, ]^2, na.rm=TRUE)
      var.b <- var(result[1, ], na.rm=TRUE)
      # amb estimador windsoritzat
      bp.se <- boxplot(result[2, ], range=range, plot=FALSE)
      if (length(bp.se$out)>0){
        var.w <- winsor.mean(result[2, ]^2, trim=0.25)
      }
      if (length(bp.est$out)>0) {
        var.b <- winsor.var(result[1, ], trim=0.25)
      }
      var <- var.w+var.b+var.b/m
      IC <- mean(result[1, ], na.rm=TRUE)+c(-1, 1)*qnorm(1-alpha/2)*sqrt(var)
      ## IC sens, spec
      # sens
      var.w.sens <- mean(result.sens[2, ]^2, na.rm=TRUE)
      var.b.sens <- var(result.sens[1, ], na.rm=TRUE)
      # amb estimador windsoritzat
      bp.se.sens <- boxplot(result.sens[2, ], range=range, plot=FALSE)
      if (length(bp.se.sens$out)>0){
        var.w.sens <- winsor.mean(result.sens[2, ]^2, trim=0.25)
      }
      if (length(bp.est.sens$out)>0) {
        var.b.sens <- winsor.var(result.sens[1, ], trim=0.25)
      }
      var.sens <- var.w.sens+var.b.sens+var.b.sens/m
      IC.sens <- mean(result.sens[1, ], na.rm=TRUE)+c(-1, 1)*qnorm(1-alpha/2)*sqrt(var.sens)
      # spec
      var.w.spec <- mean(result.spec[2, ]^2, na.rm=TRUE)
      var.b.spec <- var(result.spec[1, ], na.rm=TRUE)
      # amb estimador windsoritzat
      bp.se.spec <- boxplot(result.spec[2, ], range=range, plot=FALSE)
      if (length(bp.se.spec$out)>0){
        var.w.spec <- winsor.mean(result.spec[2, ]^2, trim=0.25)
      }
      if (length(bp.est.spec$out)>0) {
        var.b.spec <- winsor.var(result.spec[1, ], trim=0.25)
      }
      var.spec <- var.w.spec+var.b.spec+var.b.spec/m
      IC.spec <- mean(result.spec[1, ], na.rm=TRUE)+c(-1, 1)*qnorm(1-alpha/2)*sqrt(var.spec)
      # output
      out <- list(T=list(thres=est, prev=rho, costs=costs2, R=(1-rho)/rho*(costs2[1, 2]-costs2[2, 1])/(costs2[1, 1]-costs2[2, 2]), method=method),
                  CI=list(lower=IC[1], upper=IC[2], se=sqrt(var), alpha=alpha, ci.method=ifelse(method=="normal", "delta", "boot")),
                  sens=list(est=est.sens, lower=IC.sens[1], upper=IC.sens[2]),
                  spec=list(est=est.spec, lower=IC.spec[1], upper=IC.spec[2]))
    }else{
      # output
      out <- list(T=list(thres=est, prev=rho, costs=costs2, R=(1-rho)/rho*(costs2[1, 2]-costs2[2, 1])/(costs2[1, 1]-costs2[2, 2]), method=method),
                  CI=NULL,
                  sens=list(est=est.sens, lower=NULL, upper=NULL),
                  spec=list(est=est.spec, lower=NULL, upper=NULL))
    }
    out$data <- data.frame(cont.var, time, status, statusNA=status.predict.time.NA)
    
    
    # # plot -> right now no plot for this method
    # if (plot){
    #   
    # }
    
    class(out) <- "th_ICT"
    return(out)
}


print.th_ICT <- function(x, ...){
  if (x$T$method == "normal"){
    cat("\nEstimate:")
    cat("\n  Threshold: ", x$T$thres)
    cat("\n")
    if (!is.null(x$CI)){
      # if(x$CI$ci.method == "delta"){
        cat("\nConfidence interval (delta method):")
        cat("\n  Lower Limit:", x$CI$lower)
        cat("\n  Upper Limit:", x$CI$upper)
        cat("\n")
      # }
      # if(x$CI$ci.method == "boot"){
      #   cat("\nConfidence intervals (bootstrap):")
      #   cat("\n  CI based on normal distribution: ", x$CI$low.norm, " - ", x$CI$up.norm)
      #   cat("\n  CI based on percentiles: ", x$CI$low.perc, " - ", x$CI$up.perc)
      #   cat("\n  Bootstrap resamples: ", x$CI$B)
      #   cat("\n")
      # }
    }
    cat("\nSensitivity: ", x$sens$est, "\n")
    
    cat("\nSpecificity: ", x$spec$est, "\n")
    
    
    cat("\nParameters used:")
    cat("\n  Disease prevalence:", x$T$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$T$costs)
    cat("\n  R:", x$T$R)
    cat("\n  Method:", x$T$method)
    if (!is.null(x$CI)){
      cat("\n  Significance Level: ", x$CI$alpha)
    }
    cat("\n")
  }
  if (x$T$method == "empirical"){
    cat("\nEstimate:")
    cat("\n  Threshold: ", x$T$thres)
    cat("\n  Minimum Cost: ", x$T$cost)
    if (!is.null(x$CI)){
      cat("\n")
      cat("\nConfidence interval (Rubin):")
      cat("\n  Lower Limit:", x$CI$lower)
      cat("\n  Upper Limit:", x$CI$upper)
      cat("\n")
    }
    cat("\nSensitivity: ", x$sens$est, "\n")
    
    cat("\nSpecificity: ", x$spec$est, "\n")
    
    cat("\n")
    cat("\nParameters used:")
    cat("\n  Disease prevalence:", x$T$prev)
    cat("\n  Costs (Ctp, Cfp, Ctn, Cfn):", x$T$costs)
    cat("\n  R:", x$T$R)
    cat("\n  Method:", x$T$method)
    if (!is.null(x$CI)){
      cat("\n  Significance Level:", x$CI$alpha)
    }
    cat("\n")
  }

}


