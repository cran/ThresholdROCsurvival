
diagnostic_assessment_binary <- function(binary.var, time, status, predict.time, method=c("USE", "ICT"), index=c("all", "sens", "spec"), m=10, ci=TRUE, alpha=0.05, range=3){
  # checks
  method <- match.arg(method)
  index <- match.arg(index)
  if (!is.factor(binary.var)) stop("'binary.var' must be a factor with two levels: - and +")
  if (nlevels(binary.var)!=2) stop("'binary.var' must be a factor with two levels: - and +")
  if (!all(c("+", "-") %in% levels(binary.var))) stop("'binary.var' must be a factor with two levels: - and +")
  if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
  if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
  if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
  if (length(binary.var)!=length(time)) stop("'binary.var' and 'time' should have the same length")
  if (length(binary.var)!=length(status)) stop("'binary.var' and 'status' should have the same length")
  if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
  if (!is.logical(ci)) stop("'ci' should be TRUE or FALSE")
  if (!is.numeric(m) | length(m)!=1) stop("'m' should be numerical, of length 1")
  if (!is.numeric(alpha) | length(alpha)!=1) stop("'alpha' should be numerical, of length 1")
  if (alpha>1 | alpha<0) stop("'alpha' should be between 0 and 1")
  if (!is.numeric(range) | length(range)!=1) stop("'range' should be numerical, of length 1")
  if (range<0) stop("'range' should be positive")
  if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' should be numerical, of length 1")
  if (is.null(method)) stop("'method' cannot be NULL")
  if (is.null(index)) stop("'index' cannot be NULL")
  
  # status amb missings
  status.predict.time.NA <- ifelse(time>predict.time | (time==predict.time & status==0), 0, ifelse(status==1, 1, NA))
  if (sum(is.na(status.predict.time.NA))==0) warning("There are no NAs in the status at 'predict.time'.")
  
  # output
  out <- list()
  
  ##### ME
  if (method=="USE"){
    # checks
    if (sum(is.na(status.predict.time.NA))==0) warning("The estimation proceeds without unknown status exclusion.")
    if (dim(table(status.predict.time.NA))!=2) stop("USE leads to one single group. Estimation cannot proceed.")
    ## sens
    if ("all" %in% index | "sens" %in% index){
      tab <- table(binary.var, status.predict.time.NA)
      # taula 2x2 ok
      tab <- tab[c("+", "-"), c("1", "0")]
      sens.num <- tab[1, 1]
      sens.den <- tab[1, 1]+tab[2, 1]
      sens <- sens.num/sens.den
      if (ci){
        if (sens.den*sens>=5 & sens.den*(1-sens)>=5 & sens.den>=30){ # np>=5 i n(1-p)>=5 i n>=30
          result.sens <- prop.test(sens.num, sens.den, conf.level=1-alpha)
          out$sens <- c(est=sens, LL=result.sens$conf.int[1], UL=result.sens$conf.int[2])
        }else{
          result.sens <- binom.test(sens.num, sens.den, conf.level=1-alpha)
          out$sens <- c(est=sens, LL=result.sens$conf.int[1], UL=result.sens$conf.int[2])
        }
      }else{
        out$sens <- c(est=sens, LL=NA, UL=NA)
      }
    }
    ## spec
    if ("all" %in% index | "spec" %in% index){
      if (!exists("tab")){
        tab <- table(binary.var, status.predict.time.NA)
        # taula 2x2 ok
        tab <- tab[c("+", "-"), c("1", "0")]
      }
      spec.num <- tab[2, 2]
      spec.den <- tab[1, 2]+tab[2, 2]
      spec <- spec.num/spec.den
      if (ci){
        if (spec.den*spec>=5 & spec.den*(1-spec)>=5 & spec.den>=30){ # np>=5 i n(1-p)>=5 i n>=30
          result.spec <- prop.test(spec.num, spec.den, conf.level=1-alpha)
          out$spec <- c(est=spec, LL=result.spec$conf.int[1], UL=result.spec$conf.int[2])
        }else{
          result.spec <- binom.test(spec.num, spec.den, conf.level=1-alpha)
          out$spec <- c(est=spec, LL=result.spec$conf.int[1], UL=result.spec$conf.int[2])
        }
     }else{
        out$spec <- c(est=spec, LL=NA, UL=NA)
      }
    }
  ##### ICT
  }else if (method=="ICT"){
    # checks
    if (sum(is.na(status.predict.time.NA))==0) stop("ICT cannot be applied.")

    # imputation with InformativeCensoring
    df <- data.frame(id=1:length(binary.var), binary.var, time, status, status.predict.time.NA)
    df$arm <- 1
    df$arm <- factor(df$arm, levels=c("0", "1"))
    df$DCO.time <- Inf
    df$to.impute <- is.na(df$status.predict.time.NA)
    col.control <- col.headings(has.event="status", time="time", Id="id", arm="arm",
                                DCO.time="DCO.time", to.impute="to.impute")
    imputed.data.sets <- ScoreImpute(data=df, event.model=~binary.var,
                                     col.control=col.control, m=m,
                                     bootstrap.strata=df$arm,
                                     NN.control=NN.options(NN=10, w.censoring=0.2))
    # estimation for each imputated dataset
    if (!any(class(imputed.data.sets)=="ScoreImputedSet")){
      stop("Some problem has occurred in the imputation process.")
    }else{
      if ("all" %in% index | "sens" %in% index) result.sens <- matrix(rep(NA, m*2), ncol=m, nrow=2)
      if ("all" %in% index | "spec" %in% index) result.spec <- matrix(rep(NA, m*2), ncol=m, nrow=2)
      # 1a fila: estimate
      # 2a fila: se
      for (k in 1:m){
        imputed.k <- ExtractSingle(imputed.data.sets, index=k)
        status.new <- with(imputed.k$data, ifelse(impute.time>predict.time | (impute.time==predict.time & impute.event==0), 0, ifelse(impute.event==1, 1, NA)))
        # estimates
        if (dim(table(status.new))==2 & all(table(status.new)>1)){
          ### sens
          if ("all" %in% index | "sens" %in% index){
            tab <- table(binary.var, status.new)
            if (ncol(tab)==2 & nrow(tab)==2){
              tab <- tab[c("+", "-"), c("1", "0")]
              # sens
              sens.est <- tab[1, 1]/(tab[1, 1]+tab[2, 1])
              result.sens[1, k] <- sens.est
              result.sens[2, k] <- sqrt(sens.est*(1-sens.est)/(tab[1, 1]+tab[2, 1]))
            }
          }
          if ("all" %in% index | "spec" %in% index){
            if (!exists("tab")){
              tab <- table(binary.var, status.new)
            }
            if (ncol(tab)==2 & nrow(tab)==2){
              tab <- tab[c("+", "-"), c("1", "0")]
              # spec
              spec.est <- tab[2, 2]/(tab[1, 2]+tab[2, 2])
              result.spec[1, k] <- spec.est
              result.spec[2, k] <- sqrt(spec.est*(1-spec.est)/(tab[1, 2]+tab[2, 2]))
            }
          }
        }
      }
      # rubin
      if ("all" %in% index | "sens" %in% index) rubin.sens <- rubin_outliers(result.sens, range=range)
      if ("all" %in% index | "spec" %in% index) rubin.spec <- rubin_outliers(result.spec, range=range)
      # results
      ### sens
      if ("all" %in% index | "sens" %in% index){
        if (ci){
          IC.sens <- rubin.sens$est+c(-1, 1)*qnorm(1-alpha/2)*sqrt(rubin.sens$var)
          out$sens <- c(est=rubin.sens$est, LL=max(0, IC.sens[1]), UL=min(1, IC.sens[2]))
        }else{
          out$sens <- c(est=rubin.sens$est, LL=NA, UL=NA)
        }
      }
      ### spec
      if ("all" %in% index | "spec" %in% index){
        if (ci){
          IC.spec <- rubin.spec$est+c(-1, 1)*qnorm(1-alpha/2)*sqrt(rubin.spec$var)
          out$spec <- c(est=rubin.spec$est, LL=max(0, IC.spec[1]), UL=min(1, IC.spec[2]))
        }else{
          out$spec <- c(est=rubin.spec$est, LL=NA, UL=NA)
        }
      }
    }
  }
  
  # output
  out$method <- method
  out$alpha <- alpha
  out$data <- data.frame(binary.var, time, status, statusNA=status.predict.time.NA)
  class(out) <- "diagnostic_assessment"
  return(out)
}




diagnostic_assessment_continuous <- function(cont.var, time, status, predict.time, method=c("USE", "ICT", "survivalROC"), index=c("all", "AUC", "threshold", "sens", "spec"), costs=NULL, R=NULL, method.thres=c("normal", "empirical"), var.equal=FALSE, lambda=0.05, m=10, ci=TRUE, plot=FALSE, alpha=0.05, B=1000, range=3, ...){
  # checks
  method <- match.arg(method)
  method.thres <- match.arg(method.thres)
  if (predict.time>max(time)) stop("'predict.time' must be lower or equal to the greatest value in 'time'")
  if (!is.numeric(cont.var)) stop("'cont.var' must be numeric")
  if (!is.numeric(predict.time) | length(predict.time)!=1) stop("'predict.time' must be a single number")
  if (!all((names(table(status)) %in% c("0", "1")))) stop("'status' should only contain 0s and 1s")
  if (length(cont.var)!=length(time)) stop("'cont.var' and 'time' should have the same length")
  if (length(cont.var)!=length(status)) stop("'cont.var' and 'status' should have the same length")
  if (length(time)!=length(status)) stop("'time' and 'status' should have the same length")
  if (!is.logical(var.equal)) stop("'var.equal' should be TRUE or FALSE")
  if (!is.logical(ci)) stop("'ci' should be TRUE or FALSE")
  if (!is.logical(plot)) stop("'plot' should be TRUE or FALSE")
  if (!is.numeric(m) | length(m)!=1) stop("'m' should be numerical, of length 1")
  if (!is.numeric(alpha) | length(alpha)!=1) stop("'alpha' should be numerical, of length 1")
  if (alpha>1 | alpha<0) stop("'alpha' should be between 0 and 1")
  if (!is.numeric(range) | length(range)!=1) stop("'range' should be numerical, of length 1")
  if (range<0) stop("'range' should be positive")
  if (is.null(method)) stop("'method' cannot be NULL")
  if (is.null(index)) stop("'index' cannot be NULL")
  if (!all(index %in% c("all", "AUC", "threshold", "sens", "spec"))) stop("Invalid value in 'index'")
  if (("sens" %in% index | "spec" %in% index) & !("threshold" %in% index)) stop("'sens' and 'spec' require threshold estimation")
  
  # status amb missings
  status.predict.time.NA <- ifelse(time>predict.time | (time==predict.time & status==0), 0, ifelse(status==1, 1, NA))
  if (sum(is.na(status.predict.time.NA))==0) warning("There are no NAs in the status at 'predict.time'.")

  # output
  out <- list()
  
  if ("all" %in% index | "threshold" %in% index){
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
  }
  
  ##### ME
  if (method=="USE"){
    # checks
    if (sum(is.na(status.predict.time.NA))==0) warning("The estimation proceeds without unknown status exclusion.")
    if (dim(table(status.predict.time.NA))!=2) stop("USE leads to one single group. Estimation cannot proceed.")
    ## AUC
    if ("all" %in% index | "AUC" %in% index){
      rocobj <- roc(status.predict.time.NA, cont.var, quiet=TRUE)
      auc <- as.numeric(rocobj$auc)
      if (ci){
        auc.var <- var(rocobj)
        se <- sqrt(auc.var)
        CI <- expit(logit(auc)+c(-1, 1)*qnorm(1-alpha/2)*se/(auc*(1-auc)))
      }else{
        CI <- NA
      }
      out$AUC <- c(est=auc, LL=CI[1], UL=CI[2])
      if (plot){
        plot(rocobj, legacy.axes=TRUE, ...)
      }
    }
    ## threshold
    if ("all" %in% index | "threshold" %in% index){
      # estimation
      result <- suppressWarnings(minCostThresholdROC(cont.var, status.predict.time.NA, rho, costs2, method.thres, var.equal, plot, ci, alpha, B, ...))
      if (ci){
        out$threshold <- c(est=result$T$thres, LL=result$CI$lower, UL=result$CI$upper)
      }else{
        out$threshold <- c(est=result$T$thres, LL=NA, UL=NA)
      }
    }
    ## sens
    if ("all" %in% index | "sens" %in% index){
      X_cat <- ifelse(cont.var>=out$threshold[1], ">=", "<")
      if (dim(table(X_cat))!=2) stop("Dichotomized 'cont.var' does not result in a binary variable")
      tab <- table(X_cat, status.predict.time.NA)
      # taula 2x2
      tab <- tab[c(">=", "<"), c("1", "0")]
      sens.num <- tab[1, 1]
      sens.den <- tab[1, 1]+tab[2, 1]
      sens <- sens.num/sens.den
      if (ci){
        if (sens.den*sens>=5 & sens.den*(1-sens)>=5 & sens.den>=30){ # np>=5 i n(1-p)>=5 i n>=30
          result.sens <- prop.test(sens.num, sens.den, conf.level=1-alpha)
          out$sens <- c(est=sens, LL=result.sens$conf.int[1], UL=result.sens$conf.int[2])
        }else{
          result.sens <- binom.test(sens.num, sens.den, conf.level=1-alpha)
          out$sens <- c(est=sens, LL=result.sens$conf.int[1], UL=result.sens$conf.int[2])
        }
      }else{
        out$sens <- c(est=sens, LL=NA, UL=NA)
      }
    }
    ## spec
    if ("all" %in% index | "spec" %in% index){
      if (!exists("X_cat")){
        X_cat <- ifelse(cont.var>=out$threshold[1], ">=", "<")
        if (dim(table(X_cat))!=2) stop("Dichotomized 'cont.var' does not result in a binary variable")
        tab <- table(X_cat, status.predict.time.NA)
        # taula 2x2
        tab <- tab[c(">=", "<"), c("1", "0")]
      }
      spec.num <- tab[2, 2]
      spec.den <- tab[1, 2]+tab[2, 2]
      spec <- spec.num/spec.den
      if (ci){
        if (spec.den*spec>=5 & spec.den*(1-spec)>=5 & spec.den>=30){ # np>=5 i n(1-p)>=5 i n>=30
          result.spec <- prop.test(spec.num, spec.den, conf.level=1-alpha)
          out$spec <- c(est=spec, LL=result.spec$conf.int[1], UL=result.spec$conf.int[2])
        }else{
          result.spec <- binom.test(spec.num, spec.den, conf.level=1-alpha)
          out$spec <- c(est=spec, LL=result.spec$conf.int[1], UL=result.spec$conf.int[2])
        }
      }else{
        out$spec <- c(est=spec, LL=NA, UL=NA)
      }
    }
  ##### ICT
  }else if (method=="ICT"){
    # checks
    if (sum(is.na(status.predict.time.NA))==0) stop("ICT cannot be applied.")
    if (plot) warning("There are no plots implemented for the ICT method.")
    
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
    if (!any(class(imputed.data.sets)=="ScoreImputedSet")){
      stop("Some problem has occurred in the imputation process.")
    }else{
      if ("all" %in% index | "AUC" %in% index) result.AUC <- matrix(rep(NA, m*2), ncol=m, nrow=2)
      if ("all" %in% index | "threshold" %in% index) result.th <- matrix(rep(NA, m*2), ncol=m, nrow=2)
      if ("all" %in% index | "sens" %in% index) result.sens <- matrix(rep(NA, m*2), ncol=m, nrow=2)
      if ("all" %in% index | "spec" %in% index) result.spec <- matrix(rep(NA, m*2), ncol=m, nrow=2)
      # 1a fila: estimate
      # 2a fila: se
      for (k in 1:m){
        imputed.k <- ExtractSingle(imputed.data.sets, index=k)
        status.new <- with(imputed.k$data, ifelse(impute.time>predict.time | (impute.time==predict.time & impute.event==0), 0, ifelse(impute.event==1, 1, NA)))
        # estimates
        if (dim(table(status.new))==2 & all(table(status.new)>1)){
          ### AUC
          if ("all" %in% index | "AUC" %in% index){
            rocobj <- roc(status.new, cont.var, quiet=TRUE)
            result.AUC[1, k] <- as.numeric(rocobj$auc)
            result.AUC[2, k] <- sqrt(var(rocobj))
          }
          ### th
          if ("all" %in% index | "threshold" %in% index){
            result <- suppressWarnings(minCostThresholdROC(cont.var, status.new, rho, costs2, method.thres, var.equal, plot=FALSE, ci=TRUE, alpha, B))
            result.th[1, k] <- result$T$thres
            result.th[2, k] <- result$CI$se
          }
          ### sens
          if ("all" %in% index | "sens" %in% index){
            X_cat <- ifelse(cont.var>=result.th[1, k], ">=", "<")
            tab <- table(X_cat, status.new)
            if (ncol(tab)==2 & nrow(tab)==2){
              tab <- tab[c(">=", "<"), c("1", "0")]
              # sens
              sens.est <- tab[1, 1]/(tab[1, 1]+tab[2, 1])
              result.sens[1, k] <- sens.est
              result.sens[2, k] <- sqrt(sens.est*(1-sens.est)/(tab[1, 1]+tab[2, 1]))
            }
          }
          if ("all" %in% index | "spec" %in% index){
            if (!exists("X_cat")){
              X_cat <- ifelse(cont.var>=result.th[1, k], ">=", "<")
              tab <- table(X_cat, status.new)
            }
            if (ncol(tab)==2 & nrow(tab)==2){
              tab <- tab[c(">=", "<"), c("1", "0")]
              # spec
              spec.est <- tab[2, 2]/(tab[1, 2]+tab[2, 2])
              result.spec[1, k] <- spec.est
              result.spec[2, k] <- sqrt(spec.est*(1-spec.est)/(tab[1, 2]+tab[2, 2]))
            }
          }
        }
      }
      # rubin
      if ("all" %in% index | "AUC" %in% index) rubin.AUC <- rubin_outliers(result.AUC, range=range)
      if ("all" %in% index | "threshold" %in% index) rubin.th <- rubin_outliers(result.th, range=range)
      if ("all" %in% index | "sens" %in% index) rubin.sens <- rubin_outliers(result.sens, range=range)
      if ("all" %in% index | "spec" %in% index) rubin.spec <- rubin_outliers(result.spec, range=range)
      # results
      ### AUC
      if ("all" %in% index | "AUC" %in% index){
        if (ci){
          IC.auc <- expit(logit(rubin.AUC$est)+c(-1, 1)*qnorm(1-alpha/2)*sqrt(rubin.AUC$var)/(rubin.AUC$est*(1-rubin.AUC$est)))
          out$AUC <- c(est=rubin.AUC$est, LL=max(0, IC.auc[1]), UL=min(1, IC.auc[2]))
        }else{
          out$AUC <- c(est=rubin.AUC$est, LL=NA, UL=NA)
        }
      }
      ### th
      if ("all" %in% index | "threshold" %in% index){
        if (ci){
          IC.th <- rubin.th$est+c(-1, 1)*qnorm(1-alpha/2)*sqrt(rubin.th$var)
          out$threshold <- c(est=rubin.th$est, LL=IC.th[1], UL=IC.th[2])
        }else{
          out$threshold <- c(est=rubin.th$est, LL=NA, UL=NA)
        }
      }
      ### sens
      if ("all" %in% index | "sens" %in% index){
        if (ci){
          IC.sens <- rubin.sens$est+c(-1, 1)*qnorm(1-alpha/2)*sqrt(rubin.sens$var)
          out$sens <- c(est=rubin.sens$est, LL=max(0, IC.sens[1]), UL=min(1, IC.sens[2]))
        }else{
          out$sens <- c(est=rubin.sens$est, LL=NA, UL=NA)
        }
      }
      ### spec
      if ("all" %in% index | "spec" %in% index){
        if (ci){
          IC.spec <- rubin.spec$est+c(-1, 1)*qnorm(1-alpha/2)*sqrt(rubin.spec$var)
          out$spec <- c(est=rubin.spec$est, LL=max(0, IC.spec[1]), UL=min(1, IC.spec[2]))
        }else{
          out$spec <- c(est=rubin.spec$est, LL=NA, UL=NA)
        }
      }
    }
  }else if (method=="survivalROC"){
    ### AUC
    if ("all" %in% index | "AUC" %in% index){
      auc <- survivalROC(Stime=time, status=status, marker=cont.var, predict.time=predict.time, method="NNE", lambda=lambda)$AUC
      if (ci){
        df <- data.frame(cont.var, time, status)
        boot.aucs <- boot(df, statistic=survivalROCAUCboot, R=B, predict.time=predict.time, lambda=lambda)
        boot.aucs.ci <- boot.ci(boot.aucs, type=c("perc", "norm"), conf=1-alpha)
        # save
        IC <- rbind(boot.aucs.ci$normal[2:3], boot.aucs.ci$percent[4:5])
        colnames(IC) <- c("lower", "upper")
        rownames(IC) <- c("normal", "percentile")
        out$AUC <- c(est=auc, LL.norm=IC["normal", "lower"], UL.norm=IC["normal", "upper"], LL.perc=IC["percentile", "lower"], UL.perc=IC["percentile", "upper"])
      }else{
        out$AUC <- c(est=auc, LL=NA, UL=NA)
      }
    }
    if ("all" %in% index | "threshold" %in% index){
      ### threshold
      result.th <- minCostSurvivalROC(cont.var, time, status, predict.time, rho, costs2, lambda, plot, ci, alpha=alpha, B=B, ...)
      if (ci){
        out$threshold <- c(est=result.th$cut.est$optimum$cutoff, LL.norm=result.th$CI["normal", "lower"], UL.norm=result.th$CI["normal", "upper"], LL.perc=result.th$CI["percentile", "lower"], UL.perc=result.th$CI["percentile", "upper"])
      }else{
        out$threshold <- c(est=auc, LL=NA, UL=NA)
      }
    }
    if ("all" %in% index | "sens" %in% index){
      ### sens
      sens <- result.th$cut.est$sens
      if (ci){
        df <- data.frame(cont.var, time, status)
        boot.senss <- boot(df, statistic=bootminCostSurvivalROCsens, R=B, predict.time=predict.time, rho=rho, costs=costs2, lambda=lambda)
        boot.senss.ci <- boot.ci(boot.senss, type=c("perc", "norm"), conf=1-alpha)
        # save
        IC <- rbind(boot.senss.ci$normal[2:3], boot.senss.ci$percent[4:5])
        colnames(IC) <- c("lower", "upper")
        rownames(IC) <- c("normal", "percentile")
        out$sens <- c(est=sens, LL.norm=max(0, IC["normal", "lower"]), UL.norm=min(IC["normal", "upper"], 1), LL.perc=max(0, IC["percentile", "lower"]), UL.perc=min(1, IC["percentile", "upper"]))
      }else{
        out$sens <- c(est=sens, LL=NA, UL=NA)
      }
    }
    if ("all" %in% index | "spec" %in% index){
      ### spec
      spec <- result.th$cut.est$spec
      if (ci){
        df <- data.frame(cont.var, time, status)
        boot.specs <- boot(df, statistic=bootminCostSurvivalROCspec, R=B, predict.time=predict.time, rho=rho, costs=costs2, lambda=lambda)
        boot.specs.ci <- boot.ci(boot.specs, type=c("perc", "norm"), conf=1-alpha)
        # save
        IC <- rbind(boot.specs.ci$normal[2:3], boot.specs.ci$percent[4:5])
        colnames(IC) <- c("lower", "upper")
        rownames(IC) <- c("normal", "percentile")
        out$spec <- c(est=spec, LL.norm=max(0, IC["normal", "lower"]), UL.norm=min(IC["normal", "upper"], 1), LL.perc=max(0, IC["percentile", "lower"]), UL.perc=min(1, IC["percentile", "upper"]))
      }else{
        out$spec <- c(est=spec, LL=NA, UL=NA)
      }
    }
  }
  
  # output
  out$method <- method
  out$alpha <- alpha
  out$data <- data.frame(cont.var, time, status, statusNA=status.predict.time.NA)
  class(out) <- "diagnostic_assessment"
  return(out)
}


print.diagnostic_assessment <- function(x, ...){
  noms <- names(x)[-which(names(x) %in% c("method", "alpha", "data"))]
  cat("-------------------------\n")
  for (i in 1:length(noms)){
    cat(noms[i], ": ", x[[noms[i]]][1], "\n", sep="")
    if (!is.na(x[[noms[i]]][2])){
      if (length(x[[noms[i]]])==3){
        cat((1-x$alpha)*100, "% confidence interval: ", x[[noms[i]]][2], " - ", x[[noms[i]]][3], "\n", sep="")
      }else{
        cat((1-x$alpha)*100, "% confidence interval (bootstrap, normal): ", x[[noms[i]]][2], " - ", x[[noms[i]]][3], "\n", sep="")
        cat((1-x$alpha)*100, "% confidence interval (bootstrap, percentile): ", x[[noms[i]]][4], " - ", x[[noms[i]]][5], "\n", sep="")
      }
    }
    cat("-------------------------\n")
  }
}

