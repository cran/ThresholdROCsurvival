### auxiliar functions for survivalROC

minCostSurvivalROCOpt <- function(cont.var, time, status, predict.time, rho, costs, lambda, plot, ...){
  cut <- survivalROC(Stime=time, status=status, marker=cont.var, predict.time=predict.time, method="NNE", lambda=lambda)
  # optimum
  cut2 <- data.frame(cut=cut$cut.values[is.finite(cut$cut.values)], TP=cut$TP[is.finite(cut$cut.values)], FP=cut$FP[is.finite(cut$cut.values)])
  #J <- with(cut2, TP+(1-FP)-1) # particular case - Youden index
  se <- with(cut2, TP)
  sp <- with(cut2, 1-FP)
  # cost function
  n <- length(cont.var)
  c.t.pos <- costs[1, 1]
  c.t.neg <- costs[1, 2]
  c.f.pos <- costs[2, 1]
  c.f.neg <- costs[2, 2]
  D <- n*rho*(c.t.pos-c.f.neg)
  R <- (1-rho)/rho*(c.t.neg-c.f.pos)/(c.t.pos-c.f.neg)
  G <- (rho*c.f.neg+(1-rho)*c.f.pos)/(rho*(c.t.pos-c.f.neg))
  C <- D*(se+sp*R+G)
  # look for optimum
  optimum <- data.frame(cutoff=cut2$cut[which.min(C)], C=C[which.min(C)])
  out <- list(optimum=optimum, results=cut)
  # plot cont.var vs youden index
  if (plot){
    df.plot <- data.frame(cont.var=cut2$cut, C)
    colnames(df.plot) <- c("cont.var", "Cost.function")
    with(df.plot, plot(cont.var, Cost.function, type="l", ...))
    with(df.plot, points(cont.var, Cost.function, pch=16))
    points(optimum$cutoff, optimum$C, col=2, pch=16)
  }
  return(out)
}

bootminCostSurvivalROC <- function(df, indices, predict.time, rho, costs, lambda){
  # bootstrap resample
  df.boot <- df[indices, ]
  cut.boot <- with(df.boot, minCostSurvivalROCOpt(cont.var, time, status, predict.time, rho, costs, lambda, plot=FALSE))$optimum$cutoff
  return(as.numeric(cut.boot))
}

minCostSurvivalROC <- function(cont.var, time, status, predict.time, rho, costs, lambda, plot, ci, alpha, B, ...){
  cut.est <- minCostSurvivalROCOpt(cont.var, time, status, predict.time, rho, costs, lambda, plot, ...)
  if (ci){
    df <- data.frame(cont.var, time, status)
    boot.cuts <- boot(df, statistic=bootminCostSurvivalROC, R=B, predict.time=predict.time, rho=rho, costs=costs, lambda=lambda)
    boot.cuts.ci <- boot.ci(boot.cuts, type=c("perc", "norm"), conf=1-alpha)
    # save
    IC <- rbind(boot.cuts.ci$normal[2:3], boot.cuts.ci$percent[4:5])
    colnames(IC) <- c("lower", "upper")
    rownames(IC) <- c("normal", "percentile")
    se <- sd(boot.cuts$t)
  }else{
    IC <- NULL
    se <- NULL
  }
  out <- list(cut.est=cut.est, CI=IC, se=se)
  # return
  return(out)
}




### auxiliar functions for ThesholdROC (Skaltsa)

minCostThresholdROC <- function(cont.var, status, rho, costs, method, var.equal, plot, ci, alpha, B, ...){
  ngroups <- dim(table(status))
  if (ngroups<2){
    cut.est <- NA
  }else{
    df <- data.frame(cont.var, status)
    cont.var.dead <- subset(df, status==1)$cont.var
    cont.var.alive <- subset(df, status==0)$cont.var
    # cutoff estimation
    if (method=="normal"){
      if (var.equal){
        cut <- thres2(cont.var.alive, cont.var.dead, rho, costs, method="equal", ci=ci, alpha=alpha, na.rm=TRUE)
      }else{
        cut <- suppressWarnings(thres2(cont.var.alive, cont.var.dead, rho, costs, method="unequal", ci=ci, alpha=alpha, na.rm=TRUE))
      }
    }else if (method=="empirical"){
      cut <- thres2(cont.var.alive, cont.var.dead, rho, costs, method="empirical", ci=ci, ci.method="boot", B=B, alpha=alpha, na.rm=TRUE, extra.info=TRUE)
    }else if (method=="smooth"){
      cut <- thres2(cont.var.alive, cont.var.dead, rho, costs, method="smooth", ci=ci, ci.method="boot", B=B, alpha=alpha, na.rm=TRUE, extra.info=TRUE)
    }
  }
  if (plot){
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    plot(cut, ...)
    par(mfrow=c(1, 2))
    plotCostROC(cut)
  }
  return(cut)
}




### other auxiliar functions

expit <- function(x){
  exp(x)/(exp(x)+1)
}

logit <- function(p){
  log(p/(1-p))
}
