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
  out <- list(optimum=optimum, results=cut, sens=se[which.min(C)], spec=sp[which.min(C)])
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

bootminCostSurvivalROCsens <- function(df, indices, predict.time, rho, costs, lambda){
  # bootstrap resample
  df.boot <- df[indices, ]
  cut.boot <- with(df.boot, minCostSurvivalROCOpt(cont.var, time, status, predict.time, rho, costs, lambda, plot=FALSE))$sens
  return(as.numeric(cut.boot))
}

bootminCostSurvivalROCspec <- function(df, indices, predict.time, rho, costs, lambda){
  # bootstrap resample
  df.boot <- df[indices, ]
  cut.boot <- with(df.boot, minCostSurvivalROCOpt(cont.var, time, status, predict.time, rho, costs, lambda, plot=FALSE))$spec
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



#### from package survivalROC, by P Heagerty and P Saha-Chaudhuri
##
## survivalROC.R
##
## AUTHOR:  P. Heagerty
##
## DATE:  98/08/18
##
## (last modified:  98/11/17)
## Comments added by : P. SAHA
## DATE: 06/04/19
## Function name changed from roc.KM.calc() to survivalROC() as of 06/05/17
##---------------------------------------------------------------------------
survivalROC<-function( Stime, status, marker, entry=NULL,
                       predict.time, cut.values=NULL,
                       method="NNE", lambda=NULL, span=NULL,
                       window="symmetric" )
{
  ## DATE: July 12, 2006
  ## one of the option of method is changed from "smooth" to "NNE"
  ##
  ## PURPOSE:  calculations for Kaplan-Meier based ROC curve
  ##
  ##
  ###
  ### drop any missing
  ###
  ##
  times=Stime
  ## changed input from times to Stime, June 1, 2006, by Paramita Saha
  ## to make both the functions similar
  x <- marker
  
  if( is.null(entry) ) entry <- rep( 0, length(times) )
  ##
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if( sum(bad)>0 ) cat(paste("\n", sum(bad),
                             "records with missing values dropped. \n") )
  ## 
  if( is.null(cut.values) ) cut.values <- unique(x)
  cut.values <- cut.values[ order(cut.values) ]
  ncuts <- length(cut.values)
  ##
  ###
  ### sort the times
  ###
  ##
  ooo <- order( times )
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  ##
  ###
  ### overall survival probability
  ###
  ##
  s0 <- 1.0
  unique.t0 <- unique( times )
  unique.t0 <- unique.t0[ order(unique.t0) ]
  ## print( unique.t0 )
  n.times <- sum( unique.t0 <= predict.time )
  for( j in 1:n.times ){
    n<-sum( entry <= unique.t0[j] & times >= unique.t0[j] )
    d<-sum( (entry <= unique.t0[j])&(times==unique.t0[j])&(status==1) )
    if( n>0 ) s0<-s0*( 1 - d/n )	
  }
  s.pooled <- s0
  ## print( s.pooled )
  ##
  ###
  ### compute the sensitivity and specificity for each cut.value
  ###
  ##
  roc.matrix <- matrix( NA, ncuts, 2 )
  ##
  ## changed to strictly greater 98/11/17
  ##
  roc.matrix[ncuts,1] <- 0.0
  roc.matrix[ncuts,2] <- 1.0
  ##
  if( method=="KM" ){
    ##
    for( c in 1:(ncuts-1) ){
      ##
      s0 <- 1.0
      ##
      ## changed to strictly greater 98/11/17
      ##
      subset <- as.logical( x > cut.values[c] )
      ##
      e0 <- entry[subset]
      t0 <- times[subset]
      c0 <- status[subset]
      if( !is.null(t0) ){
        unique.t0 <- unique( t0 )
        unique.t0 <- unique.t0[ order(unique.t0) ]
        n.times <- sum( unique.t0 <= predict.time )
        if( n.times>0 ){
          for( j in 1:n.times ){
            n<-sum( e0 <= unique.t0[j] & t0 >= unique.t0[j] )
            d<-sum( (e0<=unique.t0[j])& (t0==unique.t0[j])&(c0==1) )
            if( n>0 ) s0<-s0*( 1 - d/n )	
          }
        }
      }
      p0 <- mean(subset) 
      roc.matrix[ c, 1 ] <-  (1-s0) * p0 / (1-s.pooled)
      roc.matrix[ c, 2 ] <-  1 - s0 * p0 / s.pooled
      ##
    }
  }## end-of method=="KM"
  ##
  ##
  if( method=="NNE" ){
    if( is.null(lambda) & is.null(span) ){
      ## cat("method = smooth requires either lambda or span! \n")
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[ order(x.unique) ]
    S.t.x <- rep( 0, length(x.unique) )
    t.evaluate <- unique( times[status==1] )
    t.evaluate <- t.evaluate[ order(t.evaluate) ]
    t.evaluate <- t.evaluate[ t.evaluate <= predict.time ]
    for( j in 1:length(x.unique) ){
      if( !is.null(span) ){
        if( window=="symmetric" ){
          ##
          ###
          ### symmetric
          ###
          ##
          ddd <- ( x - x.unique[j] )
          n <- length(x)
          ddd <- ddd[ order(ddd) ]
          index0 <- sum( ddd< 0 ) + 1
          index1 <- index0 + trunc(n*span +.5)
          
          ## Date: May 5, 2008
          ## Someone sent us an email saying that the same span
          ## produced very different results in .R vs .C
          ## code. We found out that .R was using interval
          ## length = n*span/2 on either side while .C was using
          ## length = n*span+.5 on either side, so the interval
          ## length considered in .C was double the length for
          ## .R. Today, we fixed the .R code from the following
          ## line to the line above.
          ## index1 <- index0 + trunc(n*span +.5)
          
          if( index1>n ) index1 <- n
          lambda <- ddd[ index1 ]
          wt <- as.integer( ( (x-x.unique[j]) <= lambda ) & 
                              ( (x-x.unique[j]) >= 0 ) )
          index0 <- sum( ddd<= 0 ) 
          index2 <- index0 - trunc(n*span/2)
          if( index2<1 ) index2 <- 1
          lambda <- abs( ddd[ index1 ] )
          set.index<- ( (x-x.unique[j]) >= -lambda ) & 
            ( (x-x.unique[j]) <= 0 ) 
          wt[ set.index ] <- 1
          
        }
        if( window=="asymmetric" ){
          ##
          ###
          ### asymmetric
          ###
          ##
          ddd <- ( x - x.unique[j] )
          n <- length(x)
          ddd <- ddd[ order(ddd) ]
          index0 <- sum( ddd< 0 ) + 1
          index <- index0 + trunc(n*span)
          if( index>n ) index <- n
          lambda <- ddd[ index ]
          wt <- as.integer( ( (x-x.unique[j]) <= lambda ) & 
                              ( (x-x.unique[j]) >= 0 ) )
        }
      }else{
        ## wt <- as.integer( abs(x-x.unique[j]) <= lambda )
        wt <- exp( -(x-x.unique[j])^2 / lambda^2 )
      }
      s0 <- 1.0
      for( k in 1:length(t.evaluate) ){
        n <- sum( wt*(entry<=t.evaluate[k])&(times >= t.evaluate[k]) )
        d <- sum( wt*(entry<=t.evaluate[k])&(times==t.evaluate[k])*
                    (status==1) )
        if( n > 0 ) s0 <- s0 * ( 1 - d/n )
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[ match( x, x.unique ) ]
    n <- length( times )
    S.marginal <- sum( S.all.x )/n
    for( c in 1:(ncuts-1) ){
      ##
      ## changed to strictly greater 98/11/17
      ##
      p1 <- sum( x > cut.values[c] )/n
      Sx <- sum( S.all.x[ x > cut.values[c] ] )/n
      roc.matrix[ c, 1 ] <- (p1-Sx)/(1-S.marginal)
      roc.matrix[ c, 2 ] <-  1 - Sx/S.marginal
    }
  }## end-of method=="NNE"
  ##
  sensitivity=roc.matrix[,1]
  specificity=roc.matrix[,2]
  x <- 1 - c( 0.0, specificity )
  y <- c( 1.0, sensitivity )
  ##
  n <- length( x )
  ##
  dx <- x[-n] - x[-1]
  mid.y <- ( y[-n] + y[-1] )/2
  ##
  ## print( sum(dx) )
  area <- sum( dx*mid.y )
  
  list(cut.values=c(-Inf, cut.values),
       ## renamed x as cut.values and changed from (NA,unique.x) to this one. June 1, 2006, By Paramita Saha
       ## sensitivity=roc.matrix[,1], Changed to TP and FP with 1.0 value added 06/05/15
       TP=y,
       ## specificity=roc.matrix[,2],
       FP=x,
       predict.time=predict.time, ## predict.time included in the return object, June 1, 2006
       Survival = s.pooled,
       AUC=area)
}

survivalROCAUCboot <- function(df, indices, predict.time, lambda){
  # bootstrap resample
  df.boot <- df[indices, ]
  auc.boot <- with(df.boot, survivalROC(Stime=time, status=status, marker=cont.var, predict.time=predict.time, method="NNE", lambda=lambda)$AUC)
  return(as.numeric(auc.boot))
}


### auxiliar functions for ThresholdROC (Skaltsa)

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

rubin_outliers <- function(matrix, range){
  # rubin
  var.b <- var(matrix[1, ], na.rm=TRUE)
  var.w <- mean(matrix[2, ]^2, na.rm=TRUE)
  est <- mean(matrix[1, ], na.rm=TRUE)
  # winsorize if there are outliers
  bp.se <- boxplot(matrix[2, ], range=range, plot=FALSE)
  bp.est <- boxplot(matrix[1, ], range=range, plot=FALSE)
  if (length(bp.se$out)>0){
    var.w <- winsor.mean(matrix[2, ]^2, trim=0.25)
  }
  if (length(bp.est$out)>0) {
    var.b <- winsor.var(matrix[1, ], trim=0.25)
    est <- winsor.mean(matrix[1, ], trim=0.25)
  }
  var <- var.w+var.b+var.b/ncol(matrix)
  return(list(est=est, var=var))
}

