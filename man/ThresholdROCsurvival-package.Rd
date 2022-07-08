\name{ThresholdROCsurvival-package}
\alias{ThresholdROCsurvival-package}
\docType{package}
\title{
Optimum threshold and AUC estimation based on cost function in a context of right-censored data
}
\description{
We focus on the estimation of optimal thresholds and AUCs when the outcome of interest is the status (alive or dead) of the subjects at a certain time-point t. This binary status is determined by right-censored times to event and it is missing for those subjects censored before t. Here we provide three methods (missing exclusion, imputation of censored times and using time-dependent ROC curves) to estimate optimal thresholds and AUCs in this context. Two references for the methods used here are Skaltsa et al. (2010) <doi:10.1002/bimj.200900294> and Heagerty et al. (2000) <doi:10.1111/j.0006-341x.2000.00337.x>.
}
\details{
The functions in this package are \code{th_ME()}, \code{th_ICT()} and \code{th_survivalROC()} for threshold estimation and \code{AUC_ME()} and \code{AUC_ICT()} for AUC estimation.
}
\author{
\packageAuthor{ThresholdROCsurvival}

Maintainer: \packageMaintainer{ThresholdROCsurvival}
}
\references{
Heagerty PJ, Lumley T, Pepe MS. Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic Marker. Biometrics 2000; 56(2): 337-344. doi: 10.1111/j.0006-341X.2000.00337.x

Hsu CH, Taylor JMG, Murray S, Commenges D. Survival analysis using auxiliary variables via non-parametric multiple imputation. Statistics in Medicine 2006; 25(20): 3503-3517. doi: https://doi.org/10.1002/sim.2452

Perez-Jaume S, Skaltsa K, Pallares N, Carrasco JL. ThresholdROC: Optimum Threshold Estimation Tools for Continuous Diagnostic Tests in R. Journal of Statistical Software 2017; 82(4): 1-21. doi: 10.18637/jss.v082.i04

Skaltsa K, Jover L, Carrasco JL. Estimation of the diagnostic threshold accounting for decision costs and sampling uncertainty. Biometrical Journal 2010; 52(5): 676-697. doi: 10.1002/bimj.200900294
}
\keyword{ package }