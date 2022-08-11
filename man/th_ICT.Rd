\name{th_ICT}
\alias{th_ICT}
\title{
Threshold estimation using the imputation of censored times (ICT) method
}
\description{
This function estimates optimum thresholds with survival data by using a method based on the imputation of censored times (ICT). The status of the subjects at a certain time-point of interest t (the event ocurred before or at t or not) is defined by the time-to-event variable.
}
\usage{
th_ICT(cont.var, time, status, predict.time, costs = NULL,
       R = NULL, method = c("normal", "empirical"),
       var.equal = FALSE, m = 10, ci = TRUE, alpha = 0.05,
       B = 1000, range = 3)
}
\arguments{
  \item{cont.var}{
continuous variable or biomarker to be used as predictor of the status
}
  \item{time}{
survival time
}
  \item{status}{
censoring status codified as 0=censored, 1=event
}
  \item{predict.time}{
time-point of interest to define the subjects' status as event present or absent
}
  \item{costs}{
cost matrix. Costs should be entered as a 2x2 matrix, where the first row corresponds to the true positive and true negative costs and the second row to the false positive and false negative costs. Default cost values (\code{costs=NULL}, when also \code{R=NULL}) are a combination of costs that yields R=1, which is equivalent to the Youden index method (for details about this concept, see \bold{Details} and \bold{References})
}
  \item{R}{
if the cost matrix \code{costs} is not set (that is, \code{costs=NULL}), R desired (the algorithm will choose a suitable combination of costs that leads to \code{R}). Default, \code{NULL} (which leads to R=1 using the default \code{costs}). For details about this concept, see \bold{Details} and \bold{References}
}
  \item{method}{
method used in the estimation: \code{"normal"} (default) or \code{"empirical"}. The user can specify just the initial letters. See \bold{Details} for more information about the methods available
}
  \item{var.equal}{
when \code{method="normal"}, assume equal variances? Default, \code{FALSE}. When \code{method="empirical"}, \code{var.equal} is ignored
}
  \item{m}{
the number of data sets to impute
}
  \item{ci}{
Should a confidence interval be calculated? Default, \code{TRUE}
}
  \item{alpha}{
significance level for the confidence interval. Default, 0.05
}
  \item{B}{
number of bootstrap resamples for the confidence interval when \code{method="empirical"}. Otherwise, ignored. Default, 1000
}
  \item{range}{
this value, which is passed to \code{boxplot} function from \code{graphics} package, determines the data points that are considered to be extreme in the estimates and standard errors from the multiple imputation process. We consider extreme observations those that exceed \code{range} times the interquartile range. If extreme observations are found in the estimates or standard errors from the multiple imputation process, Winsorized estimators (Wilcox, 2012) are used for the point estimate of the threshold and the between and within variances. Default, 3
}
}
\details{
First, the algorithm determines the status of the subjects at time \code{predict.time}. For those subjects whose status could not be determined because their censored time is lower than t (therefore, with missing status), we impute survival times using the method of Hsu et al (2006), that is implemented in the package \code{InformativeCensoring} (Ruau et al, 2020). The status of the subjects is then determined by these imputed times and is used to estimate the optimum threshold by minimizing the cost function using the \code{thres2} function in the \code{ThresholdROC} package (Perez-Jaume et al, 2017).

For parameter \code{method} the user can choose between \code{"normal"} (assumes binormality) or \code{"empirical"}. When \code{method="normal"}, the user can specify if the algorithm should assume equal or different variances using the parameter \code{var.equal}. For further details see the \code{thres2} function in the \code{ThresholdROC} package.

Confidence intervals are calculated using the standard error proposed by Rubin (1987).

R, mentioned in parameters \code{costs} and \code{R}, is the product of the non-disease odds and the cost ratio:
\deqn{R=((1-p)/p)((C_{TN}-C_{FP})/(C_{TP}-C_{FN})),}
where \emph{p} is the disease prevalence (estimated using Kaplan-Meier) and \emph{C_i} are the classification costs.
}
\value{
An object of class \code{th_ICT}, which is a list with the following components:

\item{T }{A list of five elements:

  \code{thres}  threshold estimate.
  
  \code{prev}  disease prevalence used.
  
  \code{costs}  cost matrix.
  
  \code{R}  R term, the product of the non-disease odds and the cost ratio (for further details about this concept, see \bold{References}).
  
  \code{method}  method used in the estimation.
  
  }
  
\item{CI }{A list of five elements (or \code{NULL} if \code{ci=FALSE}):

  \code{lower}  the lower limit of the confidence interval.

  \code{upper}  the upper limit of the confidence interval.
  
  \code{se}  the standard error. 
  
  \code{alpha}  significance level provided by the user.
  
  \code{ci.method} method used for the confidence intervals calculation.
  
  }

\item{sens}{A list of three elements:

  \code{est} estimate of sensitivity at \code{thres}.
  
  \code{lower} the lower limit of the confidence interval for the estimate of sensitivity at \code{thres} (\code{NULL} if \code{ci=FALSE}).
  
  \code{upper} the upper limit of the confidence interval for the estimate of sensitivity at \code{thres} (\code{NULL} if \code{ci=FALSE}).
  
  }
  
\item{spec}{A list of three elements:

  \code{est} estimate of specificity at \code{thres}.
  
  \code{lower} the lower limit of the confidence interval for the estimate of specificity at \code{thres} (\code{NULL} if \code{ci=FALSE}).
  
  \code{upper} the upper limit of the confidence interval for the estimate of specificity at \code{thres} (\code{NULL} if \code{ci=FALSE}).
  
  }

\item{data }{A data.frame containing the following columns previously provided by the user: \code{cont.var}, \code{time} and \code{status}, and a new column \code{statusNA}, which contains the status of the subjects at time \code{predict.time} (0=no event, 1=event, \code{NA}=missing)}
}
\references{
Hsu CH, Taylor JMG, Murray S, Commenges D. Survival analysis using auxiliary variables via non-parametric multiple imputation. Statistics in Medicine 2006; 25(20): 3503-3517. doi: https://doi.org/10.1002/sim.2452

David Ruau, Nikolas Burkoff, Jonathan Bartlett, Dan Jackson, Edmund Jones, Martin Law and Paul Metcalfe (2020). InformativeCensoring: Multiple Imputation for Informative Censoring. R package version 0.3.5. https://CRAN.R-project.org/package=InformativeCensoring

Perez-Jaume S, Skaltsa K, Pallares N, Carrasco JL. ThresholdROC: Optimum Threshold Estimation Tools for Continuous Diagnostic Tests in R. Journal of Statistical Software 2017; 82(4): 1-21. doi: 10.18637/jss.v082.i04

Rubin DB. Multiple Imputation for Nonresponse in Surveys. Wiley Series in Probability and Statistics. John Wiley & Sons (1987).

Skaltsa K, Jover L, Carrasco JL. Estimation of the diagnostic threshold accounting for decision costs and sampling uncertainty. Biometrical Journal 2010; 52(5): 676-697. doi: 10.1002/bimj.200900294

Wilcox, R. Introduction to Robust Estimation and Hypothesis Testing. 3rd Edition. Elsevier, Amsterdam (2012)
}
\seealso{
\code{\link{thres2}}
}
\examples{
data(NSCLC)
set.seed(2020)
res <- with(NSCLC, th_ICT(log(COL), OS, ST, 1095, method="normal", var.equal=FALSE, m=50))
res
exp(res$T$thres)
exp(res$CI$lower)
exp(res$CI$upper)
res$data
res$sens
res$spec
}
