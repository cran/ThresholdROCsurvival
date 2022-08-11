\name{th_ME}
\alias{th_ME}
\title{
Threshold estimation using the missing exclusion (ME) method}
\description{
This function estimates optimum thresholds with survival data by excluding subjects with missing status at the time-point of interest
}
\usage{
th_ME(cont.var, time, status, predict.time, costs = NULL,
      R = NULL, method = c("normal", "empirical"),
      var.equal = FALSE, plot = FALSE, ci = TRUE,
      alpha = 0.05, B = 1000, ...)
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
When \code{method="normal"}, assume equal variances? Default, \code{FALSE}. When \code{method="empirical"}, \code{var.equal} is ignored
}
  \item{plot}{
Should some graphs about the estimation be plotted? Default, \code{FALSE}
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
  \item{\dots}{
further arguments to be passed to \code{plot()}
}
}
\details{
First, the algorithm determines the status of the subjects at time \code{predict.time}. Those censored subjects whose status could be not be determined (therefore, with missing status) are excluded from the analysis. Then, the optimum threshold is estimated by minimizing the cost function using the \code{thres2} function in the \code{ThresholdROC} package (Perez-Jaume et al, 2017).

For parameter \code{method} the user can choose between \code{"normal"} (assumes binormality) or \code{"empirical"} (leaves out any distributional assumption). When \code{method="normal"}, the user can specify if the algorithm should assume equal or different variances using the parameter \code{var.equal}. For further details see the \code{thres2} function in the \code{ThresholdROC} package.

R, mentioned in parameters \code{costs} and \code{R}, is the product of the non-disease odds and the cost ratio:
\deqn{R=((1-p)/p)((C_{TN}-C_{FP})/(C_{TP}-C_{FN})),}
where \emph{p} is the disease prevalence (estimated using Kaplan-Meier) and \emph{C_i} are the classification costs.

To calculate sensitivity, specificity and predictive values corresponding to the estimated threshold, we suggest to use the \code{diagnostic} function in the \code{ThresholdROC} package.
}
\value{
An object of class \code{thres2}, which is a list of two components (see the help on the \code{thres2} function). Here we add a third component, \code{data}: a data.frame containing the following columns previously provided by the user: \code{cont.var}, \code{time} and \code{status}, and a new column \code{statusNA}, which contains the status of the subjects at time \code{predict.time} (0=no event, 1=event, \code{NA}=missing)
}
\references{
Perez-Jaume S, Skaltsa K, Pallares N, Carrasco JL. ThresholdROC: Optimum Threshold Estimation Tools for Continuous Diagnostic Tests in R. Journal of Statistical Software 2017; 82(4): 1-21. doi: 10.18637/jss.v082.i04

Skaltsa K, Jover L, Carrasco JL. Estimation of the diagnostic threshold accounting for decision costs and sampling uncertainty. Biometrical Journal 2010; 52(5): 676-697. doi: 10.1002/bimj.200900294
}

\seealso{
\code{\link{thres2}}
}
\examples{
data(NSCLC)
res <- with(NSCLC, th_ME(log(COL), OS, ST, 1095, method="normal",
            var.equal=FALSE, plot=TRUE, xlab="Collagen"))
res
exp(res$T$thres)
exp(res$CI$lower)
exp(res$CI$upper)
res$data
}
