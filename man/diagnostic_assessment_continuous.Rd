\name{diagnostic_assessment_continuous}
\alias{diagnostic_assessment_continuous}
\title{
Diagnostic ability assessment for continuous diagnostic tests
}
\description{
This function estimates the AUC, optimal threshold, sensitivity and specificity at a fixed time-point t for continuous diagnostic tests with survival data by using three methods: 1) missing exclusion (ME), which excludes subjects with missing status at t; 2) imputation of censored times (ICT), a method based on multiple imputation; and 3) survivalROC, which uses a method based on time-dependent ROC curves. The status of the subjects at a certain time-point of interest t (the event occurred before or at t or not) is defined by the time-to-event variable.
}
\usage{
diagnostic_assessment_continuous(cont.var, time, status, predict.time,
                                 method=c("ME", "ICT", "survivalROC"),
                                 index=c("all", "AUC", "threshold", "sens", "spec"),
                                 costs=NULL, R=NULL,
                                 method.thres=c("normal", "empirical"),
                                 var.equal=FALSE, lambda=0.05, m=10,
                                 ci=TRUE, plot=FALSE, alpha=0.05,
                                 B=1000, range=3, ...)
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
  \item{method}{
method to be used in the estimation process. The user can choose between \code{ME} (missing exclusion), \code{ICT} (imputation of censored times) or \code{survivalROC} (time-dependent ROC curves). Default, \code{ME}  
}
  \item{index}{
indices to be estimated. The user can choose one or more of the following: \code{AUC}, \code{threshold}, \code{sens} (sensitivity achieved by the optimal threshold), \code{spec} (specificity achieved by the optimal threshold). The option \code{all} (default) estimates all four indices
}
  \item{costs}{
cost matrix. Costs should be entered as a 2x2 matrix, where the first row corresponds to the true positive and true negative costs and the second row to the false positive and false negative costs. Default cost values (\code{costs=NULL}, when also \code{R=NULL}) are a combination of costs that yields R=1, which is equivalent to the Youden index method (for details about this concept, see \bold{Details} and \bold{References})
}
  \item{R}{
if the cost matrix \code{costs} is not set (that is, \code{costs=NULL}), R desired (the algorithm will choose a suitable combination of costs that leads to \code{R}). Default, \code{NULL} (which leads to R=1 using the default \code{costs}). For details about this concept, see \bold{Details} and \bold{References}
}
  \item{method.thres}{
method used in the estimation: \code{"normal"} (default) or \code{"empirical"}. The user can specify just the initial letters. See \bold{Details} for more information about the methods available
}
  \item{var.equal}{
when \code{method="normal"}, assume equal variances? Default, \code{FALSE}. When \code{method="empirical"}, \code{var.equal} is ignored
}
  \item{lambda}{
smoothing parameter for the NNE algorithm used in \code{survivalROC()} function. Default, 0.05
}
  \item{m}{
the number of data sets to impute. Default, 10
}
  \item{ci}{
Should a confidence interval be calculated? Default, \code{TRUE}
}
  \item{plot}{
Should some graphs about the estimation be plotted? Default, \code{FALSE}
}
  \item{alpha}{
significance level for the confidence interval. Default, 0.05
}
  \item{B}{
number of bootstrap resamples for the confidence interval. Only used when \code{method} is \code{survivalROC}. Otherwise, this argument is ignored. Default, 1000
}
  \item{range}{
this value, which is passed to \code{boxplot} function from \code{graphics} package, determines the data points that are considered to be extreme in the estimates and standard errors from the multiple imputation process. We consider extreme observations those that exceed \code{range} times the interquartile range. If extreme observations are found in the estimates or standard errors from the multiple imputation process, Winsorized estimators (Wilcox, 2012) are used for the point estimate and the between and within variances. Default, 3
}
  \item{...}{
extra arguments to be passed to \code{plot()}  
}
}
\details{
When \code{method} is \code{ME}: First, the algorithm determines the status of the subjects at time \code{predict.time}. Those censored subjects whose status could be not be determined (therefore, with missing status) are excluded from the analysis. Then, diagnostic ability is assessed with standard methods in the binary setting.

When \code{method} is \code{ICT}: First, the algorithm determines the status of the subjects at time \code{predict.time}. For those subjects whose status could not be determined because their censored time is lower than t (therefore, with missing status), we impute survival times using the method of Hsu et al (2006), that is implemented in the package \code{InformativeCensoring} (Ruau et al, 2020). The status of the subjects is then determined by these imputed times and is used to estimate the indices in \code{index}. Confidence intervals are calculated using the standard error proposed by Rubin (1987).

When \code{method} is \code{survivalROC}: Diagnostic ability is assessed by constructing the ROC curve at time t through time-dependent ROC curves (Heagerty et al, 2000). Confidence intervals are obtained using normal and percentile bootstrap. In normal bootstrap, the bootstrap is used to obtain an estimate of the standard error of the threshold estimate, and then the standard normal distribution is used for the confidence interval calculation. In percentile bootstrap, \code{B} bootstrap resamples are generated and the threshold is estimated in all of them. Then, the confidence interval is calculated as the empirical 1-\code{alpha}/2 and 1+\code{alpha}/2 percentiles of the \code{B} bootstrap estimates.

For parameter \code{method.thres}, the method used in the estimation of the optimal threshold, the user can choose between \code{"normal"} (assumes binormality) or \code{"empirical"} (leaves out any distributional assumption). When \code{method="normal"}, the user can specify if the algorithm should assume equal or different variances using the parameter \code{var.equal}. For further details see the \code{thres2} function in the \code{ThresholdROC} package.

R, mentioned in parameters \code{costs} and \code{R}, is the product of the non-disease odds and the cost ratio:
\deqn{R=((1-p)/p)((C_{TN}-C_{FP})/(C_{TP}-C_{FN})),}
where \emph{p} is the disease prevalence (estimated using Kaplan-Meier) and \emph{C_i} are the classification costs.
}
\value{
An object of class \code{diagnostic_assessment}, which is a list with the following components:
\item{AUC }{AUC estimate and its corresponding confidence interval (if \code{ci=TRUE}), only if AUC has been included in \code{index}}
\item{threshold }{threshold estimate and its corresponding confidence interval (if \code{ci=TRUE}), only if threshold has been included in \code{index}}
\item{sens }{Sensitivity estimate (achieved by the optimal threshold) and its corresponding confidence interval (if \code{ci=TRUE}), only if sensitivity has been included in \code{index}}
\item{spec }{Specificity estimate (achieved by the optimal threshold) and its corresponding confidence interval (if \code{ci=TRUE}), only if specificity has been included in \code{index}}
\item{method }{\code{method} used in the estimation}
\item{alpha }{significance level provided by the user}
\item{data }{A data.frame containing the following columns previously provided by the user: \code{cont.var}, \code{time} and \code{status}, and a new column \code{statusNA}, which contains the status of the subjects at time \code{predict.time} (0=no event, 1=event, \code{NA}=missing)}
}
\references{
Heagerty PJ, Lumley T, Pepe MS. Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic Marker. Biometrics 2000; 56(2): 337-344. doi: 10.1111/j.0006-341X.2000.00337.x

Heagerty PJ, Saha-Chaudhuri P (2022). survivalROC: Time-dependent ROC curve estimation from censored survival data. R package version 1.0.3.1. \url{https://CRAN.R-project.org/package=survivalROC}

Hsu CH, Taylor JMG, Murray S, Commenges D. Survival analysis using auxiliary variables via non-parametric multiple imputation. Statistics in Medicine 2006; 25(20): 3503-3517. doi: https://doi.org/10.1002/sim.2452

Kottas M, Kuss O, Zapf A. A modified Wald interval for the area under the ROC curve (AUC) in diagnostic case-control studies. BMC Medical Research Methodology 2014; 14(26). doi:10.1186/1471-2288-14-26

Perez-Jaume S, Skaltsa K, Pallares N, Carrasco JL. ThresholdROC: Optimum Threshold Estimation Tools for Continuous Diagnostic Tests in R. Journal of Statistical Software 2017; 82(4): 1-21. doi: 10.18637/jss.v082.i04

Ruau D, Burkoff N, Bartlett J, Jackson D, Jones E, Law M and Metcalfe P (2020). InformativeCensoring: Multiple Imputation for Informative Censoring. R package version 0.3.5. \url{https://CRAN.R-project.org/package=InformativeCensoring}

Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, Muller M. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 2011; 12. doi:10.1186/1471-2105-12-77

Rubin DB. Multiple Imputation for Nonresponse in Surveys. Wiley Series in Probability and Statistics. John Wiley & Sons (1987)

Skaltsa K, Jover L, Carrasco JL. Estimation of the diagnostic threshold accounting for decision costs and sampling uncertainty. Biometrical Journal 2010; 52(5): 676-697. doi: 10.1002/bimj.200900294

Wilcox, R. Introduction to Robust Estimation and Hypothesis Testing. 3rd Edition. Elsevier, Amsterdam (2012)
}

\seealso{
\code{\link{diagnostic_assessment_binary}}
}
\examples{
library(ThresholdROCsurvival)
data(NSCLC)
# missing exclusion (Youden index maximization, R=1)
with(NSCLC, diagnostic_assessment_continuous(log(COL), OS,
                                             ST, 1095, method="ME", method.thres="normal",
                                             var.equal=FALSE, ci=TRUE))
\donttest{
# multiple imputation (Youden index maximization, R=1)
set.seed(2020)
with(NSCLC, diagnostic_assessment_continuous(log(COL), OS,
                                             ST, 1095, method="ICT", method.thres="normal",
                                             var.equal=FALSE, m=50, ci=TRUE))
}

# missing exclusion (R=1.1)
with(NSCLC, diagnostic_assessment_continuous(log(COL), OS,
                                             ST, 1095, method="ME", method.thres="normal",
                                             var.equal=FALSE, ci=TRUE, R=1.1))
\donttest{
# multiple imputation (R=1.1)
set.seed(2020)
with(NSCLC, diagnostic_assessment_continuous(log(COL), OS,
                                             ST, 1095, method="ICT", method.thres="normal",
                                             var.equal=FALSE, m=50, ci=TRUE, R=1.1))
}
}
