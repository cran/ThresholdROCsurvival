\name{diagnostic_assessment_binary}
\alias{diagnostic_assessment_binary}
\title{
Diagnostic ability assessment for binary diagnostic tests
}
\description{
This function estimates sensitivity and specificity at a fixed time-point t for binary diagnostic tests with survival data by using two methods: 1) unknown status exclusion (USE), which excludes subjects with unknown status at t; and 2) imputation of censored times (ICT), a method based on multiple imputation. The status of the subjects at a certain time-point of interest t (the event occurred before or at t or not) is defined by the time-to-event variable.
}
\usage{
diagnostic_assessment_binary(binary.var, time, status, predict.time,
                             method=c("USE", "ICT"), index=c("all", "sens", "spec"),
                             m=10, ci=TRUE, alpha=0.05, range=3)
}
\arguments{
  \item{binary.var}{
binary variable to be used as predictor of the status. It should be a factor which two levels: - (negative, which indicates absence of the event) and + (positive, which indicates presence of the event)
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
method to be used in the estimation process. The user can choose between \code{USE} (unknown status exclusion) or \code{ICT} (imputation of censored times). Default, \code{USE}  
}
  \item{index}{
indices to be estimated. The user can choose one or more of the following: \code{sens} and \code{spec}. The option \code{all} (default) estimates all two indices
}
  \item{m}{
the number of data sets to impute. Default, 10
}
  \item{ci}{
Should a confidence interval be calculated? Default, \code{TRUE}
}
  \item{alpha}{
significance level for the confidence interval. Default, 0.05
}
  \item{range}{
this value, which is passed to \code{boxplot} function from \code{graphics} package, determines the data points that are considered to be extreme in the estimates and standard errors from the multiple imputation process. We consider extreme observations those that exceed \code{range} times the interquartile range. If extreme observations are found in the estimates or standard errors from the multiple imputation process, Winsorized estimators (Wilcox, 2012) are used for the point estimate and the between and within variances. Default, 3
}
}
\details{
When \code{method} is \code{USE}: First, the algorithm determines the status of the subjects at time \code{predict.time}. Those censored subjects whose status could be not be determined are excluded from the analysis. Then, diagnostic ability is assessed with standard methods in the binary setting.

When \code{method} is \code{ICT}: First, the algorithm determines the status of the subjects at time \code{predict.time}. For those subjects whose status could not be determined because their censored time is lower than t, we impute survival times using the method of Hsu et al (2006), that is implemented in the package \code{InformativeCensoring} (Ruau et al, 2020). The status of the subjects is then determined by these imputed times and is used to estimate the indices in \code{index}. Confidence intervals are calculated using the standard error proposed by Rubin (1987).
}
\value{
An object of class \code{diagnostic_assessment}, which is a list with the following components:
\item{sens }{Sensitivity estimate and its corresponding confidence interval (if \code{ci=TRUE}), only if sensitivity has been included in \code{index}}
\item{spec }{Specificity estimate and its corresponding confidence interval (if \code{ci=TRUE}), only if specificity has been included in \code{index}}
\item{method }{\code{method} used in the estimation}
\item{alpha }{significance level provided by the user}
\item{data }{A data.frame containing the following columns previously provided by the user: \code{cont.var}, \code{time} and \code{status}, and a new column \code{statusNA}, which contains the status of the subjects at time \code{predict.time} (0=no event, 1=event, \code{NA}=unknown)}
}
\references{
Heagerty PJ, Lumley T, Pepe MS. Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic Marker. Biometrics 2000; 56(2): 337-344. doi: 10.1111/j.0006-341X.2000.00337.x

Heagerty PJ, Saha-Chaudhuri P (2013). survivalROC: Time-dependent ROC curve estimation from censored survival data. R package version 1.0.3. \url{https://CRAN.R-project.org/package=survivalROC}

Hsu CH, Taylor JMG, Murray S, Commenges D. Survival analysis using auxiliary variables via non-parametric multiple imputation. Statistics in Medicine 2006; 25(20): 3503-3517. doi: https://doi.org/10.1002/sim.2452

Ruau D, Burkoff N, Bartlett J, Jackson D, Jones E, Law M and Metcalfe P (2020). InformativeCensoring: Multiple Imputation for Informative Censoring. R package version 0.3.5. \url{https://CRAN.R-project.org/package=InformativeCensoring}

Rubin DB. Multiple Imputation for Nonresponse in Surveys. Wiley Series in Probability and Statistics. John Wiley & Sons (1987)

Wilcox, R. Introduction to Robust Estimation and Hypothesis Testing. 3rd Edition. Elsevier, Amsterdam (2012)

Zhou XH, Obuchowski NA and McClish DK. Statistical methods in diagnostic medicine. John Wiley and sons (2002)
}

\seealso{
\code{\link{diagnostic_assessment_continuous}}
}
\examples{
data(NSCLC)
NSCLC$COL_cat <- factor(ifelse(NSCLC$COL>=10, "+", "-"))
set.seed(2020)
with(NSCLC, diagnostic_assessment_binary(COL_cat, OS, ST,
     1095, method="ICT", m=10, ci=TRUE))
}
