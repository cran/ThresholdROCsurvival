\name{AUC_ICT}
\alias{AUC_ICT}
\title{
AUC estimation using the imputation of censored times (ICT) method
}
\description{
This function estimates the AUC with survival data by using a method based on the imputation of censored times (ICT). The status of the subjects at a certain time-point of interest t (the event ocurred before or at t or not) is defined by the time-to-event variable.
}
\usage{
AUC_ICT(cont.var, time, status, predict.time, m = 10,
        ci = TRUE, alpha = 0.05, range = 3)
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
  \item{m}{
the number of data sets to impute
}
  \item{ci}{
Should a confidence interval be calculated? Default, \code{TRUE}
}
  \item{alpha}{
significance level for the confidence interval. Default, 0.05
}
  \item{range}{
this value, which is passed to \code{boxplot} function from \code{graphics} package, determines the data points that are considered to be extreme in the estimates and standard errors from the multiple imputation process. We consider extreme observations those that exceed \code{range} times the interquartile range. If extreme observations are found in the estimates or standard errors from the multiple imputation process, Winsorized estimators (Wilcox, 2012) are used for the point estimate of the threshold and the between and within variances. Default, 3
}
}
\details{
First, the algorithm determines the status of the subjects at time \code{predict.time}. For those subjects whose status could not be determined because their censored time is lower than t (therefore, with missing status), we impute survival times using the method of Hsu et al (2006), that is implemented in the package \code{InformativeCensoring} (Ruau et al, 2020). The status of the subjects is then determined by these imputed times and is used to estimate the AUC using the \code{roc} function from \code{pROC} package (Robin et al, 2011).

Confidence intervals are calculated using the standard error proposed by Rubin (1987).
}
\value{
An object of class \code{AUC_ICT}, which is a list with the following components:
\item{estimate }{AUC estimate}
\item{se }{Standard error of the estimate, obtained using Rubin rules (note: \code{NULL} if \code{ci=FALSE})}
\item{CI }{Confidence interval for the estimate, obtained using Rubin rules (note: \code{NULL} if \code{ci=FALSE})}
\item{data }{A data.frame containing the following columns previously provided by the user: \code{cont.var}, \code{time} and \code{status}, and a new column \code{statusNA}, which contains the status of the subjects at time \code{predict.time} (0=no event, 1=event, \code{NA}=missing)}
}
\references{
Hsu CH, Taylor JMG, Murray S, Commenges D. Survival analysis using auxiliary variables via non-parametric multiple imputation. Statistics in Medicine 2006; 25(20): 3503-3517. doi: https://doi.org/10.1002/sim.2452

David Ruau, Nikolas Burkoff, Jonathan Bartlett, Dan Jackson, Edmund Jones, Martin Law and Paul Metcalfe (2020). InformativeCensoring: Multiple Imputation for Informative Censoring. R package version 0.3.5. https://CRAN.R-project.org/package=InformativeCensoring

Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, Muller M. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 2011; 12. doi:10.1186/1471-2105-12-77

Rubin DB. Multiple Imputation for Nonresponse in Surveys. Wiley Series in Probability and Statistics. John Wiley & Sons (1987).

Wilcox, R. Introduction to Robust Estimation and Hypothesis Testing. 3rd Edition. Elsevier, Amsterdam (2012)
}

\seealso{
\code{\link{AUC_ME}}
}
\examples{
data(NSCLC)
set.seed(2020)
res <- with(NSCLC, AUC_ICT(COL, OS, ST, 1095, m=50))
res
res$data
}
