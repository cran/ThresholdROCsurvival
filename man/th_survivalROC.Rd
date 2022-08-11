\name{th_survivalROC}
\alias{th_survivalROC}
\title{
Threshold estimation using the method based on time-dependent ROC curves (survivalROC)}
\description{
This function estimates optimum thresholds with survival data by using a method based on time-dependent ROC curves
}
\usage{
th_survivalROC(cont.var, time, status, predict.time, costs = NULL,
               R = NULL, lambda = 0.05, plot = FALSE, ci = FALSE,
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
cost matrix. Costs should be entered as a 2x2 matrix, where the first row corresponds to the true positive and true negative costs and the second row to the false positive and false negative costs. Default cost values (\code{costs=NULL}, when also \code{R=NULL}) are a combination of costs that yields R=1, which is equivalent to the Youden index method (for details about this concept, see \bold{References})
}
  \item{R}{
if the cost matrix \code{costs} is not set (that is, \code{costs=NULL}), R desired (the algorithm will choose a suitable combination of costs that leads to \code{R}). Default, \code{NULL} (which leads to R=1 using the default \code{costs}).
}
  \item{lambda}{
smoothing parameter for the NNE algorithm used in \code{survivalROC()} function
}
  \item{plot}{
Should the cost function be plotted? Default, \code{FALSE}
}
  \item{ci}{
Should a confidence interval be calculated? Default, \code{FALSE}
}
  \item{alpha}{
significance level for the confidence interval. Default, 0.05
}
  \item{B}{
number of bootstrap resamples for the confidence interval. Default, 1000
}
  \item{\dots}{
further arguments to be passed to \code{plot()}
}
}
\details{
This function estimates optimal thresholds by constructing the ROC curve at time t through time-dependent ROC curves (Heagerty et al, 2000).

Confidence intervals are obtained using normal and percentile bootstrap. In normal boostrap, the bootstrap is used to obtain an estimate of the standard error of the threshold estimate, and then the standard normal distribution is used for the confidence interval calculation. In percentile boostrap, \code{B} bootstrap resamples are generated and the threshold is estimated in all of them. Then, the confidence interval is calculated as the empirical 1-\code{alpha}/2 and 1+\code{alpha}/2 percentiles of the \code{B} bootstrap estimates.
}
\value{
An object of class \code{th_survivalROC}, which is a list with the following components:

\item{T }{A list of four elements:

  \code{thres}  threshold estimate.
  
  \code{prev}  disease prevalence used.
  
  \code{costs}  cost matrix.
  
  \code{R}  R term, the product of the non-disease odds and the cost ratio (for further details about this concept, see \bold{References}).
  
  }
  
\item{CI }{A list of five elements (or \code{NULL} if \code{ci=FALSE}):

  \code{lower.norm}  the lower limit of the confidence interval using normal bootstrap.

  \code{upper.norm}  the upper limit of the confidence interval using normal bootstrap.
  
  \code{se}  the standard error.
  
  \code{lower.perc}  the lower limit of the confidence interval using percentile bootstrap.

  \code{upper.perc}  the upper limit of the confidence interval using percentile bootstrap.
  
  \code{alpha}  significance level provided by the user.
  
  \code{ci.method} method used for the confidence intervals calculation.
  
  \code{B} number of bootstrap resamples.
  
  }
}
\references{
Heagerty PJ, Lumley T, Pepe MS. Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic Marker. Biometrics 2000; 56(2): 337-344. doi: 10.1111/j.0006-341X.2000.00337.x

Skaltsa K, Jover L, Carrasco JL. Estimation of the diagnostic threshold accounting for decision costs and sampling uncertainty. Biometrical Journal 2010; 52(5): 676-697. doi: 10.1002/bimj.200900294

Heagerty PJ, Saha-Chaudhuri P (2013). survivalROC: Time-dependent ROC curve estimation from censored survival data. R package version 1.0.3. \url{https://CRAN.R-project.org/package=survivalROC}
}
\examples{
\donttest{
with(NSCLC, th_survivalROC(COL, OS, ST, 1095,
            plot=TRUE, ci=TRUE, B=500, xlab="Collagen"))
}
}
