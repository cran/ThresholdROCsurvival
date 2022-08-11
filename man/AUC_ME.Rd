\name{AUC_ME}
\alias{AUC_ME}
\title{
AUC estimation using the missing exclusion (ME) method}
\description{
This function estimates the AUC with survival data using by excluding subjects with missing status at the time-point of interest
}
\usage{
AUC_ME(cont.var, time, status, predict.time, plot = FALSE,
       ci = TRUE, alpha = 0.05, ...)
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
  \item{plot}{
Should a graph of the ROC curve be plotted? Default, \code{FALSE}
}
  \item{ci}{
Should a confidence interval be calculated? Default, \code{TRUE}
}
  \item{alpha}{
significance level for the confidence interval. Default, 0.05
}
  \item{\dots}{
further arguments to be passed to \code{plot()}
}
}
\details{
First, the algorithm determines the status of the subjects at time \code{predict.time}. Those subjects whose status could be not be determined (therefore, with missing status) are excluded from the analysis. Then, the AUC is estimated using the \code{roc} function from \code{pROC} package (Robin et al, 2011). Confidence intervals for the AUC are calculated using the logit transformation (Kottas et al, 2014).
}
\value{
An object of class \code{AUC_ME}, which is a list with three components: 
\item{AUC }{AUC estimate}
\item{CI }{2-dimensional vector containing the confidence interval}
\item{data }{a data.frame containing the following columns previously provided by the user: \code{cont.var}, \code{time} and \code{status}, and a new column \code{statusNA}, which contains the status of the subjects at time \code{predict.time} (0=no event, 1=event, \code{NA}=missing)}
}
\references{
Kottas M, Kuss O, Zapf A. A modified Wald interval for the area under the ROC curve (AUC) in diagnostic case-control studies. BMC Medical Research Methodology 2014; 14(26). doi:10.1186/1471-2288-14-26

Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, Muller M. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 2011; 12. doi:10.1186/1471-2105-12-77
}

\seealso{
\code{\link{th_ME}}
}
\examples{
data(NSCLC)
res <- with(NSCLC, AUC_ME(COL, OS, ST, 1095, plot=FALSE))
res
res$data
}
