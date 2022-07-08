\name{NSCLC}
\alias{NSCLC}
\docType{data}
\title{
Non-small cell lung cancer (NSCLC) data
}
\description{
Non-small cell lung cancer (NSCLC) is the most common lung cancer and comprises several subtypes of lung cancers. These data come from a study by Alcaraz \emph{et al.}, 2019, in which the authors investigated the prognostic value of some activation markers in NSCLC.
}
\usage{data("NSCLC")}
\format{
  A data frame with 203 observations on the following 4 variables.
  \describe{
    \item{\code{ID}}{subject's identifier}
    \item{\code{OS}}{overall survival, that is, the time from surgery until death or last follow-up, in days}
    \item{\code{ST}}{censoring status (0=censored, 1=dead)}
    \item{\code{COL}}{percentage of collagen quantified using an imaging technique from tumour samples}
  }
}
\source{
Alcaraz J, Carrasco JL, Millares L, et al. Stromal markers of activated tumor associated fibroblasts predict poor survival and are associated with necrosis in non-small cell lung cancer. Lung Cancer 2019; 135: 151 - 160. doi: 10.1016/j.lungcan.2019.07.020
}

\examples{
data(NSCLC)
summary(NSCLC)
}
\keyword{datasets}
