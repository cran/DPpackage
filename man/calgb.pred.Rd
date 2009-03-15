\name{calgb.pred}
\alias{calgb.pred}
\docType{data}
\title{Cancer and Leukemia Group B (CALGB) Data for Prediction}
\description{
  Same as 'calgb' for future patients (for prediction).
}
\usage{data(calgb.pred)}
\format{
  A data frame with 8 observations on the following 10 variables.
  \describe{
    \item{\code{Z1}}{a numeric vector giving the estimated Z1 coefficients 
                     of the logistic regression curve.}
    \item{\code{Z2}}{a numeric vector giving the estimated Z2 coefficients 
                     of the logistic regression curve.}
    \item{\code{Z3}}{a numeric vector giving the estimated Z3 coefficients 
                     of the logistic regression curve.}
    \item{\code{T1}}{a numeric vector giving the estimated time point 
                     where the horizontal line of the curve is defined,
                     i.e., the curve consists of a horizontal line up to
                     t=T1ji.}
    \item{\code{T2}}{a numeric vector giving the estimated time point where
                     the logistic component of the curve is defined, i.e.,
                     the curve consist of a logistic regression curve
                     starting at t=T2ji.}
    \item{\code{B0}}{a numeric vector giving the estimated B0 coefficients 
                     of the logistic regression curve.}
    \item{\code{B1}}{a numeric vector giving the estimated B1 coefficients 
                     of the logistic regression curve.}
    \item{\code{CTX}}{a numeric vector giving the dose level of
                      cyclophosphamide.}
    \item{\code{GM}}{a numeric vector giving the dose level GM-CSF.}
    \item{\code{AMOF}}{a numeric vector giving the dose level of
                       amifostine.}
  }
}
\examples{
data(calgb.pred)
## maybe str(calgb.pred) ; plot(calgb.pred) ...
}
\keyword{datasets}
