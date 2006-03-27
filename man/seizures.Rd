\name{seizures}
\alias{seizures}
\docType{data}
\title{Epileptic seizures}
\description{
  This data set consider information from a clinical trial of 59 epileptics,
  reported by Thall and Vail (1990). For each patient, the number of epileptic
  seizures was recorded during a baseline period of eight weeks. Patients were
  then randomized to treatment with the anti-epileptic drug progabide, or to
  aplacebo in addition to standard chemotherapy. The number of seizures was then
  recorded in four consecutive two-weeks intervals.
}

\usage{data(seizures)}

\format{
  A data frame with 295 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{an ordered factor giving a unique identifier for the subject in the study.}
    \item{\code{seize}}{a numeric vector giving the number of epileptic seizures.}
    \item{\code{visit}}{a numeric vector giving the number of the visit, 0=baseline, and
                        1,2,3, and 4 for the four consecutive two-weeks intervals.}
    \item{\code{trt}}{a numeric vector giving the treatment group.}
    \item{\code{age}}{a numeric vector giving the age at the entry.}
  }
}

\source{
  Thall, P.F., and Vail, S.C. (1990) Some covariance models for longitudinal count data with
     ovserdispersion, Biometrics, 46: 657-671.
}
\references{
  Diggle, P.J., Liang, K-Y., and Zeger, S.L. (1994) Analysis of longitudinal data. Oxford:
     Clarendon Press.
}
\examples{
data(seizures)
## maybe str(seizures) ; plot(seizures) ...
}
\keyword{datasets}
