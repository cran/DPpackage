\name{psychiatric}
\alias{psychiatric}
\docType{data}
\title{Psychiatric Clinical Trial}
\description{
  This data set consider information from a pychiatric clinical trial,
  reported by Hedeker and Gibbons (1994) and collected in the NIMH 
  Schizophrenia Collaborative Study on treatment-related changes in
  overall severity. For each patient, the Impatient Multidimensional
  Psychiatric Scale was scored. In this study, 437 patients were randomly
  assigned to receive one of four medications: placebo, chlorpromazine,
  fluphenazine, or thioridazine.
}
\usage{data(psychiatric)}
\format{
  A data frame with 1603 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{a numeric vector giving a unique identifier for the subject in the study}
    \item{\code{imps79o}}{a factor giving the Impatient Multidimensional Psuchiatric Scale
                          scored as \code{1} = normal or borderline mentally ill,
                          \code{2} = midly or moderately ill, \code{3} = markadely ill,
                          and \code{4} = severely or among the the most extremely ill}
    \item{\code{tx}}{a numeric vector giving the treatment group, 0=Placebo, 1=Drug}
    \item{\code{week}}{a numeric vector giving the week where the measurement ocurred}
    \item{\code{sweek}}{a numeric vector giving the square root transformation of time}
  }
}
\source{
  Hedeker, D. and Gibbons, R.D. (1994) A random-effects ordinal regression model for
     multilevel data. Biometrics, 50: 933-944.
}
\references{
  Hedeker, D. and Gibbons, R.D. (1994) A random-effects ordinal regression model for
     multilevel data. Biometrics, 50: 933-944.
}
\examples{
data(psychiatric)
## maybe str(psychiatric) ; plot(psychiatric) ...
}
\keyword{datasets}
