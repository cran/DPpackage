\name{schoolgirls}
\alias{schoolgirls}
\docType{data}
\title{The Heights of Schoolgirls}
\description{
  This data set consider growth information of 20 preadolescent schoolgirls reported by Goldstein 
  (1979, Table 4.3, p. 101). The height of girls was measured on a yearly basis from age 6 to 10. 
  The measurements are given at exact years of age, some having been previously adjusted to these. 
  Further, the girls were classified according to the height of their mother into three categories: 
  short mothers, medium mothers and tall mothers. 
}
\usage{data(schoolgirls)}
\format{
  A data frame with 100 observations on the following 4 variables.
  \describe{
    \item{\code{height}}{a numeric vector giving the height in cm}
    \item{\code{child}}{an ordered factor giving a unique identifier for the subject in the study}    
    \item{\code{age}}{a numeric vector giving the age of the child in years}
    \item{\code{group}}{a factor with levels 1 (short), 2 (medium), and 3 (tall) giving the 
                        mother category}
  }
}
\details{
  Measurements reported by Goldstein(1979) for one of the girls (child 5) were 114.5, 112.0, 126.4, 
  131.2, and 135.0. In this data set, the second measurement was replaced by 122.0. 
}
\source{
  Goldstein, H. (1979) The Design and Analysis of Longitudinal Studies. London: Academic Press.
}
\references{
  Verbeke, G., and Molenberghs, G. (2000) Linear Mixed Models for Longitudinal Data. New York:
  Springer-Verlag.(Section 12.7)
}
\examples{
data(schoolgirls)
## maybe str(schoolgirls) ; plot(schoolgirls) ...
}
\keyword{datasets}
