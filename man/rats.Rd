\name{rats}
\alias{rats}
\docType{data}
\title{Rats}
\description{
 This example is taken from section 6 of Gelfand and Smith (1990), and concerns 
 30 young rats whose weights were measured weekly for five weeks. 
}
\usage{data(rats)}
\format{
  A data frame with 150 observations on the following 3 variables.
  \describe{
    \item{\code{weight}}{a numeric vector giving the weight of the rat}
    \item{\code{day}}{a numeric vector giving the day of the weight evaluation}
    \item{\code{rat}}{an ordered factor giving a unique identifier for the subject in the study}
  }
}
\source{
 Gelfand, A.E. (with S. Hills, A. Racine-Poon and A.F.M. Smith) 1990. Illustration of 
 Bayesian Inference in Normal Data Models Using Gibbs Sampling. Journal Amer. Stat. 
 Assoc., 85, 972-985.
  
}
\references{
 Gelfand, A.E. (with S. Hills, A. Racine-Poon and A.F.M. Smith) 1990. Illustration of 
 Bayesian Inference in Normal Data Models Using Gibbs Sampling. Journal Amer. Stat. 
 Assoc., 85, 972-985.
 }
\examples{
data(rats)
## maybe str(rats) ; plot(rats) ...
}
\keyword{datasets}
