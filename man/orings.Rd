\name{orings}
\alias{orings}
\docType{data}
\encoding{latin1}
\title{Challenger Space Shuttle O-Ring Data}
\description{
The motivation for collecting this database was the explosion of the
USA Space Shuttle Challenger on 28 January, 1986. The Rogers commission 
concluded that the Challenger accident was caused by gas leak through 
the 6 o-ring joints of the shuttle. Dalal, Fowlkes and Hoadley (1989) 
looked at the number of distressed o-rings (among the 6) versus launch 
temperature (Temperture) and pressure (Pressure) for 23 previous 
shuttle flights. The previous shuttles were launched at temperatures 
between 53ºF and 81ºF.
}
\usage{data(orings)}
\format{
  A data frame with 138 observations on the following 4 variables.
  \describe{
    \item{ThermalDistress}{a numeric vector indicating wether the
                           o-ring experienced thermal distress}
    \item{Temperature}{a numeric vector giving the launch 
                       temperature (degrees F)}
    \item{Pressure}{a numeric vector giving the leak-check pressure (psi)}
    \item{Flight}{a numeric vector giving the temporal order of flight}
  }
}

\source{
Dalal, S.R., Fowlkes, E.B., and Hoadley, B. (1989). Risk analysis of 
  space shuttle : Pre-Challenger Prediction of Failure, Journal of the 
  American Statistical Association, 84: 945 - 957.
}
\references{
Dalal, S.R., Fowlkes, E.B., and Hoadley, B. (1989). Risk analysis of 
  space shuttle : Pre-Challenger Prediction of Failure, Journal of the 
  American Statistical Association, 84: 945 - 957.

Lavine, M. (1991). Problems in extrapolation illustrated with space 
  shuttle O-ring data.  Journal of the American Statistical Association, 
  86: 919-922.

Martz, H. F., Zimmer, W.J. (1992). The risk of catastrophic failure
  of the solid rocket boosters on the space shuttle. The American 
  Statistician, 46: 42-47. 
}
\examples{
data(orings)
## maybe str(orings) ; plot(orings) ...
}
\keyword{datasets}
