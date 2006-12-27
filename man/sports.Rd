\name{sports}
\alias{sports}
\docType{data}
\title{The Australian Athletes Data}
\description{
 This data set consider information from 202 elite Australian athletes who trained
 at the Australian Institute of Sport. The members of the sample participate in a 
 number of different sports and are about equally split between men and women.
}
\usage{data(sports)}
\format{
  A data frame with 202 observations on the following 13 variables.
  \describe{
    \item{\code{gender}}{a factor with levels \code{female} and \code{male}}
    \item{\code{sport}}{a factor with levels \code{B_Ball}, \code{Field}, \code{Gym}, \code{Netball},
                        \code{Row}, \code{Swim}, \code{T_400m}, \code{T_Sprnt}, \code{Tennis},
                        and \code{W_Polo}}
    \item{\code{rcc}}{a numeric vector giving the red cell count}
    \item{\code{wcc}}{a numeric vector giving the white cell count}
    \item{\code{Hc}}{a numeric vector giving the Hematocrit}
    \item{\code{Hg}}{a numeric vector giving the Hemoglobin level}
    \item{\code{Fe}}{a numeric vector giving the plasma ferritin concentration}
    \item{\code{bmi}}{a numeric vector giving the body mass index, weight/(height)**2}
    \item{\code{ssf}}{a numeric vector giving the sum of skin folds}
    \item{\code{Bfat}}{a numeric vector giving the body fat percentage}
    \item{\code{lbm}}{a numeric vector giving the lean body mass}
    \item{\code{Ht}}{a numeric vector giving the height (cm)}
    \item{\code{Wt}}{a numeric vector giving the weight (Kg)}
  }
}

\source{
  Cook and Weisberg (1994), An Introduction to Regression Graphics. John Wiley & Sons, New York.
}

\references{
  Cook and Weisberg (1994), An Introduction to Regression Graphics. John Wiley & Sons, New York. 
}
\examples{
data(sports)
## maybe str(sports) ; plot(sports) ...
}
\keyword{datasets}
