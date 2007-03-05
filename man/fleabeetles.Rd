\name{fleabeetles}
\alias{fleabeetles}
\docType{data}
\title{Flea-beetles}
\description{
   This data set consider physical information of 74 male flea-beetles reported by Lubischew 
   (1962). Information of three species (Ch. concinna, Ch. heptapotamica, and Ch. heikertingeri) 
   is considered and 6 measurements on each flea-beetles.
  
}
\usage{data(fleabeetles)}
\format{
  A data frame with 74 observations on the following 7 variables.
  \describe{
    \item{\code{fjft}}{a numeric vector giving the width of the first joint of 
                       the first tarsus in microns (the sum of measurements 
                       for both tarsi)}
    \item{\code{sjft}}{a numeric vector giving the width of the second joint of 
                       the first tarsus in microns (the sum of measurements 
                       for both tarsi)}
    \item{\code{mwhbee}}{a numeric vector giving the maximal width of the head 
                         between the external edges of the eyes in 0.01 mm}
    \item{\code{mwafp}}{a numeric vector giving the maximal width of the aedeagus 
                        in the fore-part in microns}
    \item{\code{faa}}{a numeric vector giving the front angle of the aedeagus 
                      (1 unit = 7.5 degrees)}
    \item{\code{awfs}}{a numeric vector giving the aedeagus width from the side 
                       in microns}
    \item{\code{species}}{a numeric vector giving the species: 1=Ch. concinna,
                          2= Ch. heptapotamica, and 3=Ch. heikertingeri}
  }
}

\source{
Lubischew, A. A. (1962) On the Use of Discriminant Functions in Taxonomy, 
   Biometrics, 18: 455-477.  
}

\references{
MacEachern, S.N., and Muller, P. (1998) Estimating Mixture of Dirichlet Process Models,
   Journal of Computational and Graphical Statistics, 7: 223-238.
}
\examples{
data(fleabeetles)
## maybe str(fleabeetles) ; plot(fleabeetles) ...
}
\keyword{datasets}
