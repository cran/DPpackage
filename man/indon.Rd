\name{indon}
\alias{indon}
\docType{data}
\title{Indonesian Children's Health Study}
\description{
  This data set consider respiratory infection information of 250 indonesian
  children reported by Sommer, Katz, and Tarwotjo (1984). The children, 
  all preschoolers, were seen quarterly for up to six quarters. At each examination, 
  the presence or absence of respiratory infection was noted. 
}

\usage{data(indon)}

\format{
  A data frame with 1200 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{an ordered factor giving a unique identifier for the subject in the study.}
    \item{\code{gender}}{a numeric vector giving the gender.}
    \item{\code{height}}{a numeric vector giving the height for age as a percentage
                         of the National Center for Health Statitics standard
                         centered at 90.}
    \item{\code{cosv}}{a numeric vector giving the seasonal cosine for the annual cycle.}
    \item{\code{sinv}}{a numeric vector giving the seasonal sine for the annual cycle.}
    \item{\code{xero}}{a numeric vector giving the presence (1) or absence (0) of
                       xeropthalmia.}
    \item{\code{baseage}}{a numeric vector giving the age at the entry.}
    \item{\code{age}}{a numeric vector giving the age of the child in months centered at 36}
    \item{\code{infect}}{a numeric vector giving the presence (1) or absence (0) of
                         respiratory infection.}
  }                         
}

\source{
  Sommer, A., Katz, J., and Tarwotjo, I. (1984) Increased risk of respiratory infection and
     diarrhea in children with pre-existing mild vitamin A deficiency, American Journal
     of Clinical Nutrition, 40: 1090-1095.
}

\references{
  Zeger, S.L., and Karim, M.R. (1991) Generalized linear models with random effects: A Gibbs
     sampling approach. Journal of the American Statistical Association, 86: 79-86.
}
\examples{
data(indon)
## maybe str(indon) ; plot(indon) ...
}
\keyword{datasets}
