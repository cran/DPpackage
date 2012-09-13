\name{igg}
\alias{igg}
\docType{data}
\title{
Immunoglobulin G concentrations
}
\description{
  This data set consider information on the serum immunoglbulin
  G (IgG) concentration from 298 children aged 6 months to 6 years old.
  The data were reported by Isaacs et al. (1983). These data were 
  further analysed by Royston and Wright (1998) and Kapitula and 
  Bedrick (2005) using the parametric exponential normal family, which includes  parameters for skew and kurtosis that can be functions of covariates.
}
\usage{data(igg)}
\format{
  A data frame with 298 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{a numeric vector giving the age of the children 
                      (in years).}
    \item{\code{igg}}{a numeric vector giving the the serum 
                      IgG concentration.}
  }
}
\source{
Isaacs, D., Altman, D. G., Tidmarsh, C. E., Valman, H. B., Webster, 
    A. D. B. (1983). Serum immunoglobulin concentrations in preschool 
    children measured by laser nephelometry: reference ranges for IgG, 
    IgA, IgM. J. Clin. Pathol. 36: 1193 - 1196.
}
\references{
Kapitula, L. R., Bedrick, E. J. (2005). Diagnostics for the exponential 
    normal growth curve model. Statist. Med. 24: 95 - 108.

Royston, P., Wright, E. M. (1998). A method for estimation age-specific 
    reference intervals ('normal ranges') based on fractional polynomials 
    and exponential transformation. J. R. Statist. Soc. A 161: 79 - 101.
}
\examples{
data(igg)
## maybe str(igg) ; plot(igg) ...
}
\keyword{datasets}
