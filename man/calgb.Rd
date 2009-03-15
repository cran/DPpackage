\name{calgb}
\alias{calgb}
\docType{data}
\title{Cancer and Leukemia Group B (CALGB) Data}
\description{
 Data from two studies carried out by the Cancer and Leukemia Group B
 (CLAGB): CALGB 8881 and CALGB 9160. In both studies, the main response was white blood cell count (WBC) 
 for each patient over time. Mueller and Rosner (1998) used a non-linear patient  
 specific regression model. The data consider the subject-specific 
 regression parameters (Z1, Z2, Z3, T1, T2, B0, B1) and information on
 covariates.

 CALGB has kindly agreed to make these data available for interested readers, subject to the following 
 conditions: i) Any paper using these data should acknowledge CALGB for the use of the data, and ii)
 the paper should reference the original papers describing the studies.
}
\usage{data(calgb)}
\format{
  A data frame with 98 observations on the following 12 variables.
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
    \item{\code{pat}}{a numeric vector giving the patient indicators.}
    \item{\code{study}}{a numeric vector giving the study indicators.}
  }
}
\source{
CALGB 8881:
Lichtman, S. M., Ratain, M. J., Echo, D. A., Rosner, G., Egorin, M. J., Budman, D. R., 
  Vogelzang,N. J., Norton, L. and Schilsky, R. L. (1993) Phase I trial and granulocyte-macrophage 
  colony-stimulating factor plus high-dose cyclophosphamide given every 2 weeks: a Cancer and 
  Leukemia Group B study. Journal of the National Cancer Institute, 85: 1319-1326.

CALGB 9160:
Budman, D., Rosner, G., Lichtman, S., Miller, A., Ratain, M. and Schilsky, R. (1998) 
   A randomized trial of wr-2721 (amifostine) as a chemoprotective agent in combination with 
   high-dose cyclophosphamide and molgramostim (GM-CSG). Cancer Therapeutics, 1: 164-167.
}
\references{
Mueller, P. and Rosner, G. (1998). Semiparametric PK/PD Models. In: Practical Nonparametric 
  and Semiparametric Bayesian Statistics, Eds: D. Dey, P. Muller, D. Sinha, New York: 
  Springer-Verlag, pp. 323-337.

Mueller, P., Quintana, F. and Rosner, G. (2004). A Method for Combining Inference 
  over Related Nonparametric Bayesian Models. Journal of the Royal Statistical 
  Society, Series B, 66: 735-749.
}
\examples{
data(calgb)
## maybe str(calgb) ; plot(calgb) ...
}
\keyword{datasets}
