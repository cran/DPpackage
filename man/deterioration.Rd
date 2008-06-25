\name{deterioration}
\alias{deterioration}
\docType{data}
\title{Time to Cosmetic Deterioration of Breast Cancer Patients}
\description{
  This data set considers information of time to cosmetic deterioration of the breast for 
  women with Stage 1 breast cancer who have undergone a lumpectomy, for two treatments, these being
  radiation, and radiation coupled with chemotherapy. There is interest in the cosmetic impact of 
  the treatments because both are considered very effective in preventing recurrence of this early 
  stage cancer. The data come from a retrospective study of 46 patients who received radiation 
  only and 48 who received radiation plus chemotherapy. Each woman made a series of visits to a 
  clinician, who determined whether or not retraction had occurred. If it had, the time of retraction 
  was known only to lie between the time of the present and last visits. The data set is presented in 
  Beadle et al. (1984a,b) and also given in Finkelstein and Wolfe (1985).
}
\usage{data(deterioration)}
\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{left}}{a numeric vector giving the left limit of the interval}
    \item{\code{right}}{a numeric vector giving the right limit of the interval}
    \item{\code{trt}}{a numeric vector giving the treatment (0 = radiation only, 
                      1 = radiation plus chemotherapy)}
  }
}

\source{
Beadle, G., Come, S., Henderson, C., Silver, B., and Hellman, S. (1984a). The effect of
  adjuvant chemotherapy on the cosmetic results after primary radiation treatment for early
  stage breast cancer. International Journal of Radiation Oncology, Biology and Physics, 
  10: 2131-2137.
   
Beadle, G., Harris, J., Silver, B., Botnick, L., and Hellman, S. (1984b). Cosmetic results
  following primary radiation therapy for early breast cancer. Cancer, 54: 2911-2918.

Finkelstein, D.M. and Wolfe, R.A. (1985). A semiparametric model for regression analysis of
  interval-censored failure time data. Biometrics, 41: 933-945.
}

\references{
Hanson, T., and Johnson, W. (2004) A Bayesian Semiparametric AFT Model for
  Interval-Censored Data. Journal of Computational and Graphical Statistics, 
  13: 341-361.
}
\examples{
data(deterioration)
## maybe str(deterioration) ; plot(deterioration) ...
}
\keyword{datasets}
