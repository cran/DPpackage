\name{hiv}
\alias{hiv}
\docType{data}
\title{
HIV-AIDS data
}
\description{
  This data set considers information from a cohort of 262 hemophiliacs at
  risk of human immunodeficiency virus (HIV) infection from infusions of   
  blood they received periodically to treat their hemophilia in two
  hospitals in France. All infected patients are believed to have become
  infected by contaminated blood factor: 105 patients received at least   
  1,000 micro grams/kg of blood factor for at least one year between 1982 and 
  1985 (heavily treated group), and 157 patients received less than 1,000  
  micro grams/kg in each year (lighter treated group). For this cohort both 
  infection with HIV and the onset of acquired immunodeficiency syndrome
  (AIDS) or other clinical symptoms could be subject to censoring. 
  Therefore, the induction time between infection and clinical AIDS are
  treated as doubly-censored.
}

\usage{data(hiv)}
\format{
  A data frame with 262 observations on the following 5 variables.
  \describe{
    \item{\code{onsetL}}{a numeric vector giving the lower limit of the 
                         HIV infection interval.}
    \item{\code{onsetU}}{a numeric vector giving the upper limit of the 
                         HIV infection interval.}
    \item{\code{failureL}}{a numeric vector giving the lower limit of the 
                           interval were clinical AIDS was observed.}
    \item{\code{failureU}}{a numeric vector giving the upper limit of the 
                           interval were clinical AIDS was observed.}
    \item{\code{trt}}{a numeric vector giving the treatment indicator: 
                      (0) indicates the lighter treated group while (1) 
                      indicates the heavily treated group.}
  }
}
\details{
This dataset was analyzed by dataset De Gruttola and Lagakos (1989).
The periodic observation of HIV infection status in these patients was possible because 
blood samples were stored and retrospectively tested for evidence of infection with the HIV. 
Note that both the distribution of chronological time of infection and induction time are of interest. 
In De Gruttola and Lagakos (1989) the proposed nonparametric maximum likelihood one-sample estimator 
was illustrated by considering the intervals for the onset and failure time, which were the results 
of a discretization of the time axis into 6-month intervals.
}
\source{
De Grutola, V. and Lagakos, S.W. (1989). Analysis of doubly-censored survival data, with
   application to AIDS. Biometrics, 45: 1-11.
}
\references{
Jara, A., Lesaffre, E., De Iorio, M., Quintana, F. (2010). Bayesian semiparametric inference 
   for multivariate doubly-interval-censored data. Annals of Applied Statistics, 4: 2126-2149.
}
\examples{
data(hiv)
## maybe str(hiv) ; plot(hiv) ...
}
\keyword{datasets}
