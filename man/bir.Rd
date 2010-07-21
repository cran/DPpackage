\name{bir}
\alias{bir}
\docType{data}
\title{British Institute of Radiology Fractionation Studies}
\description{
  The British Institute of Radiology (BIR) conducted two large-scale
  randomized clinical trials to assess the effectiveness of different
  radiotherapy treatment schedules for cancer of the larynx and pharynx.
  The cambined data come from 858 subjects with laryngeal squamous cell
  carcinomas and no involvement of asjacent organs. These data have been
  described and analyzed by Rezvani, Fowler, Hopewell, and Alcock (1993).
}
\usage{data(bir)}
\format{
  A data frame with 858 observations on the following 6 variables.
  \describe{
    \item{\code{response}}{Three-year local control 1 (0 no control).}
    \item{\code{dose}}{Total dose (grays).}
    \item{\code{df}}{Total dose x dose/fraction.}
    \item{\code{time}}{Total time of treatment (days).}
    \item{\code{kt2}}{Tumor status (indicators for 2nd level of factor).}
    \item{\code{kt3}}{Tumor status (indicators for 3rd level of factor).}
  }
}
\details{
  Three-year local control - meaning no detection of laryngeal carcinoma 
  within three years after treatment - is the binary response, coded as 1 
  if local control is achieved and 0 otherwise. For this data set, three-year
  local control is achieved for 69% of the cases. In a treatment schedule, 
  a total dose of radiation is administered in fractions over a treatment
  period. The dose per fraction df is measured in grays (Gy), the length
  of treatment period time is measured in days, and the number of fractions
  of the dose is nf. Tumors are classified by stage (i.e., the extent of 
  invasion), into three groups. This categorical covariate is coded by two
  indicator variables kt2 and kt3, which are defined by kt2=1 (kt3=1) is the
  tumor is stage II (stage III) and zero otherwise.
  
  Chappell, Nondahl and Fowler (1995) argued that the tumor stage, the total dose, the 
  total time, and the interaction of the total dose per fraction are the 
  relevant explanantory variables affecting probability of local control.
}
\references{

  Chappell, R., Nondahl, D.M., and Fowler, J.F. (1995) Modelling Dose and
  Local Control in Radiotheraphy. Journal of the American Statistical Asso-
  ciation, 90: 829 - 838.

  Newton, M.A., Czado, C., and Chappell, R. (1996) Bayesian inference 
  for semiparametric binary regression. Journal of the American Statistical 
  Association, 91, 142-153.
  
  Rezvani, M., Fowler, J., Hopewell, J., and Alcock, C. (1993) 
  Sensitivity of Human Squamous Cell Carcinoma of the Larynx to Fractionated
  Radiotherapy. British Journal of Radiology, 66: 245 - 255.
}
\examples{
data(bir)
## maybe str(bir) ; plot(bir) ...

}
\keyword{datasets}
