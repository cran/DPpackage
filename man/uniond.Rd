\name{uniond}
\alias{uniond}
\docType{data}
\title{Union Membership}
\description{
  This data set consider growth information on wages and union membership 
  for 534 workers. The datafile contains observations on 11 variables 
  sampled from the Current Population Survey of 1985.  This data set 
  demonstrates multiple regression, confounding, transformations, 
  multicollinearity, categorical variables, ANOVA, pooled tests of 
  significance, interactions and model building strategies.
}
\usage{data(uniond)}
\format{
  A data frame with 534 observations on the following 11 variables.
  \describe{
    \item{\code{education}}{a numeric vector giving the number of 
                            years of education.}
    \item{\code{south}}{a numeric vector gving an indicator variable for 
                        Southern Region (1=Person lives in South, 
                        0=Person lives elsewhere).}
    \item{\code{sex}}{a numeric vector giving an indicator variable 
                      for sex (1=Female, 0=Male).}
    \item{\code{experience}}{a numeric vector giving the number of years of 
                             work experience.}
    \item{\code{unionv}}{a numeric vector giving an indicator variable for union 
                        membership (1=Union member, 0=Not union member).}
    \item{\code{wage}}{a numeric vector giving the Wage (dollars per hour).}
    \item{\code{age}}{a numeric vector giving the Age in years.}
    \item{\code{race}}{a numeric vector giving the race (1=Other,
                       2=Hispanic, 3=White).}
    \item{\code{occupation}}{a numeric vector giving the occupational 
                             category (1=Management, 2=Sales, 3=Clerical, 
                             4=Service, 5=Professional, 6=Other).}    
   \item{\code{sector}}{a numeric vector giving the Sector (0=Other, 
                        1=Manufacturing, 2=Construction).}
    \item{\code{marr}}{a numeric vector giving the Marital Status 
                       (0=Unmarried,  1=Married).}
  }
}

\details{
  The Current Population Survey (CPS) is used to supplement census information 
  between census years. These data consist of a random sample of 534 persons 
  from the CPS, with information on wages and other characteristics of the 
  workers, including sex, number of years of education, years of work experience,
  occupational status, region of residence and union membership. We wish to 
  determine (i) whether wages are related to these characteristics and (ii) 
  whether there is a gender gap in wages. Based on residual plots, wages 
  were log-transformed to stabilize the variance. Age and work experience were 
  almost perfectly correlated (r=.98). Multiple regression of log wages against 
  sex, age, years of education, work experience, union membership, southern 
  residence, and occupational status showed that these covariates were related 
  to wages (pooled F test, p < .0001). The effect of age was not significant 
  after controlling for experience. Standardized residual plots showed no 
  patterns, except for one large outlier with lower wages than expected. 
  This was a male, with 22 years of experience and 12 years of education, 
  in a management position, who lived in the north and was not a union member. 
  Removing this person from the analysis did not substantially change the 
  results, so that the final model included the entire sample. 
  Adjusting for all other variables in the model, females earned 81% (75%, 88%)
  the wages of males (p < .0001). Wages increased 41% (28%, 56%) for every 5 
  additional years of education (p < .0001). They increased by 11% (7%, 14%) 
  for every additional 10 years of experience (p < .0001). Union members were 
  paid 23% (12%, 36%) more than non-union members (p < .0001). Northerns were 
  paid 11% (2%, 20%) more than southerns (p =.016). Management and professional
  positions were paid most, and service and clerical positions were paid least 
  (pooled F-test, p < .0001). Overall variance explained was R2 = .35. 
  In summary, many factors describe the variations in wages: occupational 
  status, years of experience, years of education, sex, union membership and 
  region of residence. However, despite adjustment for all factors that were 
  available, there still appeared to be a gender gap in wages. There is no 
  readily available explanation for this gender gap.
}

\examples{
data(uniond)
## maybe str(uniond) ; plot(uniond) ...
}
\keyword{datasets}
