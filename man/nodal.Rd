\name{nodal}
\alias{nodal}
\docType{data}
\title{Nodal Involvement Data}
\description{
  This data set consider information on the presence of
  prostatic nodal involvement collected on 53 patients with 
  prostate cancer reported by Brown (1980).
  
  For the sample of prostate cancer patients, a number of possible 
  predictor variables were measured before surgery. The patients then had 
  surgery to determine nodal involvement. It was required to see if nodal 
  involvement could be accurately predicted from the predictor variables 
  and which ones were most important. 
}
\usage{data(nodal)}
\format{
  A data frame with 53 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{an ordered factor giving a unique identifier for the subject in the study}
    \item{\code{ssln}}{a numeric vector giving the prostatic nodal involvement that
                       takes the value 1 if cancer had spread to the surrounding lymph
                       nodes and 0 otherwise}
    \item{\code{age}}{a numeric vector giving the age of the patient in years at
                      diagnosis}
    \item{\code{acid}}{a numeric vector giving the level of serum acid phosphate}
    \item{\code{xray}}{a numeric vector giving the result af an X-ray examination,
                       coded 0 if negative and 1 if positive}
    \item{\code{size}}{a numeric vector giving the size of the tumor, coded 
                       0 if small and 1 if large}
    \item{\code{grade}}{a numeric vector giving the pathological grade of the tumor, coded 
                       0 if less serious and 1 if more serious}
  }
}

\source{

Brown, B.W. (1980) Prediction analysis for binary data. In Biostatistics Casebook.
    R.G. Miller, B. Efron, B.W. Brown and L.E. Moses (editors), 3-18. John Wiley. 

}
\references{

Chib, S. (1995) Marginal Likelihood from the Gibbs output. Journal
   of the American Statistical Association, 90: 1313 - 1321.
}
\examples{
\dontrun{
    # Data
      data(nodal)
      attach(nodal)
      lacid<-log(acid)

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn<-5000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                   tune=1.1)


    # Prior distribution
      prior <- list(alpha=1,beta0=c(0,rep(0.75,5)),
                    Sbeta0=diag(c(100,rep(25,5)),6))

    # Fit the model
      fit1 <- DPbinary(ssln~age+lacid+xray+size+grade,prior=prior,mcmc=mcmc,
                       state=state,status=TRUE) 
      fit1

}

}
\keyword{datasets}
