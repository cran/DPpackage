\name{DPrandom}
\alias{DPrandom}
\alias{DPrandom.default}

\title{Extracts Random Effects}
\usage{
DPrandom(object,centered=FALSE,predictive=FALSE)
}
\arguments{
 \item{object}{DP fitted model object from which random effects
               estimates can be extracted.}

 \item{centered}{logical variable indicating whether the random effects
                 should be extracted centered, \code{bi}, or 
                 uncentered \code{thetai}.}

\item{predictive}{logical variable indicating whether actual or predictive
                  information of the random effects should be extracted.}                 
                 
}
\description{
  This generic function extracts Random Effects' information from DPpackage model
  objects.  
}

\examples{
\dontrun{
    # School Girls Data Example

      data(schoolgirls)
      attach(schoolgirls)

    # Prior information
    # Prior information

      tinv<-diag(10,2)
      prior<-list(alpha=1,nu0=4.01,tau1=0.001,tau2=0.001,
      tinv=tinv,mub=rep(0,2),Sb=diag(1000,2))

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-5000
      nsave<-25000
      nskip<-20
      ndisplay<-1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Fit the model
    
      fit1<-DPlmm(fixed=height~1,random=~age|child,prior=prior,mcmc=mcmc,
                  state=state,status=TRUE)
      fit1


    # Extract random effects
    
      DPrandom(fit1)
      DPrandom(fit1,centered=TRUE)
      
      plot(DPrandom(fit1))
      plot(DPrandom(fit1,centered=TRUE))

    # Extract predictive information of random effects
    
      DPrandom(fit1,predictive=TRUE)
      DPrandom(fit1,centered=TRUE,predictive=TRUE)
      
      plot(DPrandom(fit1,predictive=TRUE))
      plot(DPrandom(fit1,centered=TRUE,predictive=TRUE))
}
}
\author{

Alejandro Jara \email{<atjara@uc.cl>}

}
\keyword{models}
\keyword{nonparametric}