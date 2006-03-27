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
  
}

\examples{
\dontrun{
    # School Girls Data Example

      data(schoolgirls)
      attach(schoolgirls)

    # Prior information

      tinv<-diag(100,2)
      prior<-list(a0=0.01,b0=0.01,nu0=2,tau1=0.001,tau2=0.001,
      tinv=tinv,mub=rep(0,2),Sb=diag(10000,2))

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-100
      nsave<-50
      nskip<-10
      ndisplay<-10
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
\keyword{models}
