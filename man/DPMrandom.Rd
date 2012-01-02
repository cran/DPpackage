\name{DPMrandom}
\alias{DPMrandom}
\alias{DPMrandom.default}

\title{Extracts Random Effects}
\description{
    Extracts random effects from DPpackage objects: DPMlmm, DPMolmm, and DPMglmm.
}

\usage{
DPMrandom(object,centered=FALSE,predictive=FALSE,
          ngrid=1000,gridl=NULL)
}
\arguments{
 \item{object}{DPM fitted model object from which random effects
               estimates can be extracted.}

 \item{centered}{logical variable indicating whether the random effects
                 should be extracted centered, \code{bi}, or 
                 uncentered \code{thetai}.}

\item{predictive}{logical variable indicating whether actual or predictive
                  information of the random effects should be extracted.}                 

\item{ngrid}{number of grid points where the density estimate is 
             evaluated. This is only used if dimension of the random effects
             is lower or equal than 2. The default value is 1000.}

\item{gridl}{The limits of the interval or rectangle covered by the grid as 
            c(xl,xu) or c(xl, xu, yl, yu), respectively. If not specified the grid
            is defined automatically. This is only used if dimension of the random effects
             is lower or equal than 2 and if predictive=TRUE.}
                 
}

\examples{
\dontrun{

    # School Girls Data Example

      data(schoolgirls)
      attach(schoolgirls)

    # Prior information

      prior<-list(alpha=1,
                  tau1=0.01,tau2=0.01,
                  nu0=4.01,
                  tinv=diag(10,2),
                  nub=4.01,
                  tbinv=diag(10,2),
                  mb=rep(0,2),
                  Sb=diag(1000,2))

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-5000
      nsave<-10000
      nskip<-20
      ndisplay<-1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay)

    # Fitting the model
    
      fit1 <- DPMlmm(fixed=height~1,random=~age|child,
                     prior=prior,mcmc=mcmc,
                     state=state,status=TRUE)
      fit1

    # Extract random effects
    
      DPMrandom(fit1)
      DPMrandom(fit1,centered=TRUE)
      
      plot(DPMrandom(fit1))
      plot(DPMrandom(fit1,centered=TRUE))

    # Extract predictive information of random effects
    
      DPMrandom(fit1,predictive=TRUE)
      plot(DPMrandom(fit1,predictive=TRUE,gridl=c(75,89,3.8,7.5)))

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}
\keyword{models}
\keyword{nonparametric}