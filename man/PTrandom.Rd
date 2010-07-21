\name{PTrandom}
\alias{PTrandom}
\alias{PTrandom.default}

\title{Extracts Random Effects}
\description{
    Extracts random effects from PTglmm objects: PTlmm, PTolmm, and PTglmm.
}

\usage{
PTrandom(object,centered=FALSE,predictive=FALSE,ngrid=1000,gridl=NULL)
}
\arguments{
 \item{object}{PT fitted model object from which random effects
               estimates can be extracted.}

 \item{centered}{logical variable indicating whether the random effects
                 should be extracted centered, \code{bi}, or 
                 uncentered \code{thetai}. This option cannot be only used
                 to get the density estimates.}

\item{predictive}{logical variable indicating whether actual or predictive
                  information of the random effects should be extracted.}                 

\item{ngrid}{number of grid points where the density estimate is 
             evaluated. This is only used if dimension of the random effects
             is lower or equal than 2 and if predictive=TRUE. The default value is 1000.}
             
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
                  M=4,
                  tau1=0.01,tau2=0.01,
                  nu0=4.01,
                  tinv=diag(10,2),
                  mub=rep(0,2),
                  Sb=diag(1000,2))
 
    # Initial state
      state <- NULL
 
    # MCMC parameters
 
      nburn<-5000
      nsave<-5000
      nskip<-0
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                    tune1=1.5,tune2=1.1)
 
    # Fitting the model
     
      fit1<-PTlmm(fixed=height~1,random=~age|child,prior=prior,mcmc=mcmc,
                  state=state,status=TRUE)
      fit1
 
 
    # Extract random effects
    
      PTrandom(fit1)
      
      plot(PTrandom(fit1))

    # Extract predictive information of random effects
    
      PTrandom(fit1,predictive=TRUE)
      plot(PTrandom(fit1,predictive=TRUE,gridl=c(75,89,3.8,7.5)))
}
}
\author{
Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{models}
\keyword{nonparametric}