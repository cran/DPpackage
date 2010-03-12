\name{DPcaterpillar}
\alias{DPcaterpillar}

\title{Caterpillar Plots for Random Effects}
\usage{
DPcaterpillar(object, midpoint="mean", hpd=TRUE , ask=TRUE, 
              nfigr=1, nfigc=1, ...)
}
\arguments{
 \item{object}{DPrandom object from which random effects
               estimates can be extracted.}

 \item{midpoint}{variable indicating whether the mean or median 
                 of the posterior distribution of random effects
                 should be considered as "midpoint" in the caterpillar
                 plot.}

\item{hpd}{logical variable indicating whether the hpd (TRUE) or 
           pd (FALSE) of random effects should be considered in the
           caterpillar plot.}                 
           
\item{ask}{logical variable indicating whether the caterpillar plots
           should be display gradually (TRUE) or not (FALSE).}                 
                 
\item{nfigr}{integer variable indicating the number of caterpillar plots
             by row.}                 
                 
\item{nfigc}{integer variable indicating the number of caterpillar plots
             by column.}                 

\item{...}{further arguments passed to or from other methods. }             
                 
}
\description{
  This generic function produces Caterpillar Plots for Random Effects from 
  DPrandom objects.  
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
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay)

    # Fit the model
    
      fit1<-DPlmm(fixed=height~1,random=~age|child,prior=prior,
                  mcmc=mcmc,state=state,status=TRUE)
      fit1


    # Extract random effects
    
      DPcaterpillar(DPrandom(fit1))
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{plot}
\keyword{nonparametric}