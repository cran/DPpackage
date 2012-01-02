\name{DPpsBF}
\alias{DPpsBF}

\title{Computes Pseudo Bayes Factors from DPpackage output}
\description{
This function computes Pseudo Bayes Factors from DPpackage
output.

}
\usage{
DPpsBF(...)
}

\arguments{
  \item{...}{DPpackage output objects. These have to be of the same class.} 
}

\examples{
\dontrun{
    # Respiratory Data Example

      data(indon)
      attach(indon)

      baseage2 <- baseage**2
      follow <- age-baseage
      follow2 <- follow**2 

    # Prior information

      beta0 <- rep(0,9)
      Sbeta0 <- diag(1000,9)
      tinv <- diag(1,1)
      prior <- list(a0=2,b0=0.1,nu0=4,tinv=tinv,
                    mub=rep(0,1),Sb=diag(1000,1),
                    beta0=beta0,Sbeta0=Sbeta0)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5
      nsave <- 100
      nskip <- 5
      ndisplay <- 100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Fit the Probit model
      fit1 <- DPglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+
                     baseage2+follow+follow2,random=~1|id,
                     family=binomial(probit),
                     prior=prior,mcmc=mcmc,state=state,status=TRUE)

    # Fit the Logit model
      fit2 <- DPglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+
                     baseage2+follow+follow2,random=~1|id,
                     family=binomial(logit),
                     prior=prior,mcmc=mcmc,state=state,status=TRUE)

    # Model comparison
      DPpsBF(fit1,fit2)


}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}
\keyword{models}
\keyword{nonparametric}