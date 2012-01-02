\name{DPmultmeta}
\alias{DPmultmeta}
\alias{DPmultmeta.default}

\title{Bayesian analysis for a semiparametric random effects multivariate meta-analysis model using a MDP}
\description{
    This function generates a posterior density sample for a 
    semiparametric random effects multivariate meta-analysis model using a Dirichlet process
    or a Mixture of Dirichlet process prior for the distribution of the 
    random effects. 
    Support provided by the NIH/NCI R01CA75981 grant.
}
  
  
\usage{

DPmultmeta(y,asymvar,prior,mcmc,state,status,
           data=sys.frame(sys.parent()))

}

\arguments{
    \item{y}{        a vector or matrix giving the data or effects from which the density estimate 
                     is to be computed.}

    \item{asymvar}{  a vactor or matrix giving the asymptotic covariance matrix for each effect.
                     The dimension of this matrix is the number of records/studies times the
                     the half-stored elements of the study-specific covariance matrix.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing, see details
                     below), \code{m2} and \code{s2} giving the mean
                     and the covariance of the normal prior for the mean, \code{m1}, of the mean of the 
                     normal baseline distribution, respectively, \code{m1} giving the mean of the baseline
                     distribution (it must be specified if \code{m2} is missing), 
                     \code{nu} and \code{psiinv} giving the hyperparameters of the inverted Wishart distribution on the 
                     covariance matrix \code{s1} of the normal baseline 
                     distribution, and  \code{s1} giving the covariance matrix of the baseline distribution 
                     (it must be specified if \code{nu} is missing).}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, and \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{data}{     data frame.}       
    
}

\details{
  This generic function fits a semiparametric random effects multivariate meta-analysis model: 

  \deqn{y_i \sim N(\mu_i , \Sigma_i), i=1,\ldots,n}{yi ~ N(thetai, Sigmai), i=1,\ldots,n}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}

  where, \eqn{G_0=N(\theta| m_1, s_1)}{G0=N(theta| m1, s1)}.  To complete the model 
  specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{m_1 | m_2, s_2 \sim N(mu_2,s_2)}{m1 | m2, s2 ~ N(m2,s2)}
  \deqn{s_1 | \nu, \psi \sim IW(\nu,\psi)}{s1 | nu, psi ~ IW(nu,psi)}

  Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
  \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.

  To let part of the baseline distribution fixed at a particular
  value, set the corresponding hyperparameters of the prior distributions to NULL 
  in the hyperprior specification of the model.

  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for conjugate priors 
  for a collapsed state of MacEachern (1998). 
 
}

\value{
  An object of class \code{DPmultmeta} representing the random
  effects model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary}, have methods to show the results of the fit. 
  The results include \code{m1}, \code{s1}, \code{alpha}, and the 
  number of clusters.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{alpha}{ giving the value of the precision parameter} 
  
  \item{muclus}{ a matrix of dimension (nobservations+1)*(nvariables) giving the means of the clusters 
  (only the first \code{ncluster} are considered to start the chain).} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{m1}{ giving the mean of the normal baseline distributions.} 
  
  \item{s1}{ giving the covariance matrix of the normal baseline distributions.}
}

\seealso{
\code{\link{DPmeta}}
}

\references{

MacEachern, S.N. (1998) Computational Methods for Mixture of Dirichlet Process Models,
  in Practical Nonparametric and Semiparametric Bayesian Statistics, 
  eds: D. Dey, P. Muller, D. Sinha, New York: Springer-Verlag, pp. 23-44.

}

\examples{
\dontrun{
    ######################################################################   
    # Simulated Data: 
    #         mu_i ~ 0.5 N(mub1,Sigmab1) + 0.5 N(mub2,Sigmab2)
    #         y_i ~ N(mu_i,Sigma_i) 
    #         Sigma_1=...=Sigma_n=Sigma assumed to be known 
    ######################################################################
      nvar <- 2
      nrec <- 100
      Sigma <- matrix(c(0.25,0.15,0.15,0.25),nrow=nvar,ncol=nvar) 
      mub1 <- rep(-1.5,nvar)  
      mub2 <- rep( 0.5,nvar)  
      Sigmab1 <- matrix(c(0.25,-0.175,-0.175,0.25),nrow=nvar,ncol=nvar)
      Sigmab2 <- matrix(c(0.25, 0.0875, 0.0875,0.25),nrow=nvar,ncol=nvar)
      
      ind <- rbinom(nrec,1,0.5)
      z1 <- mub1+matrix(rnorm(nvar*nrec),nrow=nrec,ncol=nvar)%*%chol(Sigmab1)        
      z2 <- mub2+matrix(rnorm(nvar*nrec),nrow=nrec,ncol=nvar)%*%chol(Sigmab2)       
      mu <- ind*z1+(1-ind)*z2
      
      y <- NULL
      for(i in 1:nrec)
      {
          z <- mu[i,]+matrix(rnorm(nvar),nrow=1,ncol=nvar)%*%chol(Sigma)        
          y <- rbind(y,z)

      }
      colnames(y) <- c("y1","y2")  
 
    ######################################################################   
    # Asymptotic variance
    ######################################################################   
      z <- NULL
      for(i in 1:nvar)
      {
          for(j in i:nvar)
          {
                z <- c(z,Sigma[i,j])         
          }
      }
      asymvar <- matrix(z,nrow=nrec,ncol=nvar*(nvar+1)/2,byrow=TRUE)


    # Prior information

      s2 <-diag(100,nvar)
      m2 <-rep(0,nvar)
      nu <- 4
      psiinv <- diag(1,nvar)
     
      prior<-list(a0=1,
                  b0=1/5,
                  nu=nu,
                  m2=m2,
                  s2=s2,
                  psiinv=psiinv)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 500
      nsave <- 1000
      nskip <- 0
      ndisplay <- 100
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fitting the model
      fit1 <- DPmultmeta(y=y,asymvar=asymvar,prior=prior,
                         mcmc=mcmc,state=state,status=TRUE)

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

Peter Mueller \email{<pmueller@mdanderson.org>}

}

\keyword{models}
\keyword{nonparametric}
