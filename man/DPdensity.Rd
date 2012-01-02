\name{DPdensity}
\alias{DPdensity}
\alias{DPdensity.default}

\title{Semiparametric Bayesian density estimation using a DPM of normals}
\description{
    This function generates a posterior density sample for a 
    Dirichlet process mixture of normals model.
}
  
\usage{

DPdensity(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,
          method="neal",data=sys.frame(sys.parent()),
          na.action=na.fail)      
}

\arguments{
    \item{y}{        a vector or matrix giving the data from which the density estimate 
                     is to be computed.}

    \item{ngrid}{    number of grid points where the density estimate is 
                     evaluated. This is only used if dimension of \code{y}
                     is lower or equal than 2. The default value is 1000.}

    \item{grid}{     matrix of dimension ngrid*nvar of grid points where the density estimate is 
                     evaluated. This is only used if dimension of \code{y}
                     is lower or equal than 2. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing, see details
                     below), \code{nu2} and \code{psiinv2} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix, \code{Psi1}, 
                     of the inverted Wishart part of the baseline distribution,
                     \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     gamma prior distribution of the scale parameter \code{k0} of the normal
                     part of the baseline distribution, \code{m2} and \code{s2} giving the mean
                     and the covariance of the normal prior for the mean, \code{m1}, of the normal 
                     component of the baseline distribution, respectively, \code{nu1} and 
                     \code{psiinv1} (it must be specified if \code{nu2} is missing, see details
                     below) giving the hyperparameters of the inverted Wishart part of the baseline
                     distribution and, \code{m1} giving the mean of the normal part of the baseline 
                     distribution (it must be specified if \code{m2} is missing,
                     see details below) and, \code{k0} giving the scale parameter of the normal part
                     of the baseline distribution (it must be specified if \code{tau1} is missing,
                     see details below).}

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

    \item{method}{   the method to be used. See \code{Details}.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{DPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a Dirichlet process mixture of normal model for density
  estimation (Escobar and West, 1995):
  \deqn{y_i | \mu_i, \Sigma_i \sim N(\mu_i,\Sigma_i), i=1,\ldots,n}{yi | mui, Sigmai ~ N(mui,Sigmai), i=1,\ldots,n}
  \deqn{(\mu_i,\Sigma_i) | G \sim G}{(mui,Sigmai) | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
  
  where, the baseline distribution is the conjugate normal-inverted-Wishart,
  \deqn{G_0 = N(\mu| m_1, (1/k_0) \Sigma) IW (\Sigma | \nu_1, \psi_1)}{G0 = N(mu| m1, (1/k0) Sigma) IW (Sigma | nu1, psi1)} 
  
  To complete the model specification, independent hyperpriors are assumed (optional),
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{m_1 | m_2, s_2 \sim N(m_2,s_2)}{m1 | m2, s2 ~ N(m2,s2)}
  \deqn{k_0 | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{k0 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}
  \deqn{\psi_1 | \nu_2, \psi_2 \sim IW(\nu_2,\psi_2)}{psi1 | nu2, psi2 ~ IW(nu2,psi2)}

  Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
  \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.

  To let part of the baseline distribution fixed at a particular
  value, set the corresponding hyperparameters of the prior distributions to NULL 
  in the hyperprior specification of the model.

  Although the baseline distribution, \eqn{G_0}{G0}, is a conjugate prior in
  this model specification, the algorithms with auxiliary parameters
  described in MacEachern and Muller (1998) and Neal (2000) are
  adopted. Specifically, the no-gaps algorithm of 
  MacEachern and Muller (1998), \code{"no-gaps"}, and the algorithm 8 with \eqn{m=1} 
  of Neal (2000), \code{"neal"}, are considered in the \code{DPdensity} function. 
  The default method is the algorithm 8 of Neal.
}

\value{
  An object of class \code{DPdensity} representing the DP mixture of normals
  model fit. Generic functions such as \code{print}, \code{summary}, and \code{plot} have methods to 
  show the results of the fit. The results include the baseline parameters, \code{alpha}, and the 
  number of clusters.
  
  The function \code{DPrandom} can be used to extract the posterior mean of the 
  subject-specific means and covariance matrices.

  The MCMC samples of the parameters and the errors in the model are stored in the object 
  \code{thetasave} and \code{randsave}, respectively. Both objects are included in the 
  list \code{save.state} and are matrices which can be analyzed directly by functions 
  provided by the coda package.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{an integer giving the number of clusters.} 
  
  \item{muclus}{ a matrix of dimension (nobservations+100)*(nvariables) giving the means of the clusters 
  (only the first \code{ncluster} are considered to start the chain).} 
  
  \item{sigmaclus}{ a matrix of dimension (nobservations+100)*( (nvariables)*((nvariables)+1)/2) giving 
  the lower matrix of the covariance matrix of the clusters (only the first \code{ncluster} are 
  considered to start the chain).} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each observation belongs.}

  \item{alpha}{ giving the value of the precision parameter.} 
  
  \item{m1}{ giving the mean of the normal components of the baseline distribution.}
  
  \item{k0}{ giving the scale parameter of the normal part of the baseline distribution.}
  
  \item{psi1}{giving the scale matrix of the inverted-Wishart part of the baseline distribution.}
}

\seealso{
\code{\link{DPrandom}}, \code{\link{PTdensity}}, \code{\link{BDPdensity}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

MacEachern, S. N. and Muller, P. (1998) Estimating mixture of Dirichlet Process
  Models. Journal of Computational and Graphical Statistics, 7 (2): 223-338.

Neal, R. M. (2000). Markov Chain sampling methods for Dirichlet process mixture models.
  Journal of Computational and Graphical Statistics, 9: 249-265.

}

\examples{
\dontrun{
    ####################################
    # Univariate example
    ####################################

    # Data
      data(galaxy)
      galaxy <- data.frame(galaxy,speeds=galaxy$speed/1000) 
      attach(galaxy)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 1000
      nsave <- 10000
      nskip <- 10
      ndisplay <- 100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Example of Prior information 1
    # Fixing alpha, m1, and Psi1

      prior1 <- list(alpha=1,m1=rep(0,1),psiinv1=diag(0.5,1),nu1=4,
                     tau1=1,tau2=100)


    # Example of Prior information 2
    # Fixing alpha and m1

      prior2 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.5,1)),
                     nu1=4,nu2=4,tau1=1,tau2=100)


    # Example of Prior information 3
    # Fixing only alpha

      prior3 <- list(alpha=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)


    # Example of Prior information 4
    # Everything is random

      prior4 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)

    # Fit the models

      fit1.1 <- DPdensity(y=speeds,prior=prior1,mcmc=mcmc,
                          state=state,status=TRUE)
      fit1.2 <- DPdensity(y=speeds,prior=prior2,mcmc=mcmc,
                          state=state,status=TRUE)
      fit1.3 <- DPdensity(y=speeds,prior=prior3,mcmc=mcmc,
                          state=state,status=TRUE)
      fit1.4 <- DPdensity(y=speeds,prior=prior4,mcmc=mcmc,
                          state=state,status=TRUE)

    # Posterior means
      fit1.1
      fit1.2
      fit1.3
      fit1.4

    # Plot the estimated density
      plot(fit1.1,ask=FALSE)
      plot(fit1.2,ask=FALSE)
      plot(fit1.3,ask=FALSE)
      plot(fit1.4,ask=FALSE)

    # Extracting the density estimate
      cbind(fit1.1$x1,fit1.1$dens)
      cbind(fit1.2$x1,fit1.2$dens)
      cbind(fit1.3$x1,fit1.3$dens)
      cbind(fit1.4$x1,fit1.4$dens)
      
    # Plot the parameters (only prior 2 for illustration)
    # (to see the plots gradually set ask=TRUE)
      plot(fit1.2,ask=FALSE,output="param")

    # Plot the a specific parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1.2,ask=FALSE,output="param",param="psi1-speeds",
           nfigr=1,nfigc=2)

    # Extracting the posterior mean of the specific 
    # means and covariance matrices 
    # (only prior 2 for illustration)
      DPrandom(fit1.2) 

    # Ploting predictive information about the specific 
    # means and covariance matrices 
    # with HPD and Credibility intervals
    # (only prior 2 for illustration)
    # (to see the plots gradually set ask=TRUE)
      plot(DPrandom(fit1.2,predictive=TRUE),ask=FALSE)
      plot(DPrandom(fit1.2,predictive=TRUE),ask=FALSE,hpd=FALSE)

    # Ploting information about all the specific means 
    # and covariance matrices 
    # with HPD and Credibility intervals
    # (only prior 2 for illustration)
    # (to see the plots gradually set ask=TRUE)
      plot(DPrandom(fit1.2),ask=FALSE,hpd=FALSE)


    ####################################
    # Bivariate example
    ####################################

    # Data
      data(airquality)
      attach(airquality)

      ozone <- Ozone**(1/3)
      radiation <- Solar.R

    # Prior information

      s2 <- matrix(c(10000,0,0,1),ncol=2)
      m2 <- c(180,3)
      psiinv2 <- solve(matrix(c(10000,0,0,1),ncol=2))
     
      prior <- list(a0=1,b0=1/5,nu1=4,nu2=4,s2=s2,
                    m2=m2,psiinv2=psiinv2,tau1=0.01,tau2=0.01)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 10000
      nskip <- 10
      ndisplay <- 1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Fit the model
      fit1 <- DPdensity(y=cbind(radiation,ozone),prior=prior,mcmc=mcmc,
                        state=state,status=TRUE,na.action=na.omit)

    # Plot the estimated density
      plot(fit1)

    # Extracting the density estimate
      fit1$x1
      fit1$x2
      fit1$dens
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
