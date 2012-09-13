\name{DPMdencens}
\alias{DPMdencens}
\alias{DPMdencens.default}

\title{Bayesian density estimation for interval-censored data using a DPM of normals}
\description{
    This function generates a posterior density sample for a 
    Dirichlet process mixture of normals model for interval-censored data.
}
  
\usage{
DPMdencens(left,right,ngrid=100,grid=NULL,prior,mcmc,state,status)
}

\arguments{
    \item{left}{     a vector or matrix giving the lower limit for each response variable.
                     Note that the responses are defined on the entire real line and that
                     unknown limits should be indicated by NA.}

    \item{right}{    a vector or matrix giving the upper limit for each response variable.
                     Note that the responses are defined on the entire real line and that
                     unknown limits should be indicated by NA.}

    \item{ngrid}{    number of grid points where the density estimate is 
                     evaluated. The default value is 100.}

    \item{grid}{     matrix of dimension ngrid*nvar of grid points where  the density estimate is 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             interval limits.}
    
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
}

\details{
  This generic function fits a Dirichlet process mixture of normal model for density
  estimation (Escobar and West, 1995) based on interval-censored data: 
  \deqn{y_{ij} \in [l_{ij},u_{ij}), i=1,\ldots,n, j=1,\ldots,m}{yij in [lij,uij),i=1,\ldots,n, j=1,\ldots,m,}
  \deqn{y_i | \mu_i, \Sigma_i \sim N(\mu_i,\Sigma_i), i=1,\ldots,n,}{yi | mui, Sigmai ~ N(mui,Sigmai), i=1,\ldots,n,}
  \deqn{(\mu_i,\Sigma_i) | G \sim G,}{(mui,Sigmai) | G ~ G,}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0),}{G | alpha, G0 ~ DP(alpha G0),}
  where, \eqn{y_i=(y_{i1},\ldots,y_{im})}{yi=(yi1,\ldots,yim)}, and the baseline distribution is the conjugate 
  normal-inverted-Wishart distribution,
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
  this model specification, an algorithm based on auxiliary parameters is adopted. 
  Specifically, the algorithm 8 with \eqn{m=1} of Neal (2000) is considered in the \code{DPMdencens} function. 

  Finally, note that this function can be used to fit the DPM of normals model for ordinal data
  proposed by Kottas, Mueller and Quintana (2005). In this case, the arbitrary cut-off points must
  be specified in \code{left} and \code{right}. Samples from the predictive distribution contained
  in the (last columns) of the object randsave (please see below) can be used to obtain an estimate of 
  the cell probabilities. 
}

\value{
  An object of class \code{DPMdencens} representing the DP mixture of normals
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

  \item{y}{giving the matrix of imputed data points.}
}

\seealso{
\code{\link{DPrandom}}, \code{\link{DPdensity}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90:
  577-588.

Kottas, A., Mueller, P., Quintana, F. (2005). Nonparametric Bayesian
  modeling for multivariate ordinal data. Journal of Computational and 
  Graphical Statistics, 14: 610-625.

Neal, R. M. (2000). Markov Chain sampling methods for Dirichlet process
  mixture models. Journal of Computational and Graphical Statistics, 9:
  249-265.

}

\examples{
\dontrun{
    ####################################
    # Bivariate example:
    # Censored data is artificially
    # created 
    ####################################

      data(airquality)      
      attach(airquality)

      ozone <- Ozone**(1/3)
      radiation <- Solar.R
      y <- na.omit(cbind(radiation,ozone))

    # create censored-data 
      xxlim <- seq(0,300,50)
      yylim <- seq(1.5,5.5,1)

      left <- matrix(0,nrow=nrow(y),ncol=2)
      right <- matrix(0,nrow=nrow(y),ncol=2)

      for(i in 1:nrow(y))
      {
          left[i,1] <- NA
          right[i,1] <- NA
          if(y[i,1] < xxlim[1]) right[i,1] <- xxlim[1] 
          for(j in 1:length(xxlim))
          { 
              if(y[i,1] >= xxlim[j]) left[i,1] <- xxlim[j]	 
              if(y[i,1] >= xxlim[j]) right[i,1] <- xxlim[j+1]	 
          }	 
          left[i,2] <- NA
          right[i,2] <- NA
          if(y[i,2] < yylim[1]) right[i,2] <- yylim[1] 

          for(j in 1:length(yylim))
          { 
              if(y[i,2] >= yylim[j]) left[i,2] <- yylim[j]	 
              if(y[i,2] >= yylim[j]) right[i,2] <- yylim[j+1]	 
          }	 
      }

    # Prior information
      s2 <- matrix(c(10000,0,0,1),ncol=2)
      m2 <- c(180,3)
      psiinv2 <- diag(c(1/10000,1),2)

      prior <- list(alpha=1,nu1=4,nu2=4,s2=s2,
                    m2=m2,psiinv2=psiinv2,tau1=0.01,tau2=0.01)

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn <- 5000
      nsave <- 5000
      nskip <- 3
      ndisplay <- 1000
      mcmc <- list(nburn=nburn, 
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fitting the model
      fit1 <- DPMdencens(left=left,right=right,ngrid=100,
                         prior=prior,mcmc=mcmc,
                         state=state,status=TRUE)
      fit1
      summary(fit1)

    # Plot the estimated density
      plot(fit1)

    # Extracting the univariate density estimates
      cbind(fit1$grid[,1],fit1$funi[[1]])
      cbind(fit1$grid[,2],fit1$funi[[2]])
       
    # Extracting the bivariate density estimates
      fit1$grid[,1]
      fit1$grid[,2]
      fit1$fbiv[[1]]

    # Plot of the estimated density along with the 
    # true data points and censoring limits
      contour(fit1$grid[,1],fit1$grid[,2],fit1$fbiv[[1]])
      points(y)
      abline(v=xxlim)
      abline(h=yylim)
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
