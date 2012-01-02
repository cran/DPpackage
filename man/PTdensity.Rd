\name{PTdensity}
\alias{PTdensity}
\alias{PTdensity.default}

\title{Nonparametric Bayesian density estimation using Mixtures of Polya Trees}
\description{
    This function generates a posterior density sample for a 
    Mixture of Polya trees model.
}
  
  
\usage{

PTdensity(y,ngrid=1000,grid=NULL,prior,mcmc,state,status,
          data=sys.frame(sys.parent()),na.action=na.fail)      
      
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
                     prior distribution of the precision parameter of the Poly tree prior,
                     \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{alpha} is missing, see details
                     below), optionally \code{M} giving the finite level
                     to be considered (if \code{M} is specified, a partially specified
                     mixture of Polya trees model is fitted), \code{nu0} and \code{tinv} or 
                     \code{tau1} and \code{tau2} giving the hyperparameters of the 
                     inverted Wishart or inverted gamma prior distribution for the centering  
                     covariance or variance, respectively, \code{sigma} giving the value of the 
                     standard deviation (univariate case) or covariance
                     matrix (multivariate case) of the centering distribution (if missing and if 
                     \code{nu0} and \code{tinv} or \code{tau1} and \code{tau2} are missing, 
                     Jeffrey's prior is used for the centering (co)variance matrix,
                     \code{m0} and \code{S0} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution, and \code{mu} giving the value of the mean of the
                     centering distribution (if missing and if  
                     \code{m0} and \code{S0} are missing, Jeffery's prior is used for \code{mu}).}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), \code{tune1}, \code{tune2}, and \code{tune3} giving the 
                     positive Metropolis tuning parameter for the baseline mean, variance, and
                     precision parameter, respectively (the default value is 1.1)}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{PTdensity} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a Mixture of Polya Trees prior for the density estimation 
  (see, e.g., Lavine, 1992 and 1994; Hanson, 2006). In the univariate case, the model 
  is given by:
  \deqn{Y_1, \ldots , Y_n | G \sim G}{Y1,...,Yn | G ~ G}
  \deqn{G | \alpha,\mu,\sigma \sim PT(\Pi^{\mu,\sigma^2},\textit{A})}{G | alpha,mu,sigma2 ~ PT(Pi^{mu,sigma2},\textit{A})}
  where, the the PT is centered around a \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution, by
  taking each \eqn{m} level of the partition \eqn{\Pi^{\mu, \sigma^2}}{Pi^{mu, sigma2}} to coincide 
  with the \eqn{k/2^m, k=0,\ldots,2^m} quantile of the \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution.
  The family \eqn{\textit{A}=\{\alpha_e: e \in E^{*}\}}{\textit{A}=\{alphae: e \in E^{*}\}}, 
  where \eqn{E^{*}=\bigcup_{m=0}^{m} E^m}{E^{*}=\bigcup_{m=0}^{M} E^m} 
  and \eqn{E^m} is the \eqn{m}-fold product of \eqn{E=\{0,1\}},
  was specified as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}. 

  Analogous to the univariate model, in the multivariate case the PT prior is 
  characterized by partitions of \eqn{R^d}{R^d}, and a collection of conditional
  probabilities that link sets in adjacent tree levels, i.e., they link each
  parent set in a given level to its \eqn{2^d}{2^d} offpring stes in the subsequent
  level. The multivariate model is given by: 
  \deqn{Y_1, \ldots , Y_n | G \sim G}{Y1,...,Yn | G ~ G}
  \deqn{G | \alpha,\mu,\Sigma \sim PT(\Pi^{\mu,\Sigma},\textit{A})}{G | alpha,mu,Sigma ~ PT(Pi^{mu,Sigma},\textit{A})}
  where, the the PT is centered around a \eqn{N_d(\mu,\Sigma)}{N_d(mu,Sigma)} distribution.
  In this case, the class of partitions that we consider, starts with base sets that are Cartesian 
  products of intervals obtained as quantiles from the standard normal distribution. A multivariate 
  location-scale transformation, \eqn{Y=\mu+\Sigma^{1/2} z}{Y=mu+Sigma^{1/2} z}, is applied to each 
  base set yielding the final sets. 

  A Jeffry's prior can be specified for the centering parameters, 
  \deqn{f(\mu,\sigma^2 \propto \sigma^{-2}}{f(mu,sigma2) \propto 1/sigma2}
  and
  \deqn{f(\mu,\Sigma) \propto |\Sigma|^{-(d+1)/2}}{p(mu,Sigma) \propto |Sigma|^{-(d+1)/2}}
  in the univariate and multivariate case, respectively. Alternatively, the centering parameters
  can be fixed to user-specified values or proper priors can be assigned. In the univariate case,
  the following proper priors can be assigned:
  \deqn{\mu | m_0, S_0 \sim N(m_0,S_0)}{mu | m0, S0 ~ N(m0,S0)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim \Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}
 
  In the multivariate case, the following proper priors can be assigned: 
  \deqn{\mu | m_0, S_0 \sim N(m_0,S_0)}{mu | m0, S0 ~ N(m0,S0)}
  \deqn{\Sigma | \nu_0, T \sim IW(\nu_0,T)}{Sigma | nu0, T ~ IW(nu0,T)}
  Note that the inverted-Wishart prior is parametrized such that
  \eqn{E(\Sigma)= T^{-1}/(\nu_0-q-1)}{E(Sigma)= Tinv/(nu0-q-1)}.

  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  
  The precision parameter, \eqn{\alpha}{alpha}, of the \code{PT} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  In the computational implementation of the model, Metropolis-Hastings steps are used to 
  sample the posterior distribution of the baseline and precision parameters.
}

\value{
  An object of class \code{PTdensity} representing the Polya tree
  model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. 
  The results include \code{mu}, \code{sigma} or \code{Sigma} in the univariate
  or multivariate case, respectively, and the precision
  parameter \code{alpha}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects:
  
  \item{mu}{ giving the value of the baseline mean.} 
  
  \item{sigma}{ giving the baseline standard deviation or the baseline 
                covariance matrix in the univariate or multivariate case, 
                respectively.}
  
  \item{alpha}{ giving the value of the precision parameter.}

}

\seealso{
\code{\link{DPdensity}}, \code{\link{BDPdensity}}
}

\references{
Hanson, T. (2006) Inference for Mixtures of Finite Polya Trees.
  Journal of the American Statistical Association, 101: 1548-1565.

Lavine, M. (1992) Some aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 20: 1222-11235.

Lavine, M. (1994) More aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 22: 1161-1176.

}

\examples{
\dontrun{
    ####################################
    # Univariate example
    ####################################

    # Data
      data(galaxy)
      galaxy<-data.frame(galaxy,speeds=galaxy$speed/1000) 
      attach(galaxy)

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn <- 2000
      nsave <- 5000
      nskip <- 49
      ndisplay <- 500
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                   tune1=0.03,tune2=0.25,tune3=1.8)

    # Prior information
      prior<-list(a0=1,b0=0.01,M=6,m0=21,S0=100,sigma=20)

    # Fitting the model

      fit1 <- PTdensity(y=speeds,ngrid=1000,prior=prior,mcmc=mcmc,
                        state=state,status=TRUE)

    # Posterior means
      fit1

    # Plot the estimated density
      plot(fit1,ask=FALSE)
      points(speeds,rep(0,length(speeds)))

    # Plot the parameters
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,output="param")

    # Extracting the density estimate
      cbind(fit1$x1,fit1$dens)


    ####################################
    # Bivariate example
    ####################################

    # Data
      data(airquality)
      attach(airquality)

      ozone <- Ozone**(1/3)
      radiation <- Solar.R

    # Prior information
      prior <- list(a0=5,b0=1,M=4,
                    m0=c(0,0),S0=diag(10000,2),
                    nu0=4,tinv=diag(1,2))

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn <- 2000
      nsave <- 5000
      nskip <- 49
      ndisplay <- 500
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                   tune1=0.8,tune2=1.0,tune3=1)

    # Fitting the model
      fit1 <- PTdensity(y=cbind(radiation,ozone),prior=prior,mcmc=mcmc,
                        state=state,status=TRUE,na.action=na.omit)

      fit1

    # Plot the estimated density
      plot(fit1)

    # Extracting the density estimate
      x1 <- fit1$x1
      x2 <- fit1$x2
      z <- fit1$dens
      par(mfrow=c(1,1))
      contour(x1,x2,z)
      points(fit1$y)  

}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{models}
\keyword{nonparametric}