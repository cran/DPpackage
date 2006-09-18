\name{PTdensity}
\alias{PTdensity}
\alias{PTdensity.default}

\title{Nonparametric Bayesian density estimation using Polya Trees}
\description{
    This function generates a posterior density sample for a 
    Mixture of Polya trees model.
}
  
  
\usage{

PTdensity(y,ngrid=1000,prior,mcmc,state,status,
          data=sys.frame(sys.parent()),na.action=na.fail)      
      
}

\arguments{
    \item{y}{        a vector giving the data from which the density estimate 
                     is to be computed.}

    \item{ngrid}{    number of grid points where the density estimate is 
                     evaluated. The default value is 1000.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Poly tree prior,
                     \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{alpha} is missing, see details
                     below) and, optionally, \code{M} giving the finite level
                     to be considered. If \code{M} is specified, a Partially Specified
                     Mixture of Polya trees is fitted.}

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
                     \code{DPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a Mixture of Polya Trees prior for the density estimation 
  (see, e.g., Lavine, 1992 and 1994; Hanson, 2006):

  \deqn{Y_1, \ldots , Y_n | G \sim G}{Y1,...,Yn | G ~ G}
  \deqn{G | \alpha,\mu,\sigma^2 \sim PT(\Pi^{\mu,\sigma^2},\textit{A})}{G | alpha,mu,sigma2 ~ PT(Pi^{mu,sigma2},\textit{A})}
  \deqn{f(\mu,\sigma^{-2}) \propto \sigma^{-2}}{f(mu,sigma^-2) \propto 1/sigma2}
  
  where, the the PT is centered around a \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution, by
  taking each \eqn{m} level of the partition \eqn{\Pi^{\mu, \sigma^2}}{Pi^{mu, sigma2}} to coincide 
  with the \eqn{k/2^m, k=0,\ldots,2^m} quantile of the \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution.
  The family \eqn{\textit{A}=\{\alpha_e: e \in E^{*}\}}{\textit{A}=\{alphae: e \in E^{*}\}}, 
  where \eqn{E^{*}=\bigcup_{m=0}^{m} E^m}{E^{*}=\bigcup_{m=0}^{M} E^m} 
  and \eqn{E^m} is the \eqn{m}-fold product of \eqn{E=\{0,1\}},
  was specified as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}. 
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
  The results include \code{mu}, \code{sigma2}, and the precision
  parameter \code{alpha}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects:
  
  \item{mu}{ giving the value of the baseline mean.} 
  
  \item{sigma}{ giving the baseline standard deviation.}
  
  \item{alpha}{ giving the value of the precision parameter.}

}

\seealso{
\code{\link{DPdensity}}
}

\references{
Hanson, T. (2006) Inference for Mixtures of Finite Polya Trees.
  To appear in Journal of the American Statistical Association. 

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
      nburn<-100
      nsave<-1000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                   tune1=0.15,tune2=1.1,tune3=1.1)

    # Prior information
      prior<-list(alpha=1,M=8)


    # Fit the model

      fit1<-PTdensity(y=speeds,ngrid=1000,prior=prior,mcmc=mcmc,
                       state=state,status=TRUE)

    # Posterior means
      fit1

    # Plot the estimated density
      plot(fit1,ask=FALSE)

    # Plot the parameters
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,output="param")

    # Extracting the density estimate
      cbind(fit1$x1,fit1$dens)

}
}

\author{
Alejandro Jara \email{<Alejandro.JaraVallejos@med.kuleuven.be>}

Tim Hanson \email{<hanson@biostat.umn.edu>}
}

\keyword{models}
