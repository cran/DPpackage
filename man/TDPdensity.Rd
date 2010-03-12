\name{TDPdensity}
\alias{TDPdensity}
\alias{TDPdensity.default}

\title{Semiparametric Bayesian density estimation using DP Mixtures of Triangular Distributions}
\description{
    This function generates a posterior density sample for a 
    Triangular-Dirichlet model.
}
  
\usage{

TDPdensity(y,support=3,transform=1,ngrid=1000,prior,mcmc,state,status,
          data=sys.frame(sys.parent()),na.action=na.fail)      
      
}

\arguments{
    \item{y}{        a vector giving the data from which the density estimate 
                     is to be computed.}

    \item{support}{  an integer number giving the support of the random density,
                     1=[0,1], 2=(0, +Inf], and 3=(-In,+Inf). Depending on this,
                     the data is transformed to lie in the [0,1] interval.}

    \item{transform}{  an integer number giving the type of transformation
                      to be considered, 1=Uniform, 2=Normal,3=Logistic,4=Cauchy.
                      The types 2-4 can be only used when the support is the  
                      real line.}

    \item{ngrid}{    number of grid points where the density estimate is 
                     evaluated. This is only used if dimension of \code{y}
                     is lower or equal than 2. The default value is 1000.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{aa0} and \code{ab0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{aa0} is missing, see details
                     below), \code{a0} and \code{b0} giving the parameters of the
                     beta centering distribution of the DP prior, and
                     \code{kmax} giving the maximum value of the discrete uniform
                     prior for number of components in the Mixture of Triangular distributions.
                     Optionally, when the support of the data is the real line and the 
                     parametric transformation 2-4 are considered, the location \code{mu}
                     and the scale parameter \code{sigma2} can be included here. If not,
                     they are taked as the mean and the variance of the data, respectively.                     }

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
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{TDPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a Triangular-Dirichlet model for density
  estimation:
  \deqn{y_i | G \sim G, i=1,\ldots,n}{yi | G ~ G, i=1,\ldots,n}
  \deqn{G | kmax, \alpha, G_0 \sim TDP(kmax,\alpha G_0)}{G | kmax, alpha, G0 ~ TDP(kmax, alpha G0)}
  
  where, \eqn{y_i}{yi} is the transformed data to lie in [0,1], \code{kmax} is the upper limit of the discrete uniform prior 
  for the number of components in the Mixture of Triangular
  distributions, \eqn{\alpha}{alpha} is the total mass parameter of the Dirichlet process component,
  and \eqn{G_0}{G0} is the centering distribution of the DP. The centering distribution corresponds
  to a \eqn{G_0=Beta(a_0,b_0)}{G0=Beta(a0,b0)} distribution.
  
  Note that our representation is different to the Mixture of Triangular 
  distributions proposed by Perron and Mengersen (2001). In this function
  we consider random weights following a Dirichlet prior and we exploit the 
  underlying DP structure. By so doing, we avoid using Reversible-Jumps algorithms. 
  
  The precision or total mass parameter, \eqn{\alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.
  
}

\value{
  An object of class \code{TDPdensity} representing the Triangular-Dirichlet 
  model fit. Generic functions such as \code{print}, \code{summary}, and \code{plot} have methods to 
  show the results of the fit. The results include the degree of the polynomial \code{k}, \code{alpha}, and the 
  number of clusters.
  
  The MCMC samples of the parameters and the errors in the model are stored in the object 
  \code{thetasave} and \code{randsave}, respectively. Both objects are included in the 
  list \code{save.state} and are matrices which can be analyzed directly by functions 
  provided by the coda package.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{an integer giving the number of clusters.} 
  
  \item{yclus}{ a real vector giving the \code{y} latent variables of the clusters (only the first \code{ncluster} are 
  considered to start the chain).} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each observation belongs.}

  \item{alpha}{ giving the value of the precision parameter.} 
  
  \item{k}{ giving the number of components in the Mixture of Triangular distriutions.}
}

\seealso{
\code{\link{DPdensity}}, \code{\link{PTdensity}}, \code{\link{BDPdensity}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.


Perron, F. and Mengersen, K. (2001) Bayesian Nonparametric Modeling Using Mixtures 
  of Triangular Distributions. Biometrics, 57(2): 518-528.

}

\examples{
\dontrun{
    # Data
      data(galaxy)
      galaxy<-data.frame(galaxy,speeds=galaxy$speed/1000) 
      attach(galaxy)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-1000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Prior
      prior<-list(aa0=2.01,
                  ab0=0.01,
                  kmax=50,
                  a0=1,
                  b0=1)

    # Fitting the model

      fit<-TDPdensity(y=speeds,prior=prior,mcmc=mcmc,state=state,status=TRUE)
      
      plot(fit)
      

}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}