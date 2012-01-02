\name{BDPdensity}
\alias{BDPdensity}
\alias{BDPdensity.default}

\title{Semiparametric Bayesian density estimation using Bernstein Polynomials}
\description{
    This function generates a posterior density sample for a 
    Bernstein-Dirichlet model.
}
  
\usage{

BDPdensity(y,support=3,ngrid=1000,grid=NULL,prior,mcmc,state,status,
          data=sys.frame(sys.parent()),na.action=na.fail)      
      
}

\arguments{
    \item{y}{        a vector giving the data from which the density estimate 
                     is to be computed.}

    \item{support}{  an integer number giving the support of the random density,
                     1=[0,1], 2=(0, +Inf], and 3=(-In,+Inf). Depending on this,
                     the data is transformed to lie in the [0,1] interval.}

    \item{ngrid}{    number of grid points where the density estimate is 
                     evaluated. This is only used if dimension of \code{y}
                     is lower or equal than 2. The default value is 1000.}

    \item{grid}{     vector of grid points where the density estimate is 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{aa0} and \code{ab0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{aa0} is missing, see details
                     below), \code{a0} and \code{b0} giving the parameters of the
                     beta centering distribution of the DP prior, and
                     \code{kmax} giving the maximum value of the discrete uniform
                     prior for the degree of the Bernstein polynomial.}

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
                     \code{BDPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a Bernstein-Dirichlet model for density
  estimation (Petrone, 1999a, 1999b; Petrone and Waserman, 2002):
  \deqn{y_i | G \sim G, i=1,\ldots,n}{yi | G ~ G, i=1,\ldots,n}
  \deqn{G | kmax, \alpha, G_0 \sim BDP(kmax,\alpha G_0)}{G | kmax, alpha, G0 ~ BDP(kmax, alpha G0)}
  
  where, \eqn{y_i}{yi} is the transformed data to lie in [0,1], \code{kmax} is the upper limit of the discrete uniform prior for the degree of the Bernstein
  polynomial, \eqn{\alpha}{alpha} is the total mass parameter of the Dirichlet process component,
  and \eqn{G_0}{G0} is the centering distribution of the DP. The centering distribution corresponds
  to a \eqn{G_0=Beta(a_0,b_0)}{G0=Beta(a0,b0)} distribution.
  
  The precision or total mass parameter, \eqn{\alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.
  
}

\value{
  An object of class \code{BDPdensity} representing the Bernstein-Dirichlet 
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
  
  \item{k}{ giving the degree of the Bernstein polynomial.}
}

\seealso{
\code{\link{DPdensity}}, \code{\link{PTdensity}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Petrone, S. (1999a) Random Bernstein Polynomials. Scandinavian Journal of 
  Statistics, 26: 373-393.

Petrone, S. (1999b) Bayesian density estimation using Bernstein polynomials. 
  The Canadian Journal of Statistics, 27: 105-126.

Petrone, S. and Waserman, L. (2002) Consistency of Bernstein polynomial posterior. 
  Journal of the Royal Statistical Society, Series B, 64: 79-100.

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
                  kmax=1000,
                  a0=1,
                  b0=1)

    # Fitting the model

      fit <- BDPdensity(y=speeds,prior=prior,mcmc=mcmc,
                        state=state,status=TRUE)
      
      plot(fit)
      

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

Fernando Quintana \email{<quintana@mat.puc.cl>}

}

\keyword{models}
\keyword{nonparametric}