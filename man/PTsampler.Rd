\name{PTsampler}
\alias{PTsampler}
\alias{PTsampler.default}

\title{Polya Tree sampler function}
\description{
    This function allows a user to generate a sample from a user-defined unormalized continuos 
    distribution using the Polya tree sampler algorithm.
}
  
\usage{

PTsampler(ltarget,dim.theta,mcmc=NULL,support=NULL,pts.options=NULL,
	  status=TRUE,state=NULL)
}

\arguments{
    \item{ltarget}{  a function giving the log of the target density.}

    \item{dim.theta}{ an integer indicating the dimension of the target density.}

    \item{mcmc}{     an optional list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nsave} giving
                     the total number of scans to be saved, and \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out). Default values are 1000, 1000, and 100 for \code{nburn}, \code{nsave},
		     and \code{ndisplay}, respectively.}   

    \item{support}{  an optional matrix, of dimension dim.theta * npoints, giving the
		     initial support points. By default the function generates 400  
		     support points from a dim.theta normal distribution with mean 0 and
		     diagonal covariance matrix with 1000 in the diagonal.}

    \item{pts.options}{ an optional list of giving the parameters needed for the PTsampler 
		     algorithm. The list must include: \code{nlevel} (an integer giving the number
		     of levels of the finite Polya tree approximation; default=5), 
                     \code{tune1} (a double precision variable representing the standard deviation of the 
                     log-normal candidate distribution for the precision parameter of the Polya tree; default=1), 		     \code{delta} (a double precision number indicating the maximum distance between the
                     target and the approximation; default=0.2), \code{max.warmup} (an integer giving the 
                     maximum number of steps allowed for the warm-up phase; default=50000), \code{minc} 
 		     (a double precision variable giving the minimum value allowed for the precision
                     parameter of the Polya tree approximation; default=1), \code{cpar0} (a double
		     precision variable giving the initial value for the precision parameter of the
                     Polya tree approximation; default=1000), and \code{nadd} (an integer variable
    		     giving the number of warm-up steps after convergence; default=1000).}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{state}{    a list giving the starting points for the MCMC algorithm. The list must
		     include: \code{theta} (a vector of dimension dim.theta of parameters),
		     \code{u} (a Polya tree decomposition matrix), \code{uinv} (a matrix giving
                     the inverse of the decompositon matrix), \code{cpar} (giving the value of the
		     Polya tree precision parameter), \code{support} (a matrix giving the final support
                     points), \code{dim.theta} (an integer giving the dimension of the problem), 
		     and \code{L1} (a double precision number giving the final convergence criterion value).}       
}

\details{
  PTsampler produces a sample from a user-defined multivariate distribution using the Polya tree sampler
  algorithm. The algorithm constructs an independent proposal based on an approximation of the target density. 
  The approximation is built from a set of support points and the predictive density of a finite  multivariate Polya tree. In an initial warm-up phase, the support points are iteratively relocated to 
  regions of higher support under the target distribution to minimize the distance between 
  the target distribution and the Polya tree predictive distribution. In the sampling phase, samples from the
  final approximating mixture of finite Polya trees are used as candidates which are accepted with a standard
  Metropolis-Hastings acceptance probability. We refer to Hanson, Monteiro, and Jara (2011) for more details
  on the Polya tree sampler.
}

\value{
  An object of class \code{PTsampler} representing the MCMC sampler. Generic functions such as \code{print},
  \code{plot}, and \code{summary} have methods to show the results of the fit. 
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. 

  The object \code{thetsave} in the output list \code{save.state} contains the samples from the target density. 
}

\seealso{
\code{\link{PTdensity}}
}

\references{
Hanson, T., Monteiro, J.V.D, and Jara, A. (2011) The Polya Tree Sampler: Toward Efficient and Automatic  Independent Metropolis-Hastings Proposals. Journal of Computational and Graphical Statistics, 20: 41-62.
}

\examples{
\dontrun{

###############################
# EXAMPLE 1 (Dog Bowl)
###############################

# Target density

  target <- function(x,y)
  {
     out <- (-3/2)*log(2*pi)-0.5*(sqrt(x^2+y^2)-10)^2-
            0.5*log(x^2+y^2)
     exp(out)
  }	

  ltarget <- function(x)
  {
     out <- -0.5*((sqrt(x[1]^2+x[2]^2)-10)^2)-
             0.5*log(x[1]^2+x[2]^2)
     out
  }	

# MCMC

  mcmc <- list(nburn=5000,
               nsave=10000,
               ndisplay=500)

# Initial support points (optional)

  support <- cbind(rnorm(300,15,1),rnorm(300,15,1))

# Scanning the posterior

  fit <- PTsampler(ltarget,dim.theta=2,mcmc=mcmc,support=support)

  fit
  summary(fit)
  plot(fit,ask=FALSE)	

# Samples saved in 
# fit$save.state$thetasave
# Here is an example of how to use them
	
  par(mfrow=c(1,2))
  plot(acf(fit$save.state$thetasave[,1],lag=100))
  plot(acf(fit$save.state$thetasave[,1],lag=100))
	  
# Plotting resulting support points

  x1 <- seq(-15,15,0.2)	
  x2 <- seq(-15,15,0.2)	
  z <- outer(x1,x2,FUN="target")
  par(mfrow=c(1,1))
  image(x1,x2,z,xlab=expression(theta[1]),ylab=expression(theta[2]))
  points(fit$state$support,pch=19,cex=0.25)

# Plotting the samples from the target density

  par(mfrow=c(1,1))
  image(x1,x2,z,xlab=expression(theta[1]),ylab=expression(theta[2]))
  points(fit$save.state$thetasave,pch=19,cex=0.25)

# Re-starting the chain from the last sample

  state <- fit$state
  fit <- PTsampler(ltarget,dim.theta=2,mcmc=mcmc,
                   state=state,status=FALSE)


###############################
# EXAMPLE 2 (Ping Pong Paddle)
###############################

  bivnorm1 <- function(x1,x2)
  {
       eval <- (x1)^2+(x2)^2 
       logDET <-  0
       logPDF <- -(2*log(2*pi)+logDET+eval)/2
       out <- exp(logPDF)
       out
  }

  bivnorm2 <- function(x1,x2)
  {
       mu <- c(-3,-3)
       sigmaInv <- matrix(c(5.263158,-4.736842,
                           -4.736842,5.263158),
                            nrow=2,ncol=2)
       eval <- (x1-mu[1])^2*sigmaInv[1,1]+
               2*(x1-mu[1])*(x2-mu[2])*sigmaInv[1,2]+
               (x2-mu[2])^2*sigmaInv[2,2] 
       logDET <-  -1.660731
       logPDF <- -(2*log(2*pi)+logDET+eval)/2
       out <- exp(logPDF)
       out
  }

  bivnorm3 <- function(x1,x2)
  {
       mu <- c(2,2)
       sigmaInv <- matrix(c(5.263158,4.736842,
                            4.736842,5.263158),
                            nrow=2,ncol=2)
       eval <- (x1-mu[1])^2*sigmaInv[1,1]+
               2*(x1-mu[1])*(x2-mu[2])*sigmaInv[1,2]+
               (x2-mu[2])^2*sigmaInv[2,2] 
       logDET <-  -1.660731
       logPDF <- -(2*log(2*pi)+logDET+eval)/2
       out <- exp(logPDF)
       out
  }

  target <- function(x,y)
  {
     out <- 0.34*bivnorm1(x,y)+
	    0.33*bivnorm2(x,y)+
	    0.33*bivnorm3(x,y)
     out
  }	

  ltarget <- function(theta)
  {
     out <- 0.34*bivnorm1(x1=theta[1],x2=theta[2])+
	    0.33*bivnorm2(x1=theta[1],x2=theta[2])+
	    0.33*bivnorm3(x1=theta[1],x2=theta[2])
     log(out)
  }	


# MCMC

  mcmc <- list(nburn=5000,
               nsave=10000,
               ndisplay=500)

# Initial support points (optional)

  support <- cbind(rnorm(300,6,1),rnorm(300,6,1))

# Scanning the posterior

  fit <- PTsampler(ltarget,dim.theta=2,mcmc=mcmc,support=support)

  fit
  summary(fit)
  plot(fit,ask=FALSE)	

# Samples saved in 
# fit$save.state$thetasave
# Here is an example of how to use them
	
  par(mfrow=c(1,2))
  plot(acf(fit$save.state$thetasave[,1],lag=100))
  plot(acf(fit$save.state$thetasave[,1],lag=100))
	  
# Plotting resulting support points

  x1 <- seq(-6,6,0.05)	
  x2 <- seq(-6,6,0.05)	
  z <- outer(x1,x2,FUN="target")
  par(mfrow=c(1,1))
  image(x1,x2,z,xlab=expression(theta[1]),ylab=expression(theta[2]))
  points(fit$state$support,pch=19,cex=0.25)

# Plotting the samples from the target density

  par(mfrow=c(1,1))
  image(x1,x2,z,xlab=expression(theta[1]),ylab=expression(theta[2]))
  points(fit$save.state$thetasave,pch=19,cex=0.25)


}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{models}
\keyword{nonparametric}