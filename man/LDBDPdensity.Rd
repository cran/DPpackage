\name{LDBDPdensity}
\alias{LDBDPdensity}
\alias{LDBDPdensity.default}

\title{Bounded Density Regression using Dependent Bernstein Polynomials}
\description{
    This function generates a posterior density sample for a 
    Linear Dependent Bernstein-Dirichlet Process model for bounded conditional
    density estimation.
}
   
\usage{
LDBDPdensity(formula,xpred,prior,mcmc,state,status,ngrid=100,
             grid=NULL,compute.band=FALSE,type.band="PD",
             data=sys.frame(sys.parent()),
             na.action=na.fail,work.dir=NULL)
}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. The design matrix is used to model 
                     the distribution of the response in the LDBDP model. The response
                     is assumed to take values in [0,1].} 

    \item{xpred}{    a matrix giving the covariate values where the predictive 
                     density is evaluated.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{lambda} a double precision giving the parameter of the truncated Poisson
                     prior distribution for the degree, \code{k}, of the Bernstein polynomial,
		     \code{maxn} an integer giving the truncation of the the stick-breaking approximation to 
                     the dependent Dirichlet process, \code{alpha} giving the value of the precision parameter
                     of the dependent Dirichlet process, \code{m0} and \code{S0} giving the hyperparameters of
                     the normal prior distribution for the mean, \code{mub}, of the normal baseline 
                     distribution \code{nu} and \code{psiinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix, \code{Sb}, 
                     of the baseline distribution.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following elements: \code{nburn} an integer giving the number of burn-in 
                     scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer 
                     giving the total number of scans to be saved, \code{ndisplay} an integer giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), \code{slicebeta} a double precision giving the Slice sampling parameter
                     for the regression coefficients, and \code{slicev} a double precision giving
                     the Slice sampling parameter for the stick-breaking parameters.}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{ngrid}{    integer giving the number of grid points where the conditional density estimate is 
                     evaluated. The default is 100.}

    \item{grid}{     vector of grid points where the conditional density estimate is 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}

    \item{compute.band}{  logical variable indicating whether the credible band for the conditional density  
                          and mean function must be computed.}       

    \item{type.band}{  string indication the type of credible band to be computed; if equal to "HPD" or
                       "PD" then the 95 percent pointwise HPD or PD band is computed, respectively.}       

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{LDBDPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This generic function fits a Linear Dependent Dirichlet-Bernstein model (Barrientos, Jara and 
  Quintana, 2010), given by:
  \deqn{y_i | G_{X_i} \sim G_{X_i}}{yi | GXi ~ GXi}
  \deqn{\{G_{X}: X \in \mathcal{X} \}| k,\alpha, G_0 \sim LDBDP(k,\alpha G_0)}{{GX: X in X} | k,alpha, G0 ~ LDBDP(k,alpha G0)}
   where, \eqn{G_0 = N(\beta| \mu_b, S_b)}{G0 = N(beta| mub, Sb)}. To
   complete the model specification, independent hyperpriors are assumed,
  \deqn{k | \lambda \sim Poisson(\lambda)I(k \geq 1)}{k | lambda ~ Poisson(lambda)I(k > 0)}
  \deqn{\mu_b | m_0, S_0 \sim N(m_0,S_0)}{mub | m0, S0 ~ N(m0,S0)}
  \deqn{S_b | \nu, \Psi \sim IW(\nu,\Psi)}{Sb | nu, psi ~ IW(nu,psi)}
  
  Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
  \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.

  Note also that the LDBDP model is a extension of the Bernstein-Dirichlet model for density
  estimation (Petrone, 1999a, 1999b; Petrone and Waserman, 2002).


  The computational implementation of the model is based on the finite approximation to the dependent 
  Dirichlet process prior and on the use of conditional MCMC methods. The regression coefficients and
  stick-breaking parameters are updated jointly using multivariate Slice sampling (Neal, 2003). The 
  degree  of the Bernstein polynomial is updated using a Metropolis-Hasting algorithm.
}

\value{
  An object of class \code{LDBDPdensity} representing the LDBDP model fit. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{k}. \code{mub} and \code{Sb}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{k}{an integer giving the degree of the Bernstein polynomial.}

  \item{beta}{a matrix of dimension \code{maxn} times the 
                  number of columns in the design matrix, giving the 
                  regression coefficients for each stick-breaking component.}
                  
  \item{alpha}{ giving the value of the precision parameter.} 

  \item{mub}{ giving the mean of the normal baseline distributions.}

  \item{Sb}{ giving the covariance matrix the normal baseline distributions.}
 
  \item{v}{  giving the \code{maxn} vector of stick-breaking beta random variables.
             The last element in this vector must be equal to 1.}
}

\seealso{
\code{\link{DPcdensity}}, \code{\link{LDDPdensity}}
}

\references{
Barrientos, F., Jara, A., Quintana, F. (2010). Bounded density regression using 
   dependent Bernstein polynomials. Technical Report, Department of Statistics, 
   Pontificia Universidad Catolica de Chile.

Neal, R. (2003) Slice sampling. Anals of Statistics, 31: 705-767. 

Petrone, S. (1999a) Random Bernstein polynomials. Scandinavian Journal of 
   Statistics, 26: 373-393.

Petrone, S. (1999b) Bayesian density estimation using Bernstein polynomials. 
   The Canadian Journal of Statistics, 27: 105-126.

Petrone, S. and Waserman, L. (2002) Consistency of Bernstein polynomial posterior. 
   Journal of the Royal Statistical Society, Series B, 64: 79-100.
}

\examples{
\dontrun{

    ######################## 
    # Simulate data
    ########################
      nrec <- 500
      x <- runif(nrec)
      y <- rep(0,nrec)
      for(i in 1:nrec)
      {
           y[i] <- ifelse(runif(1) < (0.8-0.5*x[i]^2), 
                    rbeta(1,22-(x[i]^2)*20,5+x[i]*20),
                    rbeta(1,8+x[i]*5,20))
      }

    # true model

      true.dens <- function(grid,x)
      {
	  (0.8-0.5*x^2)*dbeta(grid,22-(x^2)*20,5+x*20)+
          (0.2+0.5*x^2)*dbeta(grid,8+x*5,20)
      }	

      true.mean <- function(x)
      {
          (0.8-0.5*x^2)*(22-(x^2)*20)/(22-(x^2)*20+5+x*20)+
          (0.2+0.5*x^2)*(8+x*5)/(8+x*5+20)
      }	

    # predictions
      grid <- seq(0,1,len=100)
      npred <- 25
      xpred <- matrix(1,ncol=2,nrow=npred)
      xpred[,2] <- seq(0,1,len=npred)

    # prior
      prior <- list(maxn = 25,   
                    alpha = 1, 
                    lambda = 25, 
                    nu = 4, 
                    psiinv = diag(1000,2), 
                    m0 = rep(0,2),
                    S0 = diag(1000,2))

    # mcmc
      mcmc <- list(nburn = 5000, 
                   nskip = 3, 
                   ndisplay = 100, 
                   nsave = 5000)

    # state
      state <- NULL

    # fitting the model

      fit <- LDBDPdensity(formula=y~x,xpred=xpred,
                          prior=prior,
                          mcmc=mcmc,
                          state=NULL,status=TRUE,
                          grid=grid,
                          compute.band=TRUE,type.band="PD")

      fit
      summary(fit)
      plot(fit)

    # Plots for some predictions
    # (conditional density and mean function)

      par(mfrow=c(2,2))
      plot(fit$grid,fit$densp.h[7,],type="l",lwd=2,
           xlim=c(0,1),ylim=c(0,4),xlab="y",ylab="density",lty=2) 
      lines(fit$grid,fit$densp.m[7,],lwd=2,lty=1) 
      lines(fit$grid,fit$densp.l[7,],lwd=2,lty=2)
      lines(fit$grid,true.dens(fit$grid,fit$xpred[7,2]),lwd=2,col="red")

      plot(fit$grid,fit$densp.h[13,],type="l",lwd=2,
           xlim=c(0,1),ylim=c(0,4),xlab="y",ylab="density",lty=2) 
      lines(fit$grid,fit$densp.m[13,],lwd=2,lty=1) 
      lines(fit$grid,fit$densp.l[13,],lwd=2,lty=2)
      lines(fit$grid,true.dens(fit$grid,fit$xpred[13,2]),lwd=2,col="red")

      plot(fit$grid,fit$densp.h[19,],type="l",lwd=2,
           xlim=c(0,1),ylim=c(0,4),xlab="y",ylab="density",lty=2) 
      lines(fit$grid,fit$densp.m[19,],lwd=2,lty=1) 
      lines(fit$grid,fit$densp.l[19,],lwd=2,lty=2)
      lines(fit$grid,true.dens(fit$grid,fit$xpred[19,2]),lwd=2,col="red")

      plot(x,y) 
      lines(fit$xpred[,2],fit$meanfp.m,lwd=2,lty=1) 
      lines(fit$xpred[,2],fit$meanfp.l,lwd=2,lty=2)
      lines(fit$xpred[,2],fit$meanfp.h,lwd=2,lty=2)
      lines(fit$xpred[,2],true.mean(fit$xpred[,2]),lwd=2,lty=1,col="red")
}
}

\author{
Felipe Barrientos \email{<afbarrie@.mat.puc.cl>}

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
