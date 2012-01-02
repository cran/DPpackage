\name{DPcdensity}
\alias{DPcdensity}
\alias{DPcdensity.default}

\title{Bayesian Semiparametric Conditional Density Estimation using a DPM of normals}
\description{
    This function generates a posterior density sample for a 
    Bayesian density regression model with continuous predictors using
    a Dirichlet process mixture of normals model.
}
  
\usage{
DPcdensity(y,x,xpred,ngrid=100,grid=NULL,compute.band=FALSE,
           type.band="PD",prior,mcmc,state,status,
           data=sys.frame(sys.parent()),work.dir=NULL)
}

\arguments{
    \item{y}{        a vector giving the data from which the density estimate 
                     is to be computed.}

    \item{x}{        a vector or matrix giving the continuous predictors of
                     dimension \code{nrec} times \code{nx}.}

    \item{xpred}{    a vector or matrix giving the values of the continuous 
                     predictors used for prediction.}

    \item{ngrid}{    number of grid points where the conditional density estimate is 
                     evaluated.  The default is 100.}

    \item{grid}{     vector of grid points where the conditional density estimate is 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}

    \item{compute.band}{  logical variable indicating whether the credible band for the conditional  
                          density and mean function must be computed.}       

    \item{type.band}{  string indication the type of credible band to be computed; if equal to "HPD" or
                       "PD" then the 95 percent pointwise HPD or PD band is computed, respectively.}       

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
                     see details below). Notice that the dimension of the baseline measure 
                     includes the predictor and the response, i.e., \code{nx+1}.}

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
    
     \item{work.dir}{working directory.} 
}

\details{
  This generic function fits a Dirichlet process mixture of normal model 
  (Escobar and West, 1995) for the conditional density estimation 
  \eqn{f(y \mid x)}{f(y|x)} as proposed by Muller, Erkanli and West (1996). They proposed to 
  specify a Dirichlet process mixture of normals for the joint distribution of the response
  and predictors. Although in the original paper these authors focussed on the mean regression 
  function, their method can be used to model the conditional density of the response giving the
  predictors in a semiparametric way. Indeed, their method is essentially a locally weighted 
  mixture of normal regression models with weigths predictor-dependent.

  Let \eqn{y_i}{yi} and \eqn{X_i}{Xi} be the response and the vector of predictors, respectively.
  Further, let \eqn{Z_i=(y_i,X_i)}{Zi=(yi,Xi)}. The model for the joint distribution of the response and predictors is as follows:

  \deqn{Z_i | \mu_i, \Sigma_i \sim N(\mu_i,\Sigma_i), i=1,\ldots,n}{Zi | mui, Sigmai ~ N(mui,Sigmai), i=1,\ldots,n}
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
  described in Neal (2000) is adopted. Specifically, the algorithm 8 with \eqn{m=1} 
  of Neal (2000) is considered in the \code{DPcdensity} function. 
}

\value{
  An object of class \code{DPcdensity} representing the DP mixture of normals
  model fit. Generic functions such as \code{print}, \code{summary}, and \code{plot} have methods to 
  show the results of the fit. The results include the baseline parameters, \code{alpha}, and the 
  number of clusters.
  
  The MCMC samples of the parameters and the errors in the model are stored in the object 
  \code{thetasave}. The object is included in the 
  list \code{save.state} and are matrices which can be analyzed directly by functions 
  provided by the coda package.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{an integer giving the number of clusters.} 
  
  \item{muclus}{ a matrix of dimension (nobservations+2)*(nvariables) giving the means of the clusters 
  (only the first \code{ncluster} are considered to start the chain).} 
  
  \item{sigmaclus}{ a matrix of dimension (nobservations+2)*( (nvariables)*((nvariables)+1)/2) giving 
  the lower matrix of the covariance matrix of the clusters (only the first \code{ncluster} are 
  considered to start the chain).} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each observation belongs.}

  \item{alpha}{ giving the value of the precision parameter.} 
  
  \item{m1}{ giving the mean of the normal components of the baseline distribution.}
  
  \item{k0}{ giving the scale parameter of the normal part of the baseline distribution.}
  
  \item{psi1}{giving the scale matrix of the inverted-Wishart part of the baseline distribution.}

  \item{z}{giving the matrix of response and predictors. This must be included if missing data
           (response and/or predictors) are present. Those are imputed during the MCMC.}
}

\seealso{
\code{\link{DPdensity}}, \code{\link{PTdensity}}, \code{\link{BDPdensity}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Muller, P., Erkanli, A. and West, M. (1996) Bayesian curve fitting using 
  multivariate normal mixtures. Biometrika, 83: 67-79.

Neal, R. M. (2000). Markov Chain sampling methods for Dirichlet process mixture models.
  Journal of Computational and Graphical Statistics, 9: 249-265.

}

\examples{
\dontrun{

    ########################################################## 
    # Simulated data:
    # Here we replicate the results reported with 
    # simulated data by Dunson, Pillai and Park (2007,
    # JRSS Ser. B, 69: 163-183, pag 177) where a different
    # approach is proposed. 
    ##########################################################

      dtrue <- function(grid,x)
      {
          exp(-2*x)*dnorm(grid,mean=x,sd=sqrt(0.01))+
          (1-exp(-2*x))*dnorm(grid,mean=x^4,sd=sqrt(0.04))
      } 

      nobs <- 500
      x <- runif(nobs)
      y1 <- x + rnorm(nobs, 0, sqrt(0.01))
      y2 <- x^4 + rnorm(nobs, 0, sqrt(0.04))

      u <- runif(nobs)
      prob <- exp(-2*x)
      y <- ifelse(u<prob,y1,y2)

    # Prior information
      w <- cbind(y,x)  
      wbar <- apply(w,2,mean)
      wcov <- var(w)

      prior <- list(a0=10,
                    b0=1,
                    nu1=4,
                    nu2=4,
                    s2=0.5*wcov,
                    m2=wbar,
                    psiinv2=2*solve(wcov),
                    tau1=6.01,
                    tau2=2.01)

    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=5000,
                   nsave=5000,
                   nskip=3,
                   ndisplay=100)

    # fitting the model
      xpred <- c(0.00,0.05,0.10,0.15,0.20,0.25,
                 0.30,0.35,0.40,0.45,0.49,0.55,
                 0.60,0.65,0.70,0.75,0.80,0.85,
                 0.88,0.95,1.00)     

      fit <- DPcdensity(y=y,x=x,xpred=xpred,ngrid=100, 
                        prior=prior, 
                        mcmc=mcmc, 
                        state=state, 
                        status=TRUE,
                        compute.band=TRUE,type.band="PD")

    # true model and estimates
      par(mfrow=c(2,3))      

      plot(fit$grid,fit$densp.h[3,],lwd=1,type="l",lty=2,
           main="x=0.10",xlab="values",ylab="density",ylim=c(0,4))
      lines(fit$grid,fit$densp.l[3,],lwd=1,type="l",lty=2)
      lines(fit$grid,fit$densp.m[3,],lwd=2,type="l",lty=1)
      lines(fit$grid,dtrue(fit$grid,xpred[3]),lwd=2,
            type="l",lty=1,col="red")
      
      plot(fit$grid,fit$densp.h[6,],lwd=1,type="l",lty=2,
           main="x=0.25",xlab="values",ylab="density",ylim=c(0,4))
      lines(fit$grid,fit$densp.l[6,],lwd=1,type="l",lty=2)
      lines(fit$grid,fit$densp.m[6,],lwd=2,type="l",lty=1)
      lines(fit$grid,dtrue(fit$grid,xpred[6]),lwd=2,
            type="l",lty=1,col="red")
      
      plot(fit$grid,fit$densp.h[11,],lwd=1,type="l",lty=2,
           main="x=0.49",xlab="values",ylab="density",ylim=c(0,4))
      lines(fit$grid,fit$densp.l[11,],lwd=1,type="l",lty=2)
      lines(fit$grid,fit$densp.m[11,],lwd=2,type="l",lty=1)
      lines(fit$grid,dtrue(fit$grid,xpred[11]),lwd=2,type="l",
            lty=1,col="red")
      
      plot(fit$grid,fit$densp.h[16,],lwd=1,type="l",lty=2,
           main="x=0.75",xlab="values",ylab="density",ylim=c(0,4))
      lines(fit$grid,fit$densp.l[16,],lwd=1,type="l",lty=2)
      lines(fit$grid,fit$densp.m[16,],lwd=2,type="l",lty=1)
      lines(fit$grid,dtrue(fit$grid,xpred[16]),lwd=2,type="l",
            lty=1,col="red")
      
      plot(fit$grid,fit$densp.h[19,],lwd=1,type="l",lty=2,
           main="x=0.75",xlab="values",ylab="density",ylim=c(0,4))
      lines(fit$grid,fit$densp.l[19,],lwd=1,type="l",lty=2)
      lines(fit$grid,fit$densp.m[19,],lwd=2,type="l",lty=1)
      lines(fit$grid,dtrue(fit$grid,xpred[19]),lwd=2,type="l",
            lty=1,col="red")

    # data and mean function
      plot(x,y,xlab="x",ylab="y",main="")
      lines(xpred,fit$meanfp.m,type="l",lwd=2,lty=1)
      lines(xpred,fit$meanfp.l,type="l",lwd=2,lty=2)
      lines(xpred,fit$meanfp.h,type="l",lwd=2,lty=2)

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}

}

\keyword{models}
\keyword{nonparametric}
