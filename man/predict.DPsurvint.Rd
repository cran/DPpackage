\name{predict.DPsurvint}
\alias{predict.DPsurvint}
\title{Computes the Survival Curve in a Bayesian analysis for a semiparametric AFT regression model}
\description{
    This function generates a posterior density sample
    of the Survival curve from a semiparametric AFT regression model for interval-censored 
    data.
}
    
\usage{

\method{predict}{DPsurvint}(object,grid,xnew=NULL,hpd=TRUE, ...)

}

\arguments{

    \item{object}{   \code{DPsurvint} fitted model object.}

    \item{grid}{     a vector of grid points where the survival curve should be evaluated.}
               
    \item{xnew}{     an optional matrix containing the value of the covariables 
                     with which to predict.  If omitted, the baseline survival
                     information is calculated.}   

    \item{hpd}{      a logical variable indicating whether a 95HPD interval is
                     calculated, \code{TRUE}, or a  95Credibility interval is
                     caculated, \code{FALSE}, for the survival curve at each grid point.
                     The default value is \code{TRUE}.}
    \item{...}{      further arguments to be passed.}     
}

\details{
  This function computes the survival curve based on the fit of a Mixture of Dirichlet 
  process in a AFT regression model for interval censored data (Hanson and Johnson, 2004).
  
  Given a MCMC sample of size \eqn{J}{J} of the parameters, a sample of the 
  predictive survival curve for \eqn{X}{X} is drawn as follows: for the MCMC
  scan \eqn{j}{j} of the posterior distribution, with \eqn{j=1,\ldots,J}{j=1,...,J},
  we sample from:

  \deqn{ S^{(j)}(t|X,data) \sim Beta(a^{(j)}(t),b^{(j)}(t))}{S^{(j)}(t|X,data) ~ Beta(a^{(j)}(t),b^{(j)}(t))}
  
  where,

  \deqn{a^{(j)}(t)=\alpha^{(j)} G_0^{(j)}( (t \exp(X \beta^{(j)}) , \infty) ) + \sum_{i=1}^n \delta_{V^{(j)}_i} ( (t \exp(X \beta^{(j)}) , \infty) )}{a^{(j)}(t)=\alpha^{(j)} G_0^{(j)}( (t \exp(X \beta^{(j)}) , \infty) ) + \sum_{i=1}^n \delta_{V^{(j)}_i} ( (t \exp(X \beta^{(j)}) , \infty) )}
  
  and
  
  \deqn{b^{(j)}(t)=\alpha^{(j)}+N - a^{(j)}(t)}{b^{(j)}(t)=\alpha^{(j)}+N - a^{(j)}(t)}
}

\value{
  An object of class \code{predict.DPsurvint} representing the survival information
  arising from a \code{DPsurvint} model fit. The results include the posterior mean (\code{pmean}), 
  the posterior median (\code{pmedian}), the posterior standard deviation (\code{psd}),
  the naive standard error (\code{pstd}) and the limits of the HPD or credibility intervals,
  \code{plinf} and \code{plsup}.
}



\seealso{
\code{\link{DPsurvint}}
}

\references{
Doss, H. (1994). Bayesian nonparametric estimation for incomplete data using 
  mixtures of Dirichlet priors. The Annals of Statistics, 22: 1763 - 1786.

Hanson, T., and Johnson, W. (2004) A Bayesian Semiparametric AFT Model for
  Interval-Censored Data. Journal of Computational and Graphical Statistics, 
  13: 341-361.
}

\examples{
\dontrun{
    ####################################
    # A simulated Data Set
    ####################################

     ind<-rbinom(100,1,0.5)
     vsim<-ind*rnorm(100,1,0.25)+(1-ind)*rnorm(100,3,0.25)

     x1<-rep(c(0,1),50)
     x2<-rnorm(100,0,1)

     etasim<-x1+-1*x2
     time<-vsim*exp(-etasim)

     y<-matrix(-999,nrow=100,ncol=2)

     for(i in 1:100){
        for(j in 1:15){
         if((j-1)<time[i] & time[i]<=j){
            y[i,1]<-j-1
            y[i,2]<-j
         }
     }
     if(time[i]>15)y[i,1]<-15
     }

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-5000
      nsave<-10000
      nskip<-10
      ndisplay<-50
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay,tune=0.125)

    # Prior information
      prior <- list(alpha=10,beta0=rep(0,2),Sbeta0=diag(100000,2),
                    m0=0,s0=1,tau1=0.01,tau2=0.01)


    # Fit the model

      fit1 <- DPsurvint(y~x1+x2,prior=prior,mcmc=mcmc,
                        state=state,status=TRUE) 
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)


    # Plot model parameters
      plot(fit1)
      plot(fit1,nfigr=2,nfigc=2)	

    # Plot an specific model parameter
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="x1")	
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="mu")	


    # Predictive information for baseline survival
      grid<-seq(0.00001,14,0.5)
      pred<-predict(fit1,grid=grid)
      
    # Plot Baseline information with and without Credibility band
      plot(pred)
      plot(pred,band=TRUE)
  
    # Predictive information with covariates
      npred<-10
      xnew<-cbind(rep(1,npred),seq(-1.5,1.5,length=npred))
      xnew<-rbind(xnew,cbind(rep(0,npred),seq(-1.5,1.5,length=npred)))
      grid<-seq(0.00001,14,0.5)
      pred<-predict(fit1,xnew=xnew,grid=grid)

    # Plot Baseline information
      plot(pred,band=TRUE)
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
