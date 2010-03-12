\name{DPsurvint}
\alias{DPsurvint}
\alias{DPsurvint.default}

\title{Bayesian analysis for a semiparametric AFT regression model}
\description{
    This function generates a posterior density sample
    from a semiparametric AFT regression model for interval-censored 
    data.
}
    
\usage{

DPsurvint(formula,prior,mcmc,state,status,
          data=sys.frame(sys.parent()),na.action=na.fail)
}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. In the response matrix, the unknown limits
                     should be -999.}

   \item{prior}{     a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{m0} and \code{s0} giving the mean and variance of the 
                     normal prior distribution for the mean of the log normal
                     baseline distribution, and, \code{tau1} and \code{tau2} giving the 
                     hyperparameters for the prior distribution of the variance
                     of the log normal baseline distribution, and \code{beta0} and \code{Sbeta0} 
                     giving the hyperparameters of the normal prior distribution for the regression
                     coefficients.}
               
    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on the screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), and \code{tune} giving the Metropolis tuning parameter.}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}

    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}
                     
    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{DPsurvint} to print an error message and terminate if there are any
                     incomplete observations.}       

}

\details{
  This generic function fits a Mixture of Dirichlet process in a AFT regression model
  for interval censored data (Hanson and Johnson, 2004):
  \deqn{T_i = exp(- X_i \beta) V_i, i=1,\ldots,n}{Ti = exp(- Xi beta) Vi, i=1,\ldots,n}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{\beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{V_i | G \sim G}{Vi | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
  
  where, \eqn{G_0 = Log Normal(V| \mu, \sigma)}{G0 = Log Normal(V| mu, sigma)}. To complete the model specification, independent 
  hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu | m_0, s_0 \sim N(m_0,s_0)}{mu | m0, s0 ~ N(m0,s0)}
  \deqn{\sigma^{-1} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-1 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}
  
  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  In the computational implementation of the model, \eqn{G} is considered as latent data
  and sampled partially with sufficient accuracy to be able to generate \eqn{V_1,\ldots,V_{n+1}}{V1,\ldots,Vn+1}
  which are exactly iid \eqn{G}, as proposed by Doss (1994). Both Ferguson's definition
  of DP and the Sethuraman-Tiwari (1982) representation of the process are used,
  as described in Hanson and Johnson (2004) to allow the inclusion of covariates.
  
  A Metropolis-Hastings step is used to sample the fully conditional distribution
  of the regression coefficients and errors (see, Hanson and Johnson, 2004). An extra 
  step which moves the clusters in such a way that the posterior distribution
  is still a stationary distribution, is performed in order to improve the rate of
  mixing.

}


\value{
  An object of class \code{DPsurvint} representing the semiparametric AFT regression
  model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{beta}, \code{mu}, \code{sigma}, the precision
  parameter \code{alpha}, and the number of clusters.
  
  The function \code{predict.DPsurvint} can be used to extract posterior 
  information of the survival curve.

  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects:
  
  \item{beta}{ giving the value of the regression coefficients.} 
  
  \item{v}{ giving the value of the errors (it must be consistent with the data.} 
  
  \item{mu}{ giving the mean of the lognormal baseline distribution.}
  
  \item{sigma}{ giving the variance of the lognormal baseline distribution.}
  
  \item{alpha}{ giving the value of the precision parameter.}
}

\seealso{
\code{\link{predict.DPsurvint}}
}

\references{
Doss, H. (1994). Bayesian nonparametric estimation for incomplete data using 
  mixtures of Dirichlet priors. The Annals of Statistics, 22: 1763 - 1786.

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Hanson, T., and Johnson, W. (2004) A Bayesian Semiparametric AFT Model for
  Interval-Censored Data. Journal of Computational and Graphical Statistics, 
  13: 341-361.
  
Sethuraman, J., and Tiwari, R. C. (1982) Convergence of Dirichlet Measures and 
  the Interpretation of their Parameter, in Statistical Decision Theory and Related 
  Topics III (vol. 2), eds. S. S. Gupta and J. O. Berger, New York: Academic Press, 
  pp. 305 - 315.  
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

      nburn<-20000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay,tune=0.125)

    # Prior information
      prior <- list(alpha=1,beta0=rep(0,2),Sbeta0=diag(1000,2),
                    m0=0,s0=1,tau1=0.01,tau2=0.01)


    # Fit the model

      fit1 <- DPsurvint(y~x1+x2,prior=prior,mcmc=mcmc,
                        state=state,status=TRUE) 
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

    # Plot an specific model parameter 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="x1")	
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="mu")	

    # Table of Pseudo Contour Probabilities
      anova(fit1)      

    # Predictive information with covariates
      npred<-10
      xnew<-cbind(rep(1,npred),seq(-1.5,1.5,length=npred))
      xnew<-rbind(xnew,cbind(rep(0,npred),seq(-1.5,1.5,length=npred)))
      grid<-seq(0.00001,14,0.5)
      pred1<-predict(fit1,xnew=xnew,grid=grid)

    # Plot Baseline information
      plot(pred1,all=FALSE,band=TRUE)


    #############################################################
    # Time to Cosmetic Deterioration of Breast Cancer Patients
    #############################################################

      data(deterioration)
      attach(deterioration)
      y<-cbind(left,right)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-20000
      nsave<-10000
      nskip<-20
      ndisplay<-1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay,tune=0.25)

    # Prior information
      prior <- list(alpha=10,beta0=rep(0,1),Sbeta0=diag(100,1),
                    m0=0,s0=1,tau1=0.01,tau2=0.01)

    # Fitting the model

      fit2 <- DPsurvint(y~trt,prior=prior,mcmc=mcmc,
                        state=state,status=TRUE) 
      fit2
      
    # Summary with HPD and Credibility intervals
      summary(fit2)
      summary(fit2,hpd=FALSE)

    # Plot model parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit2)

    # Table of Pseudo Contour Probabilities
      anova(fit2)      

    # Predictive information with covariates
      xnew<-matrix(c(0,1),nrow=2,ncol=1)
      grid<-seq(0.01,70,1)
      pred2<-predict(fit2,xnew=xnew,grid=grid)
      plot(pred2,all=FALSE,band=TRUE)

}      
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}