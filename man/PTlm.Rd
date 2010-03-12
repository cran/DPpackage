\name{PTlm}
\alias{PTlm}
\alias{PTlm.default}

\title{Bayesian analysis for a semiparametric linear regression model}
\description{
    This function generates a posterior density sample
    from a semiparametric linear regression model using a 
    Mixture of Polya Trees prior for the distribution of the errors.
}
    
\usage{

PTlm(formula,ngrid=200,grid=NULL,prior,mcmc,state,status,
     data=sys.frame(sys.parent()),na.action=na.fail)

}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.}
                     
   \item{ngrid}{     number of grid points where the error density estimate is 
                     evaluated. The default value is 200.}

    \item{grid}{     grid points where the density estimate is 
                     evaluated. The default is NULL.}

   \item{prior}{     a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Polya Tree
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{tau1} and \code{tau2} giving the 
                     hyperparameters for the prior distribution of the variance
                     of the normal baseline distribution, \code{beta0} and \code{Sbeta0} 
                     giving the hyperparameters of the normal prior distribution for the regression
                     coefficients, optionally \code{M} giving the finite level
                     to be considered, and \code{frstlprob} a logical variable
                     indicating whether the first level probabilities of the PT are fixed
                     defining a median regression model (the default is TRUE). Note that 
                     if \code{M} is specified, a Partially Specified
                     Mixture of Polya trees is fitted.}
               
    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, and \code{ndisplay} giving
                     the number of saved scans to be displayed on the screen (the function reports 
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
                     \code{PTlm} to print an error message and terminate if there are any
                     incomplete observations.}       

}

\details{
  By default, this generic function fits a median regression model using a Scale Mixture of Polya Trees 
  prior for the distribution of the errors (see, e.g., Lavine, 1992 and 1994, Hanson and Johnson, 2004):

  \deqn{y_i = X_i \beta + V_i, i=1,\ldots,n}{yi = Xi beta + Vi, i=1,\ldots,n}
  \deqn{V_i | G \sim G}{Vi | G ~ G}
  \deqn{G | \alpha,\sigma^2 \sim PT(\Pi^{\sigma^2},\textit{A})}{G | alpha,sigma2 ~ PT(Pi^{sigma2},\textit{A})}
  
  where, the PT is centered around a \eqn{N(0,\sigma^2)}{N(0,sigma2)} distribution, by
  taking each \eqn{m} level of the partition \eqn{\Pi^{\sigma^2}}{Pi^{sigma2}} to coincide 
  with the \eqn{k/2^m, k=0,\ldots,2^m} quantile of the \eqn{N(0,\sigma^2)}{N(0,sigma2)} distribution.
  The family \eqn{\textit{A}=\{\alpha_e: e \in E^{*}\}}{\textit{A}=\{alphae: e \in E^{*}\}}, where \eqn{E^{*}=\bigcup_{m=1}^{\infty} E^m}{E^{*}=\bigcup_{m=1}^{+infty} E^m} 
  and \eqn{E^m} is the \eqn{m}-fold product of \eqn{E=\{0,1\}},
  was specified as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=alpha m^2}. 
  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}

  Optionally, if frstlprob=FALSE (the default value is TRUE) is specified, a mean regression model is considered. In this case, the 
  following PT prior is considered:
  \deqn{G | \alpha,\mu, \sigma^2 \sim PT(\Pi^{\mu, \sigma^2},\textit{A})}{G | alpha,mu,sigma2 ~ PT(Pi^{mu,sigma2},\textit{A})}
  where, the PT is centered around a \eqn{N(0,\mu,\sigma^2)}{N(mu,sigma2)} distribution. In this case,
  the intercept term is automatically excluded from the model and the hyperparameters for the
  normal prior for \eqn{\mu}{mu} must be specified. The normal prior is given by,

  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}  
  
  The precision parameter, \eqn{\alpha}{alpha}, of the \code{PT} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  In the computational implementation of the model, Metropolis-Hastings steps are used to 
  sample the posterior distribution of the regression coefficients and 
  hyperparameters.
}


\value{
  An object of class \code{PTlm} representing the semiparametric median regression
  model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{beta}, \code{mu}, \code{sigma2}, and the precision
  parameter \code{alpha}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects:

  \item{alpha}{ giving the value of the precision parameter.}
  
  \item{beta}{ giving the value of the regression coefficients.} 

  \item{mu}{ giving the mean of the normal baseline distribution (If needed).}
  
  \item{sigma2}{ giving the variance of the normal baseline distribution.}
  
  \item{v}{ giving the value of the errors (it must be consistent with the data.} 
}


\references{

Hanson, T., and Johnson, W. (2002) Modeling regression error with a Mixture
  of Polya Trees. Journal of the American Statistical Association, 
  97: 1020 - 1033.

Lavine, M. (1992) Some aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 20: 1222-11235.

Lavine, M. (1994) More aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 22: 1161-1176.

}

\examples{
\dontrun{
    ####################################
    # A simulated Data Set
    # (Mixture of Normals)
    ####################################

      ind<-rbinom(100,1,0.5)
      vsim<-ind*rnorm(100,1,0.15)+(1-ind)*rnorm(100,3,0.15)

      x1<-rep(c(0,1),50)
      x2<-rnorm(100,0,1)

      etasim<-x1+-1*x2
      y<-etasim+vsim


    # Initial state
      state <- NULL

    # MCMC parameters
      nburn<-5000
      nsave<-10000
      nskip<-20
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay)

    # Prior information
      prior <- list(alpha=1,beta0=rep(0,3),Sbeta0=diag(1000,3),
                    tau1=0.01,tau2=0.01,M=6)

    # Fit the model

      fit1 <- PTlm(formula=y~x1+x2,prior=prior,mcmc=mcmc,state=state,
                   status=TRUE) 

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1)
      plot(fit1,nfigr=2,nfigc=2)

    # Table of Pseudo Contour Probabilities
      anova(fit1)
      

    ############################################
    # The Australian Institute of Sport's data
    # (Skew data example)
    ############################################
      data(sports)
      attach(sports)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-5000
      nsave<-10000
      nskip<-20
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay)

    # Prior information
      prior <- list(alpha=1,beta0=rep(0,3),Sbeta0=diag(1000,3),
                    tau1=0.01,tau2=0.01,M=8)

    # Fit the model

      fit2 <- PTlm(formula=bmi~lbm+gender,prior=prior,mcmc=mcmc,
                   state=state,status=TRUE) 

    # Summary with HPD and Credibility intervals
      summary(fit2)
      summary(fit2,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit2)
      plot(fit2,nfigr=2,nfigc=2)

    # Table of Pseudo Contour Probabilities
      anova(fit2)

}      
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}