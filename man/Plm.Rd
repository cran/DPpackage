\name{Plm}
\alias{Plm}
\alias{Plm.default}

\title{Bayesian analysis for a parametric linear regression model}
\description{
    This function generates a posterior density sample
    from a parametric linear regression model using a 
    normal distribution of the errors.
}
    
\usage{

Plm(formula,prior,mcmc,state,status,
     data=sys.frame(sys.parent()),na.action=na.fail)

}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.}
             
   \item{prior}{     a list giving the prior information. The list includes the following
                     parameter: \code{tau1} and \code{tau2} giving the 
                     hyperparameters for the prior distribution of the error variance,
                     \code{beta0} and \code{Sbeta0} 
                     giving the hyperparameters of the normal prior distribution for the regression
                     coefficients.}
               
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
                     \code{Plm} to print an error message and terminate if there are any
                     incomplete observations.}       

}

\details{
  This generic function fits a linear regression model:

  \deqn{y_i = X_i \beta + V_i, i=1,\ldots,n}{yi = Xi beta + Vi, i=1,\ldots,n}
  \deqn{V_i | \sigma^2 \sim N(0,\sigma^2)}{Vi | sigma2 ~ N(0,sigma2)}

  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0, S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}  
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}
}

\value{
  An object of class \code{Plm} representing the parametric linear regression
  model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{beta}, and \code{sigma2}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects:

  \item{beta}{ giving the value of the regression coefficients.} 

  \item{sigma2}{ giving the error variance.}
}

\examples{
\dontrun{

    ############################################
    # The Australian Institute of Sport's data
    ############################################
      data(sports)
      attach(sports)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 10000
      nskip <- 20
      ndisplay <- 100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                   ndisplay=ndisplay)

    # Prior information
      prior <- list(beta0=rep(0,3),
                    Sbeta0=diag(1000,3),
                    tau1=0.01,
                    tau2=0.01)

    # Fit the model

      fit <- Plm(formula=bmi~lbm+gender,prior=prior,mcmc=mcmc,
                  state=state,status=TRUE) 

    # Summary with HPD and Credibility intervals
      summary(fit)
      summary(fit,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit)
      plot(fit,nfigr=2,nfigc=2)

    # Table of Pseudo Contour Probabilities
      anova(fit)

}      
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}