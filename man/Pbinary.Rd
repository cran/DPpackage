\name{Pbinary}
\alias{Pbinary}
\alias{Pbinary.default}
\encoding{latin1}

\title{Bayesian analysis for a parametric Bernoulli regression model}
\description{
    This function generates a posterior density sample
    for a parametric binary regression model. 
}
    
\usage{

Pbinary(formula,link="logit",prior,mcmc,state,status,misc=NULL,
         data=sys.frame(sys.parent()),na.action=na.fail)
        
}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.} 
                     
   \item{link}{      a description of the link function to
                     be used in the model. The links considered by 
                     \code{Pbinary} so far are \emph{logit} (default), \emph{probit},
                     \emph{cloglog}, and \emph{cauchy}.} 

   \item{prior}{     a list giving the prior information. The list includes the following
                     parameters: \emph{beta0} and \emph{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the regression
                     coefficients.}
               
    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \emph{nburn} giving the number of burn-in 
                     scans, \emph{nskip} giving the thinning interval, \emph{nsave} giving
                     the total number of scans to be saved, \emph{ndisplay} giving
                     the number of saved scans to be displayed on the screen (the function reports 
                     on the screen when every \emph{ndisplay} iterations have been carried
                     out), and \emph{tune} giving the Metropolis tuning parameter.}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}

    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \emph{state}.}

    \item{misc}{     misclassification information. When used, this list must include
                     two objects, \emph{sens} and \emph{spec}, giving the sensitivity and
                     specificity, respectively. Both can be a vector or a scalar. 
                     This information is used to correct for misclassification in the
                     conditional bernoulli model.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{Pbinary} to print an error message and terminate if there are any
                     incomplete observations.}       

}

\details{\code{Pbinary} simulates from the posterior density of a 
  parametric Bernoulli regression model,
  
  \deqn{y_i \sim \mathcal{B}ernoulli(\pi_i)}{yi ~ Bernoulli(pii)}

  where \eqn{\pi_i = F(X_i \beta)}{pii = F(Xi beta)} and \eqn{F} is a distribution
  function on the real line known as the inverse of the link function
  in the context of generalized linear models. The links considered by 
  \code{Pbinary} so far are \emph{logit} (default), \emph{probit},
  \emph{cloglog}, and \emph{cauchy}. 
  
  To complete the model specification, the following
  prior distribution is assumed,

  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  
  A Metropolis-Hastings step is used to sample the posterior distribution
  of the regression coefficients. The Metropolis proposal distribution is 
  centered at its current value and the variance-covariance matrix 
  correspond to the variance-covariance matrix of the MLEs times the 
  tunning parameter, \emph{tune}, specified in the \emph{mcmc} list.
  
  When the model considers correction for misclassification, a modified 
  link function is used. The modified link is a function of the
  sensitivity and specificity of the classification
  (see, e.g., Jara, Garcia-Zattera and Lesaffre, 2006).
}

\value{
  An object of class \code{Pbinary} representing the parametric regression
  model fit. Generic functions such as \code{print}, \code{plot}, \code{summary}, 
  \code{predict}, and \code{anova} have methods to show the results of the fit. 
  The results include only the regression coefficients, \emph{beta}.
  
  The MCMC samples of the parameters are stored in the object \emph{thetasave}. 
  This object is included in the list \emph{save.state} and is a matrix which can 
  be analyzed directly by functions provided by the coda package.
  
  The list \emph{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \emph{status=TRUE} and create the list state based on 
  this starting values. In this case the list \emph{state} must include the following objects:
  
  \item{beta}{ giving the value of the regression coefficients.} 
}

\references{
Jara, A., Garcia-Zattera, M.J., Lesaffre, E. (2006) Semiparametric Bayesian
  Analysis of Misclassified Binary Data. XXIII International Biometric Conference,
  July 16-21, Montreal, Canada.
}

\examples{
\dontrun{
    # Bioassay Data Example
    # Cox, D.R. and Snell, E.J. (1989). Analysis of Binary Data. 2nd ed. 
    # Chapman and Hall. p. 7  
    # In this example there are 150 subjects at 5 different stimulus levels, 
    # 30 at each level.

      y<-c(rep(0,30-2),rep(1,2),
           rep(0,30-8),rep(1,8),
           rep(0,30-15),rep(1,15),
           rep(0,30-23),rep(1,23),
           rep(0,30-27),rep(1,27))

      x<-c(rep(0,30),
           rep(1,30),
           rep(2,30),
           rep(3,30),
           rep(4,30))

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn<-5000
      nsave<-5000
      nskip<-10
      ndisplay<-1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                   tune=1.1)


    # Prior distribution
      prior <- list(beta0=rep(0,2), Sbeta0=diag(10000,2))

    # Fit a logistic regression model
      fit1 <- Pbinary(y~x,link="logit",prior=prior,
                      mcmc=mcmc,state=state,status=TRUE) 
      fit1

    # Fit a probit regression model
      fit2 <- Pbinary(y~x,link="probit",prior=prior, 
                      mcmc=mcmc,state=state,status=TRUE) 
      fit2

    # Fit a cloglog regression model
      fit3 <- Pbinary(y~x,link="cloglog",prior=prior,
                      mcmc=mcmc,state=state,status=TRUE) 
      fit3

    # Fit a cauchy regression model
      fit4 <- Pbinary(y~x,link="cauchy",prior=prior,
                      mcmc=mcmc,state=state,status=TRUE) 
      fit4
}      
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
