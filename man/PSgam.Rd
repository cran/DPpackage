\name{PSgam}
\alias{PSgam}
\alias{PSgam.default}

\title{Bayesian analysis for a semiparametric generalized additive model}
\description{
    This function generates a posterior density sample for a 
    semiparametric generalized additive model, using a B-Splines
    and penalties. 
}
  
 
\usage{

PSgam(formula,family=gaussian(),offset,n,prior,mcmc,state,status,
       ngrid=20,data=sys.frame(sys.parent()),na.action=na.fail)
}

\arguments{
    \item{formula}{  a two-sided linear formula object describing the
                     linear predictor of the model, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. Built-in nonparametric
                     smoothing terms are indicated by \code{ps} for smoothing splines terms.  
                     See the documentation for \code{ps} for its arguments.}
    
    \item{family}{   a description of the error distribution and link function to
                     be used in the model. This can be a character string naming a
                     family function, a family function or the result of a call to
                     a family function. The families(links) considered by 
                     \code{PSgam} so far are gaussian(identity), binomial(logit), binomial(probit),
                      Gamma(log), and poisson(log).}
                     
    \item{offset}{   this can be used to specify an a priori known component to
                     be included in the linear predictor during the fitting (only for
                     poisson and gamma models).}

    \item{n}{        this can be used to indicate the total number of cases in a binomial
                     model (only implemented for the logistic link). If it is not specified
                     the response variable must be binary.}

    \item{prior}{    a list giving the prior information. The list include the following
                     parameters: \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the parametric part of 
                     the model, \code{taub1} and \code{taub2} giving the hyperparameters for the 
                     prior distribution of the inverse of the 
                     penalty parameters (the same for all), and 
                     \code{tau1} and \code{tau2} giving the hyperparameters for the 
		     prior distribution of the inverse of the dispersion parameter (only gaussian and
		     Gamma models).} 

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

    \item{ngrid}{    number of grid points where the smoothers are evaluated. The default value is 20.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{PSgam} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a generalized additive model (see, e.g., Hastie and Tibshirani, 1990)
  using Penalized splines (see, e.g., Eilers and Marx, 1996; Lang and Brezger, 2004). The linear
  predictor is modeled as follows:
  \deqn{\eta_i = X_i \beta + f_1(x_{1i})+...+f_p(x_{pi}), i=1,\ldots,n}{etai = Xi beta + f1(x1i)+...+fp(xpi), i=1,\ldots,n}
  where the effect \eqn{f} of the a covariate \eqn{x} is approximated by a polinomial
  spline with equally spaced knots, written in terms of a linear combination of B-spline basis functions.
  Specifically, the function \eqn{f} is aproximated by a spline of degree \eqn{l} with \eqn{r}
  equally spaced knots within the domain of \eqn{x}. It is well known that
  this spline can be written in terms of a linear combination of \eqn{q=l+r} B-spline
  basis,
  \deqn{f(x)=\sum_{j=1}^q b_j B_j(x).}{f(x)=sum bj Bj(x).}
 
  The computational implementation of the model is model-specific. 
   
  For the \code{poisson}, \code{Gamma}, and \code{binomial(logit)}, the full
  conditional distributions for fixed and random effects are generated through the Metropolis-Hastings algorithm 
  with a IWLS proposal (see, West, 1985 and Gamerman, 1997). 

  For the \code{binomial(probit)} case the following latent variable representation is used:
  \deqn{y_{ij} = I(w_{ij} > 0), j=1,\ldots,n_i.}{yij = I(wij > 0), j=1,\ldots,ni.}
}

\value{
  An object of class \code{PSgam} representing the generalized additive model
  fit. Generic functions such as \code{anova}, \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  the parametric component of the linear predictor \code{beta}, 
  the dispersion parameter of the Gamma or gaussian model, and
  the penalty parameters \code{sigmab}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{b}{ a vector of dimension q giving the value of the B-spline coefficients.} 

  \item{beta}{ giving the value of the parametric components of the linear predictor.}
  
  \item{sigmab}{ giving the penalty parameters.}
  
  \item{phi}{ giving the dispersion parameter for the Gamma or gaussian model (if needed).}  
}


\references{
Eilers, P.H.C. and Marx, B.D. (1996) Flexible Smoothing with B-splines and 
  penalties. Statistical Science, 11(2): 89-121.

Gamerman, D. (1997) Sampling from the posterior distribution in generalized linear
  mixed models. Statistics and Computing, 7: 57-68. 

Hastie, T. and Tibshirani, R. (1990) Generalized Additive Models.
  London: Chapman and Hall.

Lang, S., Brezger, A. (2004) Bayesian P-Splines
  Journal of Computational and Graphical Statistics, 13: 183-212. 

West, M. (1985) Generalized linear models: outlier accomodation, scale parameter and
  prior distributions. In Bayesian Statistics 2 (eds Bernardo et al.), 531-558, Amsterdam:
  North Holland.
}

\examples{
\dontrun{

 # Normal simulated data
   set.seed(0)
   n <- 400
   sig <- 2
   x0 <- runif(n, 0, 1)
   x1 <- runif(n, 0, 1)
   x2 <- runif(n, 0, 1)
   x3 <- runif(n, 0, 1)
   f0 <- function(x) 2 * sin(pi * x)
   f1 <- function(x) exp(2 * x)
   f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
   f3 <- function(x) 0*x
   f <- f0(x0) + f1(x1) + f2(x2)
   e <- rnorm(n, 0, sig)
   y <- f + e

 # prior
   prior <- list(taub1=2.02,
                 taub2=0.02,
                 beta0=rep(0,1),
                 Sbeta0=diag(100,1),
                 tau1=6.01,
                 tau2=2.01)

  # Initial state
    state <- NULL

  # MCMC parameters
    nburn <- 5000
    nsave <- 5000
    nskip <- 0
    ndisplay <- 100
    mcmc <- list(nburn=nburn,
                 nsave=nsave,
                 nskip=nskip,
                 ndisplay=ndisplay)


  # fitting the model
    fit1 <- PSgam(formula=y~ps(x0,x1,x2,x3,k=20,degree=3,pord=1),
                  family=gaussian(),prior=prior,
                  mcmc=mcmc,ngrid=30,
                  state=state,status=TRUE)


  # A binary example 
    g <- (f-5)/3
    g <- binomial()$linkinv(g)
    y <- rbinom(n,1,g)

  # fitting the model
    fit2 <- PSgam(formula=y~ps(x0,x1,x2,x3,k=20,degree=3,pord=1),
                  family=binomial(logit),prior=prior,
                  mcmc=mcmc,ngrid=30,
                  state=state,status=TRUE)

  # Poisson data
    g <- exp(f/4)
    y <- rpois(n,g)

  # fitting the model
    fit3 <- PSgam(formula=y~ps(x0,x1,x2,x3,k=20,degree=3,pord=1),
                  family=poisson(log),prior=prior,
                  mcmc=mcmc,ngrid=30,
                  state=state,status=TRUE)
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}
