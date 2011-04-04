\name{PTglmm}
\alias{PTglmm}
\alias{PTglmm.default}

\title{Bayesian analysis for a semiparametric generalized linear mixed model using a MMPT}
\description{
    This function generates a posterior density sample for a 
    semiparametric generalized linear mixed model, using a Mixture of
    Multivariate Polya Trees prior for the distribution of the 
    random effects. 
}
  
  
\usage{

PTglmm(fixed,random,family,offset,n,prior,mcmc,state,status,
      data=sys.frame(sys.parent()),na.action=na.fail)
}

\arguments{
    \item{fixed}{    a two-sided linear formula object describing the
                     fixed-effects part of the model, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.}
    
    \item{random}{   a one-sided formula of the form \code{~z1+...+zn | g}, with 
                     \code{z1+...+zn} specifying the model for the random effects and 
                     \code{g} the grouping variable. The random effects formula will be
                     repeated for all levels of grouping.}

    \item{family}{   a description of the error distribution and link function to
                     be used in the model. This can be a character string naming a
                     family function, a family function or the result of a call to
                     a family function. The families(links) considered by 
                     \code{PTglmm} so far are binomial(logit), binomial(probit),
                      Gamma(log), and poisson(log). The gaussian(identity) case is 
                     implemented separately in the function \code{PTlmm}.}
                     
    \item{offset}{   this can be used to specify an a priori known component to
                     be included in the linear predictor during the fitting (only for
                     poisson and gamma models).}

    \item{n}{        this can be used to indicate the total number of cases in a binomial
                     model (only implemented for the logistic link). If it is not specified
                     the response variable must be binary.}

    \item{prior}{    a list giving the prior information. The list include the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Polya Tree (PT)
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{nu0} and \code{tinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix of the normal
                     baseline distribution, \code{sigma} giving the value of the covariance
                     matrix of the centering distribution (it must be specified if 
                     \code{nu0} and \code{tinv} are missing),
                     \code{mub} and \code{Sb} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution, \code{mu} giving the value of the mean of the
                     centering distribution (it must be specified if 
                     \code{mub} and \code{Sb} are missing), \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the fixed effects
                     (must be specified only if fixed effects are considered in the model), 
                     \code{M} giving the finite level
                     of the PT prior to be considered, and \code{frstlprob} a logical variable
                     indicating whether the first level probabilities of the PT are fixed
                     or not (the default is FALSE), \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution for the inverse of the precision parameter of the 
                     Gamma model (they must be specified only if the Gamma model is considered), and
                     \code{typep} indicating whether the type of decomposition of the centering 
                     covariance matrix is random (1) or not (0).} 

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), \code{nbase} giving the number scans to be performed before the 
                     parameters of the centering distribution and the precision parameter are
                     updated (i.e., the update of this parameters is invoked only once in every 
                     \code{nbase} scans) (the default value is 1), \code{tune1}, \code{tune2}, \code{tune3}, 
                     \code{tune4} and \code{tune5} giving the Metropolis tuning parameter for the baseline mean, 
                     variance, precision parameter, partition and dispersion parameter (only for the Gamma mode),
                     respectively. If \code{tune1}, \code{tune2}, \code{tune3} or \code{tune4} are not 
                     specified or negative, an adpative Metropolis algorithm is performed. 
                     If \code{tune5} is not specified, a default value 
                     of 1.1 is assumed. Finally, the integer \code{samplef} indicates whether
                     the functional parameters must be sample (1) or not (0).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{PTglmm} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a generalized linear mixed-effects model using 
  a Mixture of Multivariate Polya Trees prior (see, Lavine 1992; 1994, for details about univariate PT) 
  for the distribution of the random effects as described in Jara, Hanson and Lesaffre (2009). The linear
  predictor is modeled as follows:
  \deqn{\eta_i = X_i \beta_F + Z_i \beta_R + Z_i b_i, i=1,\ldots,n}{etai = Xi betaF + Zi betaR + Zi bi, i=1,\ldots,n}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha,\mu,\Sigma,O \sim PT^M(\Pi^{\mu,\Sigma,O},\mathcal{A})}{G | alpha,mu,Sigma,O ~ PT^M(Pi^{mu,Sigma,O},A)}

  where, \eqn{\theta_i = \beta_R + b_i}{thetai = betaR + bi}, \eqn{\beta = \beta_F}{beta = betaF}, and \eqn{O}{O}
  is an orthogonal matrix defining the decomposition of the centering covariance matrix. As in Hanson (2006), the PT prior is centered around the 
  \eqn{N_d(\mu,\Sigma)}{N_d(mu,Sigma)} distribution. However, we consider the class of partitions \eqn{\Pi^{\mu,\Sigma, O}}{Pi^{mu,Sigma, O}}. The
  partitions starts with base sets that are Cartesian products of intervals obtained as quantiles from the standard normal distribution. A multivariate 
  location-scale transformation \eqn{\theta=\mu+\Sigma^{1/2} z}{theta=mu+Sigma^{1/2} z} is applied to each 
  base set yielding the final sets. Here \eqn{\Sigma^{1/2}=T'O'}{Sigma^{1/2}=T'O'}, where \eqn{T}{T}
  is the unique upper triangular Cholesky matrix of \eqn{\Sigma}{Sigma}. The family \eqn{\mathcal{A}=\{\alpha_e: e \in E^{*}\}}{A=\{alphae: e \in E*\}}, 
  where \eqn{E^{*}=\bigcup_{m=0}^{M} E_d^m}{E*=U_{m=0}^{M} E_d^m}, 
  with \eqn{E_d}{E_d} and \eqn{E_d^m}{E_m} the \eqn{d}{d}-fold product of \eqn{E=\{0,1\}}
  and the the \eqn{m}-fold product of \eqn{E_d}{E_d}, respectively. The family \eqn{\mathcal{A}}{A}
  was specified 
  as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}. 

  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\Sigma | \nu_0, T \sim IW(\nu_0,T)}{Sigma | nu0, T ~ IW(nu0,T)}
  \deqn{O \sim Haar(q)}{O ~ Haar(q)}

  Note that the inverted-Wishart prior is parametrized such that
  \eqn{E(\Sigma)= T^{-1}/(\nu_0-q-1)}{E(Sigma)= Tinv/(nu0-q-1)}.

  The precision parameter, \eqn{\alpha}{alpha}, of the \code{PT} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. 

  The inverse of the dispersion parameter of the Gamma model
  is modeled using \code{gamma} distribution, \eqn{\Gamma(\tau_1/2,\tau_2/2)}{Gamma(tau1/2,tau2/2)}. 

  The computational implementation of the model is based on the marginalization of
  the \code{PT} as discussed in Jara, Hanson and Lesaffre (2009). 
}

\value{
  An object of class \code{PTglmm} representing the generalized linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{betaR}, \code{betaF}, \code{mu}, the elements of \code{Sigma}, the precision parameter 
  \code{alpha}, the dispersion parameter of the Gamma model, and \code{ortho}.
  
  The function \code{PTrandom} can be used to extract the posterior mean of the 
  random effects.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{alpha}{ giving the value of the precision parameter.} 
  
  \item{b}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the random effects
  for each subject.} 

  \item{beta}{ giving the value of the fixed effects.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.}
  
  \item{sigma}{ giving the variance matrix of the normal baseline distributions.}
  
  \item{phi}{ giving the precision parameter for the Gamma model (if needed).}  

  \item{ortho}{ giving the orthogonal matrix \code{H}, used in the decomposition of the covariance matrix.} 
}

\seealso{
\code{\link{PTrandom}},
\code{\link{PTglmm}}  , \code{\link{PTolmm}}, 
\code{\link{DPMglmm}}, \code{\link{DPMlmm}}, \code{\link{DPMolmm}},
\code{\link{DPlmm}}  , \code{\link{DPglmm}}, \code{\link{DPolmm}}
}

\references{

Hanson, T. (2006) Inference for Mixtures of Finite Polya Trees.
  Journal of the American Statistical Association, 101: 1548-1565.

Jara, A., Hanson, T., Lesaffre, E. (2009) Robustifying Generalized Linear 
   Mixed Models using a New Class of Mixtures of Multivariate Polya Trees. 
   Journal of Computational and Graphical Statistics, 18(4): 838-860.

Lavine, M. (1992) Some aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 20: 1222-11235.

Lavine, M. (1994) More aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 22: 1161-1176.
}

\examples{
\dontrun{
    # Respiratory Data Example
      data(indon)
      attach(indon)

      baseage2 <- baseage**2
      follow <- age-baseage
      follow2 <- follow**2 

    # Prior information

      prior <- list(alpha=1,
                    M=4, 
                    frstlprob=FALSE,
                    nu0=4,
                    tinv=diag(1,1),
                    mub=rep(0,1),
                    Sb=diag(1000,1),
                    beta0=rep(0,9),
                    Sbeta0=diag(10000,9))

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 5000
      nskip <- 20
      ndisplay <- 100
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay,
                   tune1=0.5,tune2=0.5,
                   samplef=1)

    # Fitting the Logit model
      fit1 <- PTglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+baseage2+
                     follow+follow2,random=~1|id,family=binomial(logit),
                     prior=prior,mcmc=mcmc,state=state,status=TRUE)

      fit1 

      plot(PTrandom(fit1,predictive=TRUE))

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

    # Extract random effects
      PTrandom(fit1)
      PTrandom(fit1,centered=TRUE)

    # Extract predictive information of random effects
      PTrandom(fit1,predictive=TRUE)

    # Predictive marginal and joint distributions      
      plot(PTrandom(fit1,predictive=TRUE))

    # Fitting the Probit model
      fit2 <- PTglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+baseage2+
                     follow+follow2,random=~1|id,family=binomial(probit),
                     prior=prior,mcmc=mcmc,state=state,status=TRUE)
      fit2 

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit2,ask=FALSE)
      plot(fit2,ask=FALSE,nfigr=2,nfigc=2)	

    # Extract random effects
      PTrandom(fit2)

    # Extract predictive information of random effects
      PTrandom(fit2,predictive=TRUE)

    # Predictive marginal and joint distributions      
      plot(PTrandom(fit2,predictive=TRUE))
}
}
\author{
Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{models}
\keyword{nonparametric}
